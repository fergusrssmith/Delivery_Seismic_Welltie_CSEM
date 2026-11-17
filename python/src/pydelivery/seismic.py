"""Convolutional seismic forward model and the seismic likelihood set-up.

The synthetic trace is ``wavelet * reflectivity``.  Each reflection spike is
spread over the four nearest samples with a Lagrange interpolation, which
keeps the synthetic smooth in the layer times.  Because convolution is
linear, the synthetic at the error-sampling points is ``G @ R`` for a matrix
``G`` built once per trace; no FFT is needed inside the sampler.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from .config import ModelConfig, Stack
from .su import Trace, read_su

SAMPLE_TIME_MAGIC_CONSTANT = 0.253  # error-sampling interval = this / f_peak
MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER = 1.124


@dataclass
class Wavelet:
    """``samples[i]`` is the wavelet at time ``(i - n_neg) * dt`` seconds."""

    samples: np.ndarray
    dt: float
    n_neg: int  # number of samples at negative time

    @property
    def reach_before(self) -> float:
        """How far before a reflection time the wavelet extends (seconds)."""
        return self.n_neg * self.dt

    @property
    def reach_after(self) -> float:
        """How far after a reflection time the wavelet extends (seconds)."""
        return (len(self.samples) - 1 - self.n_neg) * self.dt

    def peak_frequency(self) -> float:
        """Peak frequency (Hz) from the first moment of the power spectrum,
        assuming a Ricker-like shape."""
        n = _next_pow2(len(self.samples))
        spec = np.abs(np.fft.rfft(self.padded(n))) ** 2
        k = np.arange(1, n // 2)
        mean_index = np.sum(k * spec[1 : n // 2]) / np.sum(spec[1 : n // 2])
        return mean_index / MEAN_TO_PEAK_FREQUENCY_RATIO_RICKER / (n * self.dt)

    def padded(self, n: int) -> np.ndarray:
        """Wavelet wrapped onto a length-``n`` periodic grid, time zero at index 0."""
        out = np.zeros(n)
        idx = (np.arange(len(self.samples)) - self.n_neg) % n
        out[idx] = self.samples
        return out


def load_wavelet(path) -> Wavelet:
    tr = read_su(path)[0]
    if tr.ns < 2 or tr.header["dt"] == 0:
        raise ValueError(f"bad wavelet file {path}")
    dt = tr.dt
    f1 = tr.header["delrt"] * 1e-3 if tr.header["delrt"] != 0 else tr.header["f1"]
    n_neg = -int(round(f1 / dt))
    if not 0 <= n_neg < tr.ns:
        raise ValueError(f"wavelet {path} does not straddle time zero (f1={f1})")
    return Wavelet(np.asarray(tr.data, dtype=float), dt, n_neg)


def _next_pow2(n: int) -> int:
    return 1 << (n - 1).bit_length()


def error_sampling_interval(wavelet: Wavelet) -> tuple[float, float]:
    """Spacing of approximately independent error samples.

    Returns ``(dt_used, correction)``: ``dt_used`` is the largest multiple of
    the wavelet sampling interval not exceeding ``0.253 / f_peak``, and
    ``correction = sqrt(dt_used / dt_ideal)`` rescales the residuals to what
    sampling at the ideal interval would give.
    """
    dt_ideal = SAMPLE_TIME_MAGIC_CONSTANT / wavelet.peak_frequency()
    if dt_ideal < 0.3 * wavelet.dt:
        raise ValueError("wavelet spectrum is odd: independent sample time << dt")
    dt_used = max(wavelet.dt * math.floor(dt_ideal / wavelet.dt), wavelet.dt)
    return dt_used, math.sqrt(dt_used / dt_ideal)


@dataclass
class StackModel:
    """Everything needed to score one stack against one seismic trace."""

    G: np.ndarray  # (n_err, ns) synthetic at error points = G @ R
    observed: np.ndarray  # (n_err,) seismic interpolated to the error points
    scale: float  # correction / noise_rms
    theta_sq_denominator: float  # Dix: theta^2 = vp^2 / this
    wavelet: Wavelet
    stack: Stack


@dataclass
class SeismicModel:
    t_start: float
    t_end: float
    dt: float
    ns: int
    stacks: list[StackModel]
    n_err: int
    correction: float

    @property
    def effective_n(self) -> float:
        """Equivalent number of independent measurements."""
        return self.n_err * self.correction

    def reflectivity(self, t_ref: np.ndarray, rc: np.ndarray) -> np.ndarray:
        """Project spikes ``rc`` at times ``t_ref`` (both (batch, n)) on the grid."""
        batch = t_ref.shape[0]
        R = np.zeros((batch, self.ns))
        j = 1 + np.floor((t_ref - self.t_start) / self.dt).astype(np.int64)
        t = self.t_start + j * self.dt
        t = np.minimum(np.maximum(t, t_ref), t_ref + self.dt)
        ok = (j < self.ns - 1) & (t >= self.t_start) & (t <= self.t_end) & (rc != 0.0)
        wt = np.clip((t - t_ref) / self.dt, 0.0, 1.0)
        # four-point Lagrange weights (Abramowitz & Stegun 25.2.13)
        tmp1 = wt - 2.0
        tmp2 = wt * wt - 1.0
        a_m1 = -wt * (wt - 1.0) * tmp1 / 6.0
        a_0 = tmp2 * tmp1 / 2.0
        a_1 = -wt * (wt + 1.0) * tmp1 / 2.0
        a_2 = wt * tmp2 / 6.0
        spread = ok & (j > 1)
        single = ok & ~spread
        b = np.broadcast_to(np.arange(batch)[:, None], j.shape)
        for off, a in ((-2, a_2), (-1, a_1), (0, a_0), (1, a_m1)):
            np.add.at(R, (b[spread], j[spread] + off), (a * rc)[spread])
        np.add.at(R, (b[single], j[single]), rc[single])
        return R

    def synthetic(self, stack_index: int, R: np.ndarray) -> np.ndarray:
        """Full synthetic trace(s) on the model grid for reflectivity ``R``."""
        w = self.stacks[stack_index].wavelet
        n = _next_pow2(self.ns)
        Rp = np.zeros(R.shape[:-1] + (n,))
        Rp[..., : self.ns] = R
        s = np.fft.irfft(np.fft.rfft(Rp) * np.fft.rfft(w.padded(n)), n)
        return s[..., : self.ns]


def build_seismic_model(
    cfg: ModelConfig, traces: list[Trace], wavelets: list[Wavelet]
) -> SeismicModel:
    """Set up the likelihood for one trace location.

    ``cfg`` must already carry this trace's prior means (the error window is
    placed from the prior mean reflector times); ``traces`` and ``wavelets``
    are given per stack, in the order of ``cfg.stacks``.
    """
    nl = cfg.n_layers
    tops = [lay.time for lay in cfg.layers]
    t_first = min(tops[1:]) if nl > 1 else tops[0]
    t_last = max(tops)
    dt = wavelets[0].dt
    if any(abs(w.dt - dt) > 1e-12 for w in wavelets):
        raise NotImplementedError("wavelets with differing sample intervals")

    t_start = min(t_first - w.reach_after for w in wavelets)
    t_start = math.ceil(t_start * 1e3 - 1e-9) / 1e3
    t_end = max(t_last + w.reach_before for w in wavelets)
    ns = 1 + int(math.ceil((t_end - t_start) / dt - 1e-9))
    n = _next_pow2(ns)

    dt_err, correction = error_sampling_interval(wavelets[0])
    delta_j = int(round(dt_err / dt))

    stack_models = []
    n_err = 0
    for stack, tr, w in zip(cfg.stacks, traces, wavelets):
        t_above = w.reach_after if stack.time_above is None else stack.time_above
        t_below = w.reach_before if stack.time_below is None else stack.time_below
        if t_above > w.reach_after + 1e-12 or t_below > w.reach_before + 1e-12:
            raise ValueError("time_above/time_below exceed the wavelet length")
        t0 = tops[1] - t_above if nl > 1 else tops[0] - t_above
        t1 = tops[-1] + t_below
        n_coarse = int(round((t1 - t0) / dt_err)) + 1
        if n_coarse <= 0:
            raise ValueError("model interval too narrow for a meaningful seismic misfit")
        j_start = max(0, int(round((t0 - t_start) / dt)))
        j = j_start + delta_j * np.arange(n_coarse)
        t = t_start + j_start * dt + dt_err * np.arange(n_coarse)
        # observed seismic, linearly interpolated to the error points
        i_cont = (t - tr.t0) / tr.dt
        if i_cont[0] < 0:
            raise ValueError("seismic trace starts too late for the model interval")
        if math.floor(i_cont[-1]) + 2 >= tr.ns:
            raise ValueError("seismic trace ends too early for the model interval")
        i = i_cont.astype(int)
        wt = i_cont - i
        observed = tr.data[i] * (1.0 - wt) + tr.data[i + 1] * wt
        wp = w.padded(n)
        G = wp[(j[:, None] - np.arange(ns)[None, :]) % n]
        denom = stack.stack_velocity**4 * stack.reflector_time**2 / stack.offset_mean_square
        stack_models.append(
            StackModel(G, observed, correction / stack.noise_rms, denom, w, stack)
        )
        n_err += n_coarse
    return SeismicModel(t_start, t_end, dt, ns, stack_models, n_err, correction)
