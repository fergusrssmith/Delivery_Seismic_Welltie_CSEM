import numpy as np
import pytest

from helpers import CONTACTS, NG_WEDGE
from pydelivery import Inversion
from pydelivery.io import (TraceSet, ibm_to_float, import_java_project, read_prior_table,
                           read_results, read_su, read_traces, write_traces)


@pytest.mark.parametrize("ext", ["h5", "npz", "sgy"])
def test_trace_roundtrip(tmp_path, ext):
    tr = read_su(NG_WEDGE / "near-stack.su")
    write_traces(tmp_path / f"x.{ext}", tr)
    back = read_traces(tmp_path / f"x.{ext}")
    np.testing.assert_allclose(back.data, tr.data, rtol=1e-6)
    np.testing.assert_allclose(back.t0, tr.t0)
    assert back.dt == pytest.approx(tr.dt)
    assert back.headers["tracl"][3] == tr.headers["tracl"][3]


def test_ibm_float():
    words = np.array([0xC276A000, 0x42640000, 0x00000000, 0x41100000], dtype=np.uint32)
    np.testing.assert_allclose(ibm_to_float(words), [-118.625, 100.0, 0.0, 1.0])


def test_segy_rejects_sub_ms_start(tmp_path):
    with pytest.raises(ValueError):
        write_traces(tmp_path / "x.sgy", TraceSet(np.zeros((1, 4)), 0.004, 0.0005))


@pytest.mark.parametrize("fmt", ["h5", "npz", "sgy"])
def test_import_java_project(tmp_path, fmt):
    new_xml = import_java_project(CONTACTS / "ModelDescription.xml", tmp_path,
                                  CONTACTS / "prior_traces.su", fmt)
    a = Inversion.from_files(CONTACTS / "ModelDescription.xml", CONTACTS / "prior_traces.su")
    b = Inversion.from_files(new_xml, tmp_path / ("priors.npz" if fmt == "npz" else "priors.h5"))
    np.testing.assert_allclose(a.stacks[0].data, b.stacks[0].data, rtol=1e-6)
    np.testing.assert_allclose(a.wavelets[0].samples, b.wavelets[0].samples, rtol=1e-6)
    assert a.wavelets[0].n_neg == b.wavelets[0].n_neg
    for k in a.priors:
        np.testing.assert_allclose(a.priors[k], b.priors[k], equal_nan=True)


def test_java_priors_respect_varies_areally():
    # GasOilFlank links probability_of_oil but marks it varies_areally="false":
    # Java ignores the prior-trace values, so must we
    from pydelivery import load_config
    cfg = load_config(NG_WEDGE / "ModelDescription.xml")
    table = read_prior_table(NG_WEDGE / "prior_traces.su", cfg)
    assert set(table) == {"time", "net_to_gross"}
    assert np.isnan(table["net_to_gross"][0, 0]) and table["net_to_gross"][7, 1] > 0


def test_results_file(tmp_path):
    inv = Inversion.from_files(NG_WEDGE / "ModelDescription.xml", NG_WEDGE / "prior_traces.su")
    res = inv.run(20, traces=[0, 5], seed=0, burn_in=50, thin=2, n_chains=10)
    for ext in ("h5", "npz"):
        inv.save(tmp_path / f"r.{ext}", res)
        back = read_results(tmp_path / f"r.{ext}")
        assert back["samples"]["NG"].shape == (2, 20, 3)
        assert back["samples"]["t_base"].shape == (2, 20)
        assert back["fluid_states"] == [["bbb"], ["bbb"]]
        assert back["layer_names"][1] == "reservoir"
        np.testing.assert_allclose(back["samples"]["t"][1], res[1].samples["t"])
