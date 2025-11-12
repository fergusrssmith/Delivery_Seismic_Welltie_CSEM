package au.csiro.JamesGunning.lib.Minpack.optimization;

interface Lmdif_fcn {

   void fcn(int m, int n, double x[], double fvec[],
            int iflag[]);

}