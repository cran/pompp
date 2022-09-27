#include "include/PolyaGammaSampler.h"
#include <R.h>
#include <Rmath.h>

// Algorithms copied from https://github.com/jwindle/BayesLogit/tree/master/src

double draw_from_PolyaGamma(double Z)
{
  // Change the parameter.
  Z = fabs(Z) * 0.5;

  // Now sample 0.25 * J^*(1, Z := Z/2).
  double fz = 0.125 * __PI*__PI + 0.5 * Z*Z;
  // ... Problems with large Z?  Try using q_over_p.
  // double p  = 0.5 * __PI * exp(-1.0 * fz * __TRUNC) / fz;
  // double q  = 2 * exp(-1.0 * Z) * pigauss(__TRUNC, Z);

  double X;
  double S;
  double Y;
  // int iter = 0; If you want to keep track of iterations.

  while (1) {

    // if (unif() < p/(p+q))
    if ( log(runif(0,1)) < mass_texpon(Z, fz) )
      X = __TRUNC + rexp(1.) / fz;
    else
      X = rtigauss(Z);

    S = a(0, X);
    Y = runif(0,1) * S;
    int n = 0;
    int go = 1;

    // Cap the number of iterations?
    while (go) {

      // Break infinite loop.  Put first so it always checks n==0.
      // Has problem with OpenMP
      // if (n % 1000 == 0) R_CheckUserInterrupt();

      ++n;
      if (n%2==1) {
        S -= a(n, X);
        if ( Y<=S ) return 0.25 * X;
      }
      else {
        S += a(n, X);
        if ( Y>S ) go = 0;
      }

    }
    // Need Y <= S in event that Y = S, e.g. when X = 0.

  }
} // draw_like_devroye

