#ifndef __POLYAGAMMA_SAMPLER_POMPP_H__
#define __POLYAGAMMA_SAMPLER_POMPP_H__

#include <R.h>
#include <Rmath.h>

double logspace_add(double, double);
double pnorm(double, double, double, int, int);
double runif(double, double);
double rexp(double);
double rnorm(double, double);

#define __PI 3.141592653589793238462643383279502884197
#define LOG4OVERPI 0.24156447527049048367275521675704965061
#define __TRUNC 0.64
#define __TRUNC_RECIP 1.5625

double draw_from_PolyaGamma(double Z);

static inline double mass_texpon(double Z, double fz)
{
  double x0 = log(fz) + fz * __TRUNC;

  return -log1pexp(LOG4OVERPI + logspace_add(
      x0 - Z + pnorm(sqrt(1.0 / __TRUNC) * (__TRUNC * Z - 1), 0, 1, 1, 1),
      x0 + Z + pnorm(sqrt(1.0 / __TRUNC) * (__TRUNC * Z + 1) * -1.0, 0, 1, 1, 1)
  ));
}

static inline double rtigauss(double Z)
{
  Z = fabs(Z);
  double X = __TRUNC + 1.0;
  if (__TRUNC_RECIP > Z) { // mu > t
    double alpha = 0.0;
    while (runif(0,1) > alpha) {
      // X = t + 1.0;
      // while (X > t)
      // 	X = 1.0 / gamma_rate(0.5, 0.5);
      // Slightly faster to use truncated normal.
      double E1 = rexp(1.0); double E2 = rexp(1.0);
      while ( E1*E1 > 2 * E2 / __TRUNC) {
        E1 = rexp(1.0); E2 = rexp(1.0);
      }
      X = 1 + E1 * __TRUNC;
      X = __TRUNC / (X * X);
      alpha = exp(-0.5 * X * (Z*Z));
    }
  }
  else {
    double mu = 1.0 / Z;
    while (X > __TRUNC) {
      double Y = rnorm(0., 1.); Y *= Y;
      double half_mu = 0.5 * mu;
      double mu_Y    = mu  * Y;
      X = mu + half_mu * mu_Y - half_mu * sqrt(4 * mu_Y + mu_Y * mu_Y);
      if (runif(0,1) > mu / (mu + X))
        X = mu*mu / X;
    }
  }
  return X;
}

static inline double a(int n, double x)
{
  double K = (n + 0.5) * __PI;
  if (x > __TRUNC) {
    return K * exp( -0.5 * K*K * x );
  }
  else if (x > 0) {
    double expnt = -1.5 * (log(0.5 * __PI)  + log(x)) + log(K) - 2.0 * (n+0.5)*(n+0.5) / x;
    return exp(expnt);
  }
  return 0.;
}

#endif
