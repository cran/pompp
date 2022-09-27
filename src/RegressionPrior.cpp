#include <RcppEigen.h>
#include "include/RegressionPrior.hpp"
#ifdef _OPENMP
#include <omp.h>
#endif
#include "include/safeR.hpp"

Eigen::VectorXd NormalPrior::sample(const Eigen::VectorXd& mean,
                                    const Eigen::MatrixXd& precision) {
  sigmaSolver.compute(precision + priorPrecision);

  return sigmaSolver.matrixU().solve(rnorm(mean.size())) +
    sigmaSolver.solve(mean + precisionTimesMean);
}

inline double NormalPrior::logPrior(const Eigen::VectorXd& betas) {
  Eigen::VectorXd m = betas - priorMean;
  return m.transpose() * sigmaSolver.solve(m);
}
