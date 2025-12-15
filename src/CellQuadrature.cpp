#include "CellQuadrature.h"

namespace discreteQuadrature
{
double
integrate(const std::vector<double> & fs,
          const GaussLegendreRule & quad,
          const double a,
          const double b)
{
  if (fs.size() != quad.abscissa.size())
    throw std::runtime_error("Number of gridpoints given dows not match quadrature values!");

  const double front_coeff = (b - a) / 2.0;
  double running_sum = 0;
  for (auto i = 0; i < fs.size(); i++)
    running_sum += fs[i] * quad.weights[i];

  return front_coeff * running_sum;
}
}