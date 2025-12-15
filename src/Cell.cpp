#include "Cell.h"
#include <cmath>
#include <vector>
#include <stdexcept>
#include "RNG.h"

Cell::Cell(const double xmin,
           const double xmax,
           const double Sigma_a,
           const double Sigma_s,
           const double nu_Sigma_f)
  : CellBase(xmin, xmax, Sigma_a, Sigma_s, nu_Sigma_f)
{
  _n_per_abs = _nu_Sigma_f / _Sigma_a;
  _abs_ratio = _Sigma_a / _Sigma_t;
}

double
Cell::randomPositionInCell(const Cell cell, UniformRNG rng)
{
  double xmax = cell.xMax();
  double xmin = cell.xMin();
  double rn = rng.generateRN();
  return rn * (xmax * xmin) + xmin;
}

void
Cell::setWeight(const double weight, const double upper_weight, const double lower_weight)
{
  _center_weight = weight;
  _upper_weight = upper_weight;
  _lower_weight = lower_weight;
}
