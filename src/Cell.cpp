#include "Cell.h"
#include <cmath>
#include <vector>
#include <stdexcept>
#include "RNG.h"

Cell::Cell(const double xmin,
           const double xmax,
           const double Sigma_a,
           const double Sigma_s,
           const double nu_Sigma_f,
           const double source)
  : CellBase(xmin, xmax, Sigma_a, Sigma_s, nu_Sigma_f, source)
{
  _n_per_abs = _nu_Sigma_f / _Sigma_a;
  _abs_ratio = _Sigma_a / _Sigma_t;
  _abs_prob = _Sigma_a / _Sigma_t;
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
Cell::createWeightWindow(const double fwcadis_adjoint_flux, const double window_width)
{
  double denom = fwcadis_adjoint_flux * (window_width + 1.0) / 2.0;
  _lower_weight = 1.0 / denom;
  _upper_weight = window_width * _lower_weight;
  _center_weight = (_lower_weight + _upper_weight) / 2.0;
}
