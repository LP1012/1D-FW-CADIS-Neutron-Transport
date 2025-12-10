#include "Cell.h"
#include <cmath>
#include <vector>
#include <stdexcept>

Cell::Cell(const double xmin,
           const double xmax,
           const double Sigma_a,
           const double Sigma_s,
           const double nu_Sigma_f)
  : _xmin(xmin), _xmax(xmax), _Sigma_a(Sigma_a), _Sigma_s(Sigma_s), _nu_Sigma_f(nu_Sigma_f)
{
  if (_xmax < _xmin)
    throw std::runtime_error("Max cell size is less then minimum cell size!");
  _cell_width = std::abs(xmax - xmin);
  _cell_center = (_xmax + _xmin) / 2.0;
  _Sigma_t = _Sigma_a + _Sigma_s;
  _n_per_abs = _nu_Sigma_f / _Sigma_a;
  _abs_ratio = _Sigma_a / _Sigma_t;
}

unsigned int
Cell::cellIndex(const double position, const double mu, const std::vector<Cell> cells)
{
  for (auto i = 0; i < cells.size(); i++)
  {
    const double cell_max = cells[i].xMax();
    const double cell_min = cells[i].xMin();
    if (position < cell_max)
      return i;
    else if (isOnBoundary(position, cell_max) && mu < 0)
      return i;
    else if (isOnBoundary(position, cell_min) && mu > 0)
      return i;
  }
  throw std::runtime_error("Position not located within given cells!");
}

void
Cell::setWeight(const double weight, const double upper_weight, const double lower_weight)
{
  _center_weight = weight;
  _upper_weight = upper_weight;
  _lower_weight = lower_weight;
}

bool
Cell::isOnBoundary(const double pos, const double boundary)
{
  const double error = std::abs(pos - boundary);
  const double tol = 1e-15;
  if (error < tol)
    return true;
  else
    return false;
}