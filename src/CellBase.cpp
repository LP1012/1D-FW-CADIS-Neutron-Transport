#include "CellBase.h"
#include <cmath>
#include <vector>
#include <stdexcept>

CellBase::CellBase(const double xmin,
                   const double xmax,
                   const double Sigma_a,
                   const double Sigma_s,
                   const double nu_Sigma_f,
                   const double volumetric_source)
  : _xmin(xmin),
    _xmax(xmax),
    _Sigma_a(Sigma_a),
    _Sigma_s(Sigma_s),
    _nu_Sigma_f(nu_Sigma_f),
    _vol_source(volumetric_source)
{
  if (_xmax < _xmin)
    throw std::runtime_error("Max cell size is less then minimum cell size!");
  _cell_width = std::abs(xmax - xmin);
  _cell_center = (_xmax + _xmin) / 2.0;
  _Sigma_t = _Sigma_a + _Sigma_s;
}

template <typename T>
unsigned int
CellBase::cellIndex(const double position, const double mu, const std::vector<T> & cells)
{
  for (auto i = 0; i < cells.size(); i++)
  {
    const double cell_max = cells[i].xMax();
    const double cell_min = cells[i].xMin();
    if (position < cell_max && position > cell_min)
      return i;
    else if (isOnBoundary(position, cell_max) && mu < 0)
      return i;
    else if (isOnBoundary(position, cell_min) && mu > 0)
      return i;
  }
  throw std::runtime_error("Position not located within given cells!");
}

bool
CellBase::isOnBoundary(const double pos, const double boundary)
{
  const double error = std::abs(pos - boundary);
  const double tol = 1e-15;
  if (error < tol)
    return true;
  else
    return false;
}
