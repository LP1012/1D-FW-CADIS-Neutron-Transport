#include "CellBase.h"
#include <cmath>
#include <vector>

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
