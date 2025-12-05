#include "Cell.h"
#include <cmath>
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
}

void
Cell::setWeight(const double weight, const double upper_weight, const double lower_weight)
{
  _center_weight = weight;
  _upper_weight = upper_weight;
  _lower_weight = lower_weight;
}

void
Cell::setForwardFlux(const double forward_flux)
{
  _forward_flux = forward_flux;
}

void
Cell::setAdjointFlux(const double adjoint_flux)
{
  _adjoint_flux = adjoint_flux;
}