#include "Cell.h"

Cell::Cell(const double xmin, const double xmax) : _xmin(xmin), _xmax(xmax) {}

void Cell::setWeight(const double weight)
{
  _weight = weight;
}

void Cell::setForwardFlux(const double forward_flux)
{
  _forward_flux = forward_flux;
}

void Cell::setAdjointFlux(const double adjoint_flux)
{
  _adjoint_flux = adjoint_flux;
}