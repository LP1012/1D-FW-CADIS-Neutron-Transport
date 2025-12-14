#pragma once

#include "Cell.h"

class SNCell : public CellBase
{
public:
  using CellBase::CellBase;

protected:
  double _scalar_flux;
  std::vector<double> _angular_fluxes;
  double _source;
};