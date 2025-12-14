#pragma once

#include "Cell.h"

class SNCell : public CellBase
{
public:
  SNCell(const Cell & cell);
  SNCell(const double xmin,
         const double xmax,
         const double Sigma_a,
         const double Sigma_s,
         const double nu_Sigma_f);

protected:
  double _scalar_flux;
  std::vector<double> _angular_fluxes;
  double _source;
};