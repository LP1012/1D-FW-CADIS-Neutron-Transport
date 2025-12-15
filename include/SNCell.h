#pragma once

#include "Cell.h"

template <std::size_t N>
class SNCell : public CellBase
{
public:
  SNCell(const Cell & cell);
  SNCell(const double xmin,
         const double xmax,
         const double Sigma_a,
         const double Sigma_s,
         const double nu_Sigma_f);

  void setScalarFlux(const double new_flux) { _scalar_flux = new_flux; }
  void setSource(const double new_source) { _source = new_source; }
  void setAngularFlux(const double new_angular_flux, const unsigned int index)
  {
    _angular_fluxes[index] = new_angular_flux;
  }

  const double source() const { return _source; }
  const double scalarFlux() const { return _scalar_flux; }
  const std::vector<double> angularFluxes() const { return _angular_fluxes; }

protected:
  double _scalar_flux;
  std::vector<double> _angular_fluxes;
  double _source;
};