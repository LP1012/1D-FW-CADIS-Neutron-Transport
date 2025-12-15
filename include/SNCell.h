#pragma once

#include "Cell.h"

class SNCell : public CellBase
{
public:
  SNCell(const Cell & cell, const unsigned int quadrature_order);
  SNCell(const double xmin,
         const double xmax,
         const double Sigma_a,
         const double Sigma_s,
         const double nu_Sigma_f,
         const unsigned int quadrature_order);

  void computeScalarFlux();

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
  const unsigned int _quad;
  double _scalar_flux;
  std::vector<double> _angular_fluxes;
  double _source;
};