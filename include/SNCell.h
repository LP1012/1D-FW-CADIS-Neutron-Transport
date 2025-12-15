#pragma once

#include "Cell.h"

class SNCell : public CellBase
{
public:
  SNCell(const Cell & cell, const unsigned int gl_quad);
  SNCell(const double xmin,
         const double xmax,
         const double Sigma_a,
         const double Sigma_s,
         const double nu_Sigma_f,
         const unsigned int gl_quad);

  void setScalarFlux(const double new_flux) { _scalar_flux = new_flux; }
  void setSource(const double new_source) { _source = new_source; }
  void setAngularFlux(const double new_angular_flux, const unsigned int index)
  {
    _angular_fluxes[index] = new_angular_flux;
  }

  const double source() const { return _source; }
  const double scalarFlux() const { return _scalar_flux; }

protected:
  const unsigned int _gl_quad;
  double _scalar_flux;
  std::vector<double> _angular_fluxes;
  double _source;

  void initializeAngularFluxVector();
};