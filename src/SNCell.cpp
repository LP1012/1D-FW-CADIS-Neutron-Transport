#include "SNCell.h"
#include "CellQuadrature.h"

SNCell::SNCell(const Cell & cell, const unsigned int quadrature_order)
  : CellBase(cell.xMin(), cell.xMax(), cell.sigmaA(), cell.sigmaS(), cell.nuSigmaF()),
    _quad(quadrature_order),
    _angular_fluxes(_quad, 0.0),
    _scalar_flux(0)
{
  if (_Sigma_s + _nu_Sigma_f > 1e-15)
    _source = 1.0; // set dummy source value if scattering or fission is present
}

SNCell::SNCell(const double xmin,
               const double xmax,
               const double Sigma_a,
               const double Sigma_s,
               const double nu_Sigma_f,
               const unsigned int quadrature_order)
  : CellBase(xmin, xmax, Sigma_a, Sigma_s, nu_Sigma_f),
    _quad(quadrature_order),
    _angular_fluxes(_quad, 0.0),
    _scalar_flux(0)
{
  if (_Sigma_s + _nu_Sigma_f > 1e-15)
    _source = 1.0; // set dummy source value if scattering or fission is present
}

void
SNCell::computeScalarFlux()
{
  const auto & quad = discreteQuadrature::getGaussLegendreRule(_quad);
  _scalar_flux = 2.0 * M_PI * discreteQuadrature::integrate(_angular_fluxes, quad, -1.0, 1.0);
}