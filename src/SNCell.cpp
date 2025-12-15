#include "SNCell.h"

SNCell::SNCell(const Cell & cell, const unsigned int gl_quad)
  : CellBase(cell.xMin(), cell.xMax(), cell.sigmaA(), cell.sigmaS(), cell.nuSigmaF()),
    _gl_quad(gl_quad)
{
  initializeAngularFluxVector();
}

SNCell::SNCell(const double xmin,
               const double xmax,
               const double Sigma_a,
               const double Sigma_s,
               const double nu_Sigma_f,
               const unsigned int gl_quad)
  : CellBase(xmin, xmax, Sigma_a, Sigma_s, nu_Sigma_f), _gl_quad(gl_quad)
{
  initializeAngularFluxVector();
  if (Sigma_s + nu_Sigma_f > 1e-15)
    _source = 1.0; // set dummy source value if scattering or fission is present
}

void
SNCell::initializeAngularFluxVector()
{
  for (auto i = 0; i < _gl_quad; i++)
    _angular_fluxes.push_back(0.0); // just initialize angular fluxes to 0s
}