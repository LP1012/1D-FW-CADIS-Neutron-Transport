#include "SNCell.h"

SNCell::SNCell(const Cell & cell)
  : CellBase(cell.xMin(), cell.xMax(), cell.sigmaA(), cell.sigmaS(), cell.nuSigmaF())
{
}

SNCell::SNCell(const double xmin,
               const double xmax,
               const double Sigma_a,
               const double Sigma_s,
               const double nu_Sigma_f)
  : CellBase(xmin, xmax, Sigma_a, Sigma_s, nu_Sigma_f)
{
  if (Sigma_s + nu_Sigma_f > 1e-15)
    _source = 1.0; // set dummy source value if scattering or fission is present
}