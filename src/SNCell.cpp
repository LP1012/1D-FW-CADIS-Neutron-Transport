#include "SNCell.h"
#include "CellQuadrature.h"

template <std::size_t N>
SNCell<N>::SNCell(const Cell & cell)
  : CellBase(cell.xMin(), cell.xMax(), cell.sigmaA(), cell.sigmaS(), cell.nuSigmaF()),
    _angular_fluxes(N, 0.0)
{
  if (_Sigma_s + _nu_Sigma_f > 1e-15)
    _source = 1.0; // set dummy source value if scattering or fission is present
}

template <std::size_t N>
SNCell<N>::SNCell(const double xmin,
                  const double xmax,
                  const double Sigma_a,
                  const double Sigma_s,
                  const double nu_Sigma_f)
  : CellBase(xmin, xmax, Sigma_a, Sigma_s, nu_Sigma_f), _angular_fluxes(N, 0.0)
{
  if (_Sigma_s + _nu_Sigma_f > 1e-15)
    _source = 1.0; // set dummy source value if scattering or fission is present
}

template <std::size_t N>
void
SNCell<N>::computeScalarFlux()
{
  _scalar_flux = 2.0 * M_PI * discreteQuadrature::integrate<N>(_angular_fluxes, -1.0, 1.0);
}