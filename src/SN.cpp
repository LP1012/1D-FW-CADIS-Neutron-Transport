
#include "SN.h"
#include "Cell.h"
#include "SNCell.h"
#include "CellQuadrature.h"

#include <stdexcept>
#include <vector>
#include <cmath>

template <std::size_t N>
SN<N>::SN(const std::vector<Cell> & cells) : _num_cells(cells.size())
{
  if (GQ_order % 2 == 1.0)
    throw std::runtime_error("Gauss quadrature order cannot be odd!");
  populateSNCells(cells);
  _mus = discreteQuadrature::getAbscissaAsVector<N>();
  _k = 1.0; // assume critical just for initial guess
}

template <std::size_t N>
void
SN<N>::populateSNCells(const std::vector<Cell> & cells)
{
  for (auto & cell : cells)
    _sn_cells.push_back(SNCell(cell));
  normalizeSources(); // perform initial normalization
}

template <std::size_t N>
void
SN<N>::computeScalarFluxAll()
{
  for (auto & cell : _sn_cells)
    cell.computeScalarFlux();
}

template <std::size_t N>
void
SN<N>::sweepLeft(const unsigned int mu_index)
{
  const double mu = _mus[mu_index];
  static_assert(mu > 0);
  double left_flux = 0; // vacuum BC on left side
  for (auto i = 0; i < _num_cells; i++)
  {
    double new_angular_flux = computeAngularFlux(_sn_cells[i], left_flux, mu);
    _sn_cells[i].setAngularFlux(new_angular_flux, mu_index);
    double right_flux = computeCellFlux(new_angular_flux, left_flux);
    left_flux = right_flux;
  }
}

template <std::size_t N>
void
SN<N>::sweepRight(const unsigned int mu_index)
{
  const double mu = _mus[mu_index];
  static_assert(mu < 0);
  double right_flux = 0; // vacuum BC on right side
  for (auto i = _num_cells - 1; i > -1; i--)
  {
    double new_angular_flux = computeAngularFlux(_sn_cells[i], right_flux, mu);
    _sn_cells[i].setAngularFlux(new_angular_flux, right_flux);
    double left_flux = computeCellFlux(new_angular_flux, right_flux);
    right_flux = left_flux;
  }
}

template <std::size_t N>
double
SN<N>::computeAngularFlux(const SNCell<N> & cell, const double cell_flux, const double mu)
{
  double sigma_t = cell.sigmaT();
  double delta = cell.cellWidth();
  double source = cell.source();
  double numerator = cell_flux + delta * source / (2.0 * std::abs(mu));
  double denom = 1.0 + sigma_t * delta / (2.0 * std::abs(mu));
  return numerator / denom;
}

template <std::size_t N>
double
SN<N>::computeCellFlux(const double cell_centered_flux, const double known_flux)
{
  return 2.0 * cell_centered_flux - known_flux;
}

template <std::size_t N>
void
SN<N>::normalizeSources()
{
  // integrate all sources
  double total_source = 0;
  for (auto & cell : _sn_cells)
    total_source += cell.cellWidth() * cell.source(); // approximate with a midpoint rule

  for (auto & cell : _sn_cells)
  {
    double old_source = cell.source();
    cell.setScalarFlux(old_source / total_source);
  }
}

template <std::size_t N>
std::vector<double>
SN<N>::getScalarFlux()
{
  std::vector<double> scalar_flux;
  for (auto & cell : _sn_cells)
    scalar_flux.push_back(cell.scalarFlux());
  return scalar_flux;
}

template <std::size_t N>
void
SN<N>::updateK()
{
  double old_fission_contrib = integrateFissionSource<N>(_sn_cells);

  computeScalarFluxAll<N>(); // updates to new scalar flux values in _sn_cells
  double new_fission_contrib = integrateFissionSource<N>(_sn_cells) double old_k = _k;
  _k = old_k * new_fission_contrib / old_fission_contrib;
}

template <std::size_t N>
double
SN<N>::integrateFissionSource(const std::vector<SNCell<N>> & cells)
{
  double running_sum = 0.0;
  for (auto i = 0; i < cells.size(); i++)
  {
    double dx = cells[i].cellWidth();
    double nu_sigma_f = cells[i].nuSigmaF();
    double scalar_flux = cell[i].scalarFlux();
    running_sum += dx * scalar_flux * nu_sigma_f; // approximate midpoint rule
  }
  return running_sum;
}

template <std::size_t N>
void
SN<N>::updateSource()
{
  for (auto & cell : _sn_cells)
  {
    double flux = cell.scalarFlux();
    double scattering_xs = cell.sigmaS();
    double nu_sigma_f = cell.nuSigmaF();
    double new_source = 0.5 * flux * (scattering_xs + 1.0 / _k * nu_sigma_f);
    cell.setSource(new_source);
  }
  normalizeSources(); // normalize to single source particle
}