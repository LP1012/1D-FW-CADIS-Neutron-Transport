
#include "SN.h"
#include "Cell.h"
#include "SNCell.h"
#include "CellQuadrature.h"

#include <stdexcept>
#include <vector>
#include <cmath>
#include <fstream>
#include <iostream>
#include <cassert>

SN::SN(const std::vector<Cell> & cells, const unsigned int GQ_order)
  : _num_cells(cells.size()), _gq_order(GQ_order)
{
  if (_gq_order % 2 == 1.0)
    throw std::runtime_error("Gauss quadrature order cannot be odd!");
  populateSNCells(cells);
  _gauss_legendre_rule = discreteQuadrature::getGaussLegendreRule(_gq_order);
  _mus = _gauss_legendre_rule.abscissa;
  _k = 1.0; // assume critical just for initial guess
  _is_converged = false;
}

void
SN::run()
{
  unsigned int safety = 0;
  while (!_is_converged && safety <= 100)
  {
    for (auto i = 0; i < _mus.size(); i++)
    {
      assert(std::abs(_mus[i]) > 1e-15);
      if (_mus[i] < 0)
        sweepLeft(i); // sweep left
      else
        sweepRight(i); // sweep right
    }
    // save old values for convergence tests
    std::vector<double> old_scalar_flux = getScalarFlux();
    const double old_k = _k;

    updateK();          // update k and scalar fluxes
    updateSource();     // update source
    normalizeSources(); // normalize to single source particle

    // set new values to variables for convenience
    const double new_k = _k;
    std::vector<double> new_scalar_flux = getScalarFlux();

    safety++;

    printf("|   %d    |", safety);
    // check for convergence
    _is_converged = isConverged(old_scalar_flux, new_scalar_flux, old_k, new_k);
  }

  printf("\nFinal k-eff: %.6f\n", _k);

  printf("\nExporting results... ");
  exportToCsv(); // export results to csv
  printf("Done.\n");
}

void
SN::exportToCsv()
{
  std::string output = "SN_output";
  std::ofstream outfile;
  outfile.open(output);
  if (!outfile.is_open())
    throw std::runtime_error("SN output file not opened successfully!");

  outfile << "position,scalar_flux" << std::endl;
  for (auto cell : _sn_cells)
    outfile << cell.cellCenter() << "," << cell.scalarFlux() << std::endl;

  outfile.close();
}

bool
SN::isConverged(const std::vector<double> & old_flux,
                const std::vector<double> & new_flux,
                const double old_k,
                const double new_k)
{
  double relative_k = std::abs(new_k - old_k) / new_k;

  std::vector<double> error_vector;
  for (auto i = 0; i < old_flux.size(); i++)
    error_vector.push_back(new_flux[i] - old_flux[i]);
  double error_vector_norm = L2Norm(error_vector);
  double new_flux_norm = L2Norm(new_flux);
  double relative_flux = error_vector_norm / new_flux_norm;

  printf("    |   %.4e    |   %.4e    |\n", relative_k, relative_flux);

  double tol = 1e-6;
  if (relative_k < tol && relative_flux < tol)
    return true;
  else
    return false;
}

double
SN::L2Norm(const std::vector<double> & vector)
{
  double running_sum = 0.0;
  for (auto val : vector)
    running_sum += val * val;
  return std::sqrt(running_sum);
}

void
SN::populateSNCells(const std::vector<Cell> & cells)
{
  for (auto & cell : cells)
    _sn_cells.push_back(SNCell(cell, _gq_order));
  normalizeSources(); // perform initial normalization

  // set initial flux guess
  for (auto & cell : _sn_cells)
    cell.setScalarFlux(1.0);
}

void
SN::computeScalarFluxAll()
{
  for (auto & cell : _sn_cells)
    cell.computeScalarFlux();
}

void
SN::sweepRight(const unsigned int mu_index)
{
  const double mu = _mus[mu_index];
  assert(mu > 0);
  double left_flux = 0; // vacuum BC on left side
  for (auto i = 0; i < _num_cells; i++)
  {
    double new_angular_flux = computeAngularFlux(_sn_cells[i], left_flux, mu);
    _sn_cells[i].setAngularFlux(new_angular_flux, mu_index);
    double right_flux = computeCellFlux(new_angular_flux, left_flux);
    left_flux = right_flux;
  }
}

void
SN::sweepLeft(const unsigned int mu_index)
{
  const double mu = _mus[mu_index];
  assert(mu < 0);
  double right_flux = 0; // vacuum BC on right side
  for (int i = _num_cells - 1; i > -1; i--)
  {
    double new_angular_flux = computeAngularFlux(_sn_cells[i], right_flux, mu);
    _sn_cells[i].setAngularFlux(new_angular_flux, mu_index);
    double left_flux = computeCellFlux(new_angular_flux, right_flux);
    right_flux = left_flux;
  }
}

double
SN::computeAngularFlux(const SNCell & cell, const double cell_flux, const double mu)
{
  double sigma_t = cell.sigmaT();
  double delta = cell.cellWidth();
  double source = cell.source();
  double numerator = cell_flux + delta * source / (2.0 * std::abs(mu));
  double denom = 1.0 + sigma_t * delta / (2.0 * std::abs(mu));
  return numerator / denom;
}

double
SN::computeCellFlux(const double cell_centered_flux, const double known_flux)
{
  return 2.0 * cell_centered_flux - known_flux;
}

void
SN::normalizeSources()
{
  // integrate all sources
  double total_source = 0;
  for (auto & cell : _sn_cells)
    total_source += cell.cellWidth() * cell.source(); // approximate with a midpoint rule

  for (auto & cell : _sn_cells)
  {
    double old_source = cell.source();
    cell.setSource(old_source / total_source);
  }
}

std::vector<double>
SN::getScalarFlux()
{
  std::vector<double> scalar_flux;
  for (auto & cell : _sn_cells)
    scalar_flux.push_back(cell.scalarFlux());
  return scalar_flux;
}

void
SN::updateK()
{
  double old_fission_contrib = integrateFissionSource(_sn_cells);

  computeScalarFluxAll(); // updates to new scalar flux values in _sn_cells
  double new_fission_contrib = integrateFissionSource(_sn_cells);
  double old_k = _k;
  _k = old_k * new_fission_contrib / old_fission_contrib;
}

double
SN::integrateFissionSource(const std::vector<SNCell> & cells)
{
  double running_sum = 0.0;
  for (auto i = 0; i < cells.size(); i++)
  {
    double dx = cells[i].cellWidth();
    double nu_sigma_f = cells[i].nuSigmaF();
    double scalar_flux = cells[i].scalarFlux();
    running_sum += dx * scalar_flux * nu_sigma_f; // approximate midpoint rule
  }
  return running_sum;
}

void
SN::updateSource()
{
  for (auto & cell : _sn_cells)
  {
    double flux = cell.scalarFlux();
    double scattering_xs = cell.sigmaS();
    double nu_sigma_f = cell.nuSigmaF();
    double new_source = 0.5 * flux * (scattering_xs + 1.0 / _k * nu_sigma_f);
    cell.setSource(new_source);
  }
}

std::vector<double>
SN::getSources()
{
  std::vector<double> sources;
  for (auto & cell : _sn_cells)
    sources.push_back(cell.source());
  return sources;
}