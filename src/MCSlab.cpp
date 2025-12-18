#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "XMLUtils.h"
#include "tinyxml2.h"
#include "utils.h"

#include <algorithm> // for std::sort
#include <array>
#include <cassert>
#include <cmath>
#include <fstream>
#include <iostream>
#include <random>
#include <vector>

MCSlab::MCSlab(const std::string input_file_name,
               const unsigned int n_particles,
               const unsigned int n_generations,
               const unsigned int n_inactive,
               const std::vector<Region> regions,
               const std::vector<Cell> cells,
               const bool implicit_capture,
               const bool use_weight_windows)
  : _input_file_name(input_file_name),
    _n_particles(n_particles),
    _n_generations(n_generations),
    _n_inactive(n_inactive),
    _regions(regions),
    _cells(cells),
    _implicit_capture(implicit_capture),
    _use_wws(use_weight_windows),
    _rng()
{
  _n_total_cells = _cells.size();
  _n_fissionable_regions = 0;
  createFissionCells();
}

void
MCSlab::createFissionCells()
{
  for (auto & cell : _cells)
  {
    if (cell.nuSigmaF() > 1e-15)
      _fissionable_cells.push_back(cell);
  }
}

void
MCSlab::initializeOutput()
{
  std::string outfile_name = _input_file_name;
  removeSuffix(outfile_name, ".xml");

  std::string collision_outfile_name = outfile_name + "_col.csv";
  std::string pl_outfile_name = outfile_name + "_pl.csv";

  printf("Creating output files...\n");
  printf("  Collision:  %s\n", collision_outfile_name.c_str());
  printf("  Pathlength: %s\n\n", pl_outfile_name.c_str());
  _collision_outfile.open(collision_outfile_name);
  _pathlength_outfile.open(pl_outfile_name);

  if (!_collision_outfile.is_open())
    throw std::runtime_error("Collision output file not opened successfully!");
  if (!_pathlength_outfile.is_open())
    throw std::runtime_error("Pathlength output file not opened successfully!");

  _collision_outfile << "position,weight" << std::endl;
  _pathlength_outfile << "start,end,mu,pathlength,weight" << std::endl;

  exportRegionsToCsv(outfile_name);
}

void
MCSlab::k_eigenvalue()
{
  // Create code head
  printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
         "\n");
  printf("* _____ ______   ________  ________  ___       ________  ________     "
         "*\n"
         "*|\\   _ \\  _   \\|\\   ____\\|\\   ____\\|\\  \\     |\\   __  \\|\\  "
         " "
         "__  \\    *\n"
         "*\\ \\  \\\\\\__\\ \\  \\ \\  \\___|\\ \\  \\___|\\ \\  \\    \\ \\  "
         "\\|\\  \\ \\  \\|\\ /_   *\n*"
         " \\ \\  \\\\|__| \\  \\ \\  \\    \\ \\_____  \\ \\  \\    \\ \\   __  "
         "\\ \\   __  \\  *\n*"
         "  \\ \\  \\    \\ \\  \\ \\  \\____\\|____|\\  \\ \\  \\____\\ \\  \\ "
         "\\  \\ \\  \\|\\  \\ *\n*"
         "   \\ \\__\\    \\ \\__\\ \\_______\\____\\_\\  \\ \\_______\\ \\__\\ "
         "\\__\\ \\_______\\*\n*"
         "    \\|__|     "
         "\\|__|\\|_______|\\_________\\|_______|\\|__|\\|__|\\|_______|*\n*"
         "                            \\|_________|                             "
         "*\n*"
         "                                                                     "
         "*\n");
  printf("* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *"
         "\n");

  printf("\nSimulation Specification:\n");
  printf("    Number of Regions:        %lu\n", _regions.size());
  printf("    Particles per Generation: %d\n", _n_particles);
  printf("    Number of Generations:    %d\n", _n_generations);
  printf("        Inactive Cycles:      %d\n", _n_inactive);
  printf("        Active Cycles:        %d\n\n", _n_generations - _n_inactive);

  initializeOutput();

  printf("----------------------------------------------------------\n");
  printf("| Generation | Shannon Entropy |    keff    |  std_dev   |\n");
  printf("----------------------------------------------------------\n");

  // this is where the simulation will be run

  for (auto i = 0; i < _n_generations; i++)
  {
    std::vector<unsigned long int> source_bins(_n_total_cells, 0);

    unsigned int fissions_in_old_bank = _old_fission_bank.size();

    _n_neutrons_born = 0; // initialize to 0

    for (auto j = 0; j < _n_particles; j++)
    {
      // generate neutrons from fission bank positions first
      if (fissions_in_old_bank > 0 && j < fissions_in_old_bank)
      {
        Neutron neutron = _old_fission_bank[j];
        runHistory(neutron, i, source_bins);
      }
      // if no more neutrons in fission bank, create new ones in fuel
      else
      {
        Cell start_cell = randomFissionCell();
        double random_start_pos = Cell::randomPositionInCell(start_cell, _rng);
        double start_mu = Neutron::randomIsoAngle(_rng);
        Neutron neutron(random_start_pos, start_mu, &start_cell, 1.0);
        runHistory(neutron, i, source_bins);
      }
    }

    while (_split_bank.size() > 0)
    {
      Neutron neutron = _split_bank.front();
      runHistory(neutron, i, source_bins);
      _split_bank.pop_front();
    }

    _shannon_entropy = shannonEntropy(source_bins);

    // need to calculate both generational k and simulation k
    if (i >= _n_inactive)
      calculateK();

    // if too many fission sites are stored, remove the old ones
    while (_new_fission_bank.size() > _n_particles)
    {
      _new_fission_bank.pop_front();
    }
    std::mt19937 gen;
    std::shuffle(_new_fission_bank.begin(), _new_fission_bank.end(), gen);

    assert(_new_fission_bank.size() <= _n_particles);

    _old_fission_bank = _new_fission_bank; // reassign fission bank

    // spit out results
    if (i < _n_inactive)
    {
      // k-eff not calculated for inactive cycles
      printf("|    %3d     |    %.4e   |            |            |\n", i + 1, _shannon_entropy);
    }
    else if (i == _n_inactive)
      printf("|    %3d     |    %.4e   |  %.6f  |            |\n", i + 1, _shannon_entropy, _k_eff);
    else
    {
      // only print std-dev after 3 approximations have been made
      printf("|    %3d     |    %.4e   |  %.6f  |  %.6f  |\n",
             i + 1,
             _shannon_entropy,
             _k_eff,
             _k_std);
    }
  }

  printf("----------------------------------------------------------\n");
  printf("Final k-eff = %.6f +/- %.6f\n\n", _k_eff, _k_std);
  printf("Simulation complete. :-)\n\n");
}

Cell
MCSlab::randomFissionCell()
{
  unsigned int n_fissionable_cells = _fissionable_cells.size();
  const double cell_prob_width = 1.0 / static_cast<double>(n_fissionable_cells);

  auto rn = _rng.generateRN();
  rn = std::min(rn, 1.0 - 1e-15); // safety factor

  unsigned int cell_id = std::floor(rn * static_cast<double>(n_fissionable_cells));
  return _fissionable_cells[cell_id];
}

void
MCSlab::runHistory(Neutron & neutron,
                   const unsigned int generation_num,
                   std::vector<unsigned long int> & source_bins)
{
  neutron.checkNeutron();
  while (neutron.isAlive())
  {
    double distanceToCollision =
        neutron.distanceToCollision();                // store in a variable because this
                                                      // calculation involves a random number
    double distanceToEdge = neutron.distanceToEdge(); // find distance to nearest edge

    if (distanceToEdge < distanceToCollision) // neutron has reached edge of region
    {
      neutron.checkNeutron();
      neutronEscapesCell(neutron, generation_num);
    }
    else // neutron collides in cell
    {
      double collision_location =
          neutron.pos() + distanceToCollision * neutron.mu(); // calculate where collision occurred

      // tally collision position and path traveled
      recordPathLenTally(
          generation_num, neutron.pos(), collision_location, neutron.mu(), neutron.weight());
      recordCollisionTally(generation_num, collision_location, neutron.weight());

      // shift neutron position
      neutron.movePositionWithinCell(collision_location);

      if (_implicit_capture)
        implicitCapture(neutron);
      else
      {
        bool absorbed = testAbsorption(neutron); // did an absorption occur?
        if (absorbed)
        {
          source_bins[collisionIndex(neutron)] += 1;
          absorption(neutron);
        }
        else
          scatter(neutron);
      }
    }
  }
}

bool
MCSlab::testAbsorption(const Neutron & neutron)
{
  double rn = _rng.generateRN();
  return (rn < neutron.cell().absorptionRatio());
}

void
MCSlab::absorption(Neutron & neutron)
{

  unsigned int n_born = nNeutronsBorn(neutron);
  addFissionsToBank(n_born, neutron);
  _n_neutrons_born += n_born;

  neutron.kill(); // kill current neutron
}

void
MCSlab::addFissionsToBank(const unsigned int n_neutrons_born, const Neutron & neutron)
{
  if (n_neutrons_born == 0)
    return;

  const double weight_lost_in_collision =
      !_implicit_capture ? 1.0 : neutron.weight() * neutron.cell().absorptionProbability();

  for (auto i = 0; i < n_neutrons_born; i++)
  {
    double fission_neutron_mu = Neutron::randomIsoAngle(_rng);
    Cell * fission_cell = &_cells[Cell::cellIndex(neutron.pos(), fission_neutron_mu, _cells)];
    Neutron fission_neutron = Neutron(neutron.pos(), fission_neutron_mu, fission_cell, 1.0);
    _new_fission_bank.push_back(fission_neutron); // add to fission bank
  }
}

unsigned int
MCSlab::nNeutronsBorn(const Neutron & neutron)
{
  // skip over unnecessary steps
  if (neutron.cell().nPerAbsorption() < 1e-15)
    return static_cast<unsigned int>(0); // not sure if the cast is necessary...

  // sample number of fission
  double production_rn = _rng.generateRN(); // generate random number
  unsigned int n_born = 0;                  // initialize number of neutrons born
  double neutrons_expected =
      neutron.cell().nPerAbsorption() * neutron.weight(); // expected neutrons from collision

  if (production_rn < neutrons_expected - std::floor(neutrons_expected))
    n_born = std::ceil(neutrons_expected);
  else
    n_born = std::floor(neutrons_expected);

  return static_cast<unsigned int>(n_born);
}

void
MCSlab::scatter(Neutron & neutron)
{
  // this would be where we could calculate the energy lost in scatter.
  // however, because the material properties are energy-indendent, we only
  // need the new angle to proceed.
  double new_mu = neutron.randomIsoAngle(_rng);
  neutron.setMu(new_mu);
}

void
MCSlab::neutronEscapesCell(Neutron & neutron, const unsigned int generation)
{
  // record starting position
  const double start_x = neutron.pos();
  const double mu = neutron.mu();
  const Cell start_cell = neutron.cell();
  const unsigned int start_cell_index = Cell::cellIndex(start_x, mu, _cells);

  const int dir = (mu > 0) ? +1 : -1;
  double x_edge = (dir > 0) ? start_cell.xMax() : start_cell.xMin();     //  put on edge
  recordPathLenTally(generation, start_x, x_edge, mu, neutron.weight()); // tally

  //  if at end of domain, kill
  if ((dir < 0 && start_cell_index == 0) || (dir > 0 && start_cell_index == _cells.size() - 1))
  {
    neutron.kill();
    return;
  }
  neutron.movePositionAndCell(x_edge, _cells);
  if (_use_wws)
  {
    if (!neutron.weightIsOkay())
      splitOrRoulette(neutron);
  }
  neutron.checkNeutron();
}

void
MCSlab::splitOrRoulette(Neutron & neutron)
{
  if (neutron.weight() < neutron.cell().lowerWeight())
    neutron.roulette();
  else
    neutron.split(_split_bank);
}

void
MCSlab::implicitCapture(Neutron & neutron)
{
  assert(_implicit_capture == true);
  const unsigned int n_born = nNeutronsBorn(neutron);
  addFissionsToBank(n_born, neutron);
  _n_neutrons_born += n_born;

  double new_weight = neutron.weight() * (1.0 - neutron.cell().absorptionProbability());
  neutron.changeWeight(new_weight);
  scatter(neutron);
}

unsigned int
MCSlab::collisionIndex(const Neutron & neutron)
{
  double collision_location = neutron.pos();
  for (auto i = 0; i < _n_total_cells; i++)
  {
    if (collision_location < _cells[i].xMax())
      return i;
  }
  throw std::runtime_error("Collision location not within domain of problem!");
}

double
MCSlab::shannonEntropy(const std::vector<unsigned long int> & collision_bins)
{
  double shannon_entropy = 0; // initialize

  // calculate total number of collisions
  uint64_t n_total_collisions = 0ULL;
  n_total_collisions =
      std::accumulate(collision_bins.begin(), collision_bins.end(), static_cast<uint64_t>(0));

  if (n_total_collisions == 0ULL)
    return 0.0; // no collisions

  // normalize collision_bins
  std::vector<double> normalized_col_bins(collision_bins.size(), 0);
  for (auto i = 0; i < collision_bins.size(); i++)
    normalized_col_bins[i] =
        static_cast<double>(collision_bins[i]) / static_cast<double>(n_total_collisions);

  // calculate shannon entropy
  for (auto collision_frac : normalized_col_bins)
  {
    if (collision_frac > 0)
      shannon_entropy -= collision_frac * std::log2(collision_frac);
  }

  return shannon_entropy;
}

void
MCSlab::calculateK()
{
  _k_gen = static_cast<double>(_n_neutrons_born) /
           static_cast<double>(_n_particles); // calculate multiplication factor

  _k_gen_vec.push_back(_k_gen); // add value to running list

  double running_sum = 0;
  for (auto k : _k_gen_vec)
    running_sum += k;
  _k_eff = running_sum / static_cast<double>(_k_gen_vec.size()); // update simulation k-eff

  double sum_sqrd_error = 0;
  for (auto k : _k_gen_vec)
  {
    double error = k - _k_eff;
    sum_sqrd_error += std::pow(error, 2); // sum of squared errors
  }
  _k_std = std::sqrt(sum_sqrd_error /
                     static_cast<double>(_k_gen_vec.size() *
                                         (_k_gen_vec.size() -
                                          1))); // calcuate SAMPLE standard deviation of the mean
}

void
MCSlab::recordCollisionTally(const int current_generation,
                             const double location,
                             const double weight)
{
  if (current_generation > _n_inactive - 1)
  {
    _collision_outfile << location << "," << weight << std::endl;
  }
}

void
MCSlab::recordPathLenTally(const int current_generation,
                           const double start_pos,
                           const double end_pos,
                           const double mu,
                           const double weight)
{
  if (current_generation > _n_inactive - 1)
  {
    double pathlength = (end_pos - start_pos) / mu; // check this!
    _pathlength_outfile << start_pos << "," << end_pos << "," << mu << "," << pathlength << ","
                        << weight << std::endl;
  }
}

void
MCSlab::exportRegionsToCsv(const std::string & outfile)
{
  std::string region_outfile_name = outfile + "_regions.csv";
  std::ofstream region_outfile;
  region_outfile.open(region_outfile_name);
  if (!region_outfile.is_open())
    throw std::runtime_error("Region output file not opened successfully!");

  region_outfile << "region_num,xmin,xmax,Sigma_a,Sigma_s,Sigma_t,nuSigma_f" << std::endl;
  for (auto & region : _regions)
    region_outfile << region.regionIndex() << "," << region.xMin() << "," << region.xMax() << ","
                   << region.SigmaA() << "," << region.SigmaS() << "," << region.SigmaT() << ","
                   << region.nuSigF() << std::endl;

  printf("Regions exported to %s\n\n", region_outfile_name.c_str());
}