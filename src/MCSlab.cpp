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
               const std::vector<Cell> cells)
  : _n_particles(n_particles),
    _n_generations(n_generations),
    _n_inactive(n_inactive),
    _regions(regions),
    _cells(cells),
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

  _collision_outfile << "position,region,weight,type" << std::endl;
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

  printf("--------------------------------------------------------\n");
  printf("|Generation| Shannon Entropy |    keff    |  std_dev   |\n");
  printf("--------------------------------------------------------\n");

  // this is where the simulation will be run

  for (auto i = 0; i < _n_generations; i++)
  {
    std::vector<unsigned long int> source_bins(_n_total_cells, 0);

    unsigned int fissions_in_old_bank = _old_fission_bank.size();

    _n_neutrons_born = 0; // initialize to 0

    for (auto j = 0; j < _n_particles; j++)
    {
      // generate neutrons from fission bank positions first
      double safe_start_pos =
          (_fissionable_cells.front().xMax() + _fissionable_cells.front().xMin()) / 2.0;
      double starting_mu = Neutron::randomIsoAngle(_rng);
      unsigned int starting_neutron_cell_indedx =
          Cell::cellIndex(safe_start_pos, starting_mu, _cells);
      Neutron neutron(
          safe_start_pos, starting_mu, _cells[starting_neutron_cell_indedx]); // initialize neutron

      // adjust neutron start position based on randomness or fission bank
      if (fissions_in_old_bank > 0 && j < fissions_in_old_bank)
        neutron.movePositionAndCell(_old_fission_bank[j].pos(), _cells);
      else
      {
        neutron.setRandomStartPosition(_fissionable_cells,
                                       _cells); // set location in fuel
      }

      if (!(neutron.pos() <= neutron.cell().xMax() && neutron.pos() >= neutron.cell().xMin()))
        throw std::runtime_error("Neutron cell not set correctly! Position "
                                 "not located within bounds!");

      // begin random walk
      while (neutron.isAlive())
      {
        double distanceToCollision =
            neutron.distanceToCollision();                // store in a variable because this
                                                          // calculation involves a random number
        double distanceToEdge = neutron.distanceToEdge(); // find distance to nearest edge

        if (distanceToEdge < distanceToCollision) // neutron has reached edge of region
          neutronEscapesCell(neutron, i);
        else // neutron collides in region
        {
          double collision_location =
              neutron.pos() +
              distanceToCollision * neutron.mu();  // calculate where collision occurred
          bool absorbed = testAbsorption(neutron); // did an absorption occur?

          // tally collision position and path traveled
          recordPathLenTally(i, neutron.pos(), collision_location, neutron.mu(), neutron.weight());
          recordCollisionTally(i, collision_location, neutron.weight(), absorbed);

          // shift neutron position
          neutron.movePositionWithinRegion(collision_location);

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
    _shannon_entropy = shannonEntropy(source_bins);

    // need to calculate both generational k and simulation k
    if (i >= _n_inactive)
      calculateK();

    // if too many fission sites are stored, remove the old ones
    while (_new_fission_bank.size() > _n_particles)
    {
      _new_fission_bank.pop_front();
    }
    // if (_new_fission_bank.size() > _n_particles)
    // {
    //   _new_fission_bank.erase(_new_fission_bank.begin(),
    //                           _new_fission_bank.begin() +
    //                               (_new_fission_bank.size() - _n_particles));
    // }

    assert(_new_fission_bank.size() <= _n_particles);

    _old_fission_bank = _new_fission_bank; // reassign fission bank

    // spit out results
    if (i < _n_inactive)
    {
      // k-eff not calculated for inactive cycles
      printf("|    %d     |    %.4e   |            |            |\n", i + 1, _shannon_entropy);
    }
    else if (i == _n_inactive)
      printf("|    %d     |    %.4e   |  %.6f  |            |\n", i + 1, _shannon_entropy, _k_eff);
    else
    {
      // only print std-dev after 3 approximations have been made
      printf(
          "|    %d     |    %.4e   |  %.6f  |  %.6f  |\n", i + 1, _shannon_entropy, _k_eff, _k_std);
    }
  }

  printf("--------------------------------------------------------\n");
  printf("Final k-eff = %.6f +/- %.6f\n\n", _k_eff, _k_std);
  printf("Simulation complete. :-)\n\n");
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
  // sample number of fission
  double production_rn = _rng.generateRN(); // generate random number
  Cell neutron_cell = neutron.cell();       // hold region absorption occurs

  unsigned int n_born; // initialize number of neutrons born

  if (production_rn < neutron_cell.nPerAbsorption() - std::floor(neutron_cell.nPerAbsorption()))
    n_born = std::ceil(neutron_cell.nPerAbsorption());
  else
    n_born = std::floor(neutron_cell.nPerAbsorption());

  _n_neutrons_born += n_born;

  for (auto i = 0; i < n_born; i++)
  {
    double fission_neutron_mu = Neutron::randomIsoAngle(_rng);
    Cell fission_cell = _cells[Cell::cellIndex(neutron.pos(), fission_neutron_mu, _cells)];
    Neutron fission_neutron =
        Neutron(neutron.pos(),
                fission_neutron_mu,
                fission_cell);                    // create neutron at location with isotropic angle
    _new_fission_bank.push_back(fission_neutron); // add to fission bank
  }

  neutron.kill(); // kill current neutron
}

void
MCSlab::scatter(Neutron & neutron)
{
  // this would be where we could calculate the energy lost in scatter.
  // however, because the material properties are energy-indendent, we only
  // need the new angle to proceed.
  neutron.randomIsoAngle(_rng);
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
  // //  if next region(s) is/are voids, jump over them
  // unsigned int counting_index = start_cell_index + dir;
  // while (_regions[counting_index].SigmaT() < 1e-15)
  // {
  //   if ((dir < 0 && start_cell_index == 0) ||
  //       (dir > 0 && start_cell_index == _regions.size() - 1))
  //   {
  //     double x_jump = (dir > 0) ? _regions[counting_index].xMax() :
  //     _regions[counting_index].xMin(); recordPathLenTally(generation, start_x, x_jump, mu,
  //     neutron.weight()); neutron.kill(); return;
  //   }
  //   counting_index += dir;
  // }

  // double x_jump = (dir > 0) ? _regions[counting_index].xMin() : _regions[counting_index].xMax();
  // recordPathLenTally(generation, start_x, x_jump, mu, neutron.weight());
  // neutron.movePositionAndRegion(x_jump, _regions);
}

// void
// MCSlab::neutronEscapesRegion(Neutron & neutron, const unsigned int
// generation)
// {
//   unsigned int start_index = neutron.region().regionIndex();
//   unsigned int idx = start_index;

//   const double mu = neutron.mu();
//   const int dir = (mu > 0) ? +1 : -1;

//   // First: record the starting position
//   const double x0 = neutron.pos();

//   double x_edge = (dir > 0) ? _regions[start_index].xMax() :
//   _regions[start_index].xMin();

//   idx += dir;

//   // Skip all void regions (defined by SigmaT == 0)
//   while (_regions[idx].SigmaT() < 1e-15)
//   {
//     idx += dir;
//   }

//   // If at outer edge, neutron escapes the problem
//   if ((dir > 0 && idx == _regions.size() - 1) || (dir < 0 && idx == 0))
//   {
//     // Escape position is boundary of the current (void) region
//     double x_escape = (dir > 0) ? _regions[idx].xMax() :
//     _regions[idx].xMin();

//     recordPathLenTally(generation, x0, x_escape, mu);
//     neutron.kill();
//     return;
//   }

//   // Now idx is the first non-void region in flight direction
//   double new_x = (dir > 0) ? _regions[idx].xMax()  // entering from left
//                            : _regions[idx].xMin(); // entering from right

//   recordPathLenTally(generation, x0, new_x, mu);

//   // Move into the new region
//   // printf("We are somehow moving a region...\n");
//   neutron.movePositionAndRegion(new_x, _regions);
// }

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
                             const double weight,
                             const bool absorbed)
{
  if (current_generation > _n_inactive - 1)
  {
    _collision_outfile << location << "," << weight << ",";
    if (absorbed)
      _collision_outfile << "absorption" << std::endl;
    else
      _collision_outfile << "scatter" << std::endl;
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