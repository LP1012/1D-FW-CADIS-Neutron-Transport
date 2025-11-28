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

MCSlab::MCSlab(const std::string input_file_name) : _input_file_name(input_file_name), _rng()
{
  _n_total_cells = 0;
  _n_fissionable_regions = 0;
  readInput();
  fissionRegions();
};

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

  _collision_outfile << "position,region,type" << std::endl;
  _pathlength_outfile << "start,end,mu,pathlength,region" << std::endl;

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
          (_fissionable_regions.front().xMax() + _fissionable_regions.front().xMin()) / 2.0;
      Neutron neutron(safe_start_pos, _regions); // initialize neutron

      // adjust neutron start position based on randomness or fission bank
      if (fissions_in_old_bank > 0 && j < fissions_in_old_bank - 1)
        neutron.movePositionAndRegion(_old_fission_bank[j].pos(), _regions);
      else
      {

        neutron.setRandomStartPosition(_fissionable_regions,
                                       _regions); // set location in fuel
      }

      if (!(neutron.pos() <= neutron.region().xMax() && neutron.pos() >= neutron.region().xMin()))
        throw std::runtime_error(
            "Neutron region not set correctly! Position not located within bounds!");

      // begin random walk
      while (neutron.isAlive())
      {
        double distanceToCollision =
            neutron.distanceToCollision(); // store in a variable because this calculation involves
                                           // a random number
        double distanceToEdge = neutron.distanceToEdge(); // find distance to nearest edge

        if (distanceToEdge < distanceToCollision) // neutron has reached edge of region
          neutronEscapesRegion(neutron, i);
        else // neutron collides in region
        {
          double collision_location =
              neutron.pos() +
              distanceToCollision * neutron.mu();  // calculate where collision occurred
          bool absorbed = testAbsorption(neutron); // did an absorption occur?

          // tally collision position and path traveled
          recordPathLenTally(
              i, neutron.pos(), collision_location, neutron.mu(), neutron.region().regionIndex());
          recordCollisionTally(i, collision_location, neutron.region().regionIndex(), absorbed);

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
      _new_fission_bank.front() = _new_fission_bank.back();
      _new_fission_bank.pop_back();
    }
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
  return (rn < neutron.region().absorptionRatio());
}

void
MCSlab::absorption(Neutron & neutron)
{
  // sample number of fission
  double production_rn = _rng.generateRN(); // generate random number
  Region neutron_region = neutron.region(); // hold region absorption occurs

  unsigned int n_born; // initialize number of neutrons born

  if (production_rn < neutron_region.nPerAbsorption() - std::floor(neutron_region.nPerAbsorption()))
    n_born = std::ceil(neutron_region.nPerAbsorption());
  else
    n_born = std::floor(neutron_region.nPerAbsorption());

  _n_neutrons_born += n_born;

  for (auto i = 0; i < n_born; i++)
  {
    Neutron fission_neutron = Neutron(neutron.pos(),
                                      _regions);  // create neutron at location with isotropic angle
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
  neutron.randomIsoAngle();
}

void
MCSlab::neutronEscapesRegion(Neutron & neutron, const unsigned int generation)
{
  unsigned int start_index = neutron.region().regionIndex();
  unsigned int idx = start_index;

  const double mu = neutron.mu();
  const int dir = (mu > 0) ? +1 : -1;

  // First: record the starting position
  const double x0 = neutron.pos();

  // If at outer edge, neutron escapes the problem
  if ((dir > 0 && idx == _regions.size() - 1) || (dir < 0 && idx == 0))
  {
    // Escape position is boundary of the current (void) region
    double x_escape = (dir > 0) ? _regions[idx].xMax() : _regions[idx].xMin();

    recordPathLenTally(generation, x0, x_escape, mu, start_index);
    neutron.kill();
    return;
  }

  // Skip all void regions (defined by SigmaT == 0)
  while (_regions[idx].SigmaT() < 1e-15)
  {
    // If at outer edge, neutron escapes the problem
    if ((dir > 0 && idx == _regions.size() - 1) || (dir < 0 && idx == 0))
    {
      // Escape position is boundary of the current (void) region
      double x_escape = (dir > 0) ? _regions[idx].xMax() : _regions[idx].xMin();

      recordPathLenTally(generation, x0, x_escape, mu, start_index);
      neutron.kill();
      return;
    }
    idx += dir;
  }

  // Now idx is the first non-void region in flight direction
  double new_x = (dir > 0) ? _regions[idx].xMax()  // entering from left
                           : _regions[idx].xMin(); // entering from right

  recordPathLenTally(generation, x0, new_x, mu, start_index);

  // Move into the new region
  // printf("We are somehow moving a region...\n");
  neutron.movePositionAndRegion(new_x, _regions);
}

void
MCSlab::readInput()
{
  tinyxml2::XMLDocument input_file;

  // check if the input file exists
  if (input_file.LoadFile(_input_file_name.c_str()) != tinyxml2::XML_SUCCESS)
    throw std::runtime_error("Input file not found!");

  auto * root = input_file.FirstChildElement("simulation");
  if (!root)
    throw std::runtime_error("No <simulation> root element");

  // loop over regions
  auto * regionsElement = root->FirstChildElement("regions");
  auto * region = regionsElement->FirstChildElement("region");
  while (region)
  {
    // auto id = getAttributeOrThrow<unsigned int>(region, "id");
    auto xmin = getAttributeOrThrow<double>(region, "xmin");
    auto xmax = getAttributeOrThrow<double>(region, "xmax");
    auto n_cells = getAttributeOrThrow<unsigned int>(region, "n_cells");
    auto Sigma_a = getAttributeOrThrow<double>(region, "Sigma_a");
    auto Sigma_s = getAttributeOrThrow<double>(region, "Sigma_s");
    auto nuSigma_f = getAttributeOrThrow<double>(region, "nuSigma_f");

    Region region_obj(xmin, xmax, n_cells, Sigma_a, Sigma_s,
                      nuSigma_f); // create region

    // add checks for overlap and void regions here
    if (_regions.size() > 0)
    {
      auto prev_region = _regions.back();

      if (prev_region.xMax() < region_obj.xMin())
      {
        // add a void region between separated regions
        Region void_region = Region::voidRegion(prev_region.xMax(), region_obj.xMin(), 10);
        void_region.setIndex(_regions.size()); // check this
        _regions.push_back(void_region);
      }
      else if (prev_region.xMax() > region_obj.xMin())
        throw std::runtime_error("Error! Regions overlap."); // check if overlap
      else if (prev_region.xMin() > region_obj.xMin())
        throw std::runtime_error("Error! Regions are not sorted"); // check regions are sorted
    }

    region_obj.setIndex(_regions.size());
    _regions.push_back(region_obj); // add region to list of regions

    region = region->NextSiblingElement("region"); // move to next
  }

  // add cell bounds and widths
  unsigned int count = 0;
  for (auto region : _regions)
  {
    _n_total_cells += region.nCells();
    if (count == 0)
    {
      for (auto i = 0; i < region.cellBounds().size(); i++)
      {
        // each region holds its own set of bounds, so we need an inner
        // loop over each region's bounds
        _all_cell_bounds.push_back(region.cellBounds()[i]);
      }
    }
    else
    {
      for (auto i = 1; i < region.cellBounds().size(); i++)
      {
        // each region holds its own set of bounds, so we need an inner
        // loop over each region's bounds
        _all_cell_bounds.push_back(region.cellBounds()[i]);
      }
    }
    for (auto i = 0; i < region.cellLocs().size(); i++)
      _cell_widths.push_back(region.cellLocs()[i][1] - region.cellLocs()[i][0]);

    for (auto center : region.cellCenters())
      _all_cell_centers.push_back(center);

    for (auto i = 0; i < region.nCells(); i++)
      _Sigma_t_vals.push_back(region.SigmaT());

    count++;
  }

  // load settings
  auto * settings = root->FirstChildElement("settings");
  if (!settings)
    throw std::runtime_error("<settings> element not set!");
  // test for existence of data, then set private attributes
  auto * n_part_attrib = settings->FindAttribute("n_particles");
  auto * n_gen_attrib = settings->FindAttribute("n_generations");
  auto * n_inactive_attrib = settings->FindAttribute("n_inactive");

  _n_particles = getAttributeOrThrow<unsigned int>(settings, "n_particles");
  _n_generations = getAttributeOrThrow<unsigned int>(settings, "n_generations");
  _n_inactive = getAttributeOrThrow<unsigned int>(settings, "n_inactive");
}

void
MCSlab::fissionRegions()
{
  for (auto region : _regions)
  {
    if (region.nuSigF() > 1e-15)
      _fissionable_regions.push_back(region);
  }
  _n_fissionable_regions = _fissionable_regions.size();
}

unsigned int
MCSlab::collisionIndex(const Neutron & neutron)
{
  double collision_location = neutron.pos();
  for (auto i = 1; i <= _n_total_cells + 1; i++)
  {
    if (collision_location < _all_cell_bounds[i])
      return i - 1;
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
                             const unsigned int region_num,
                             const bool absorbed)
{
  if (current_generation > _n_inactive - 1)
  {
    _collision_outfile << location << "," << region_num << ",";
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
                           const unsigned int region_num)
{
  if (current_generation > _n_inactive - 1)
  {
    double pathlength = (end_pos - start_pos) / mu; // check this!
    _pathlength_outfile << start_pos << "," << end_pos << "," << mu << "," << pathlength << ","
                        << region_num << std::endl;
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