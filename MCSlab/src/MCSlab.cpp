#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "XMLUtils.h"
#include "tinyxml2.h"

#include <algorithm> // for std::sort
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <vector>

MCSlab::MCSlab(const std::string input_file_name)
    : _input_file_name(input_file_name), _rng() {
  readInput();
  setMinMax();
};

void MCSlab::k_eigenvalue() {
  // this is where the simulation will be run
  for (auto i = 0; i < _n_generations; i++) {
    // define bins
    std::vector<unsigned long int> collision_bins(_n_total_cells, 0);

    // put something in to not count tallies for (i-1)<n_inactive
    unsigned int fissions_in_old_bank = _old_fission_bank.size();
    for (auto j = 0; j < _n_particles; j++) {
      // generate neutrons from fission bank positions first
      Neutron neutron(0, _regions); // initialize neutron

      // adjust neutron start position based on randomness or fission bank
      if (j < fissions_in_old_bank - 1)
        neutron.movePositionAndRegion(_old_fission_bank[j].pos(), _regions);
      else
        neutron.setRandomStartPosition(
            _fissionable_regions); // set location in fuel

      // begin random walk
      while (neutron.isAlive()) {

        double distanceToCollision = neutron.distanceToCollision();

        // find distance to nearest edge
        double distanceToEdge = neutron.distanceToEdge();

        if (distanceToCollision > distanceToEdge) {
          // neutron has reached edge of region
          unsigned int current_index = neutron.region().regionIndex();
          if (neutron.mu() > 0) {
            if (current_index == _regions.size() - 1)
              // neutron escapes on right side
              neutron.kill();
            else {
              // move neutron to a region on the right
              unsigned int new_index = current_index + 1;
              while (_regions[new_index].SigmaA() < 1e-8)
                new_index++; // skip over void regions

              double new_position = _regions[new_index].xMin();
              neutron.movePositionAndRegion(new_position, _regions);
            }
          } else {
            if (current_index == 0)
              // neutron escapes on left side
              neutron.kill();
            else {
              // move neutron to a region on the left
              unsigned int new_index = current_index - 1;
              while (_regions[new_index].SigmaA() < 1e-8)
                new_index--; // skip over void regions

              double new_position = _regions[new_index].xMax();
              neutron.movePositionAndRegion(new_position, _regions);
            }
          }

        } else {
          collision_bins[collisionIndex(neutron)] +=
              1; // add one collision to bin
          MCSlab::collision(neutron);
        }
      }
    }
    _k = static_cast<double>(_new_fission_bank.size()) /
         static_cast<double>(_n_particles); // calculate multiplication factor
    _old_fission_bank = _new_fission_bank;  // reassign fission bank
    _new_fission_bank.clear(); // clear new bank for next generation
  }
}

void MCSlab::collision(Neutron neutron) {
  // test if collision or absorption
  double rn = _rng.generateRN();
  if (rn < neutron.region().absorptionRatio()) {
    // sample number of fission
    double production_rn = _rng.generateRN(); // generate random number
    Region neutron_region = neutron.region(); // hold region absorption occurs

    unsigned int n_born; // initialize number of neutrons born

    (production_rn < neutron_region.nPerAbsorption() -
                         std::floor(neutron_region.nPerAbsorption()))
        ? n_born = std::floor(neutron_region.nPerAbsorption())
        : n_born = std::floor(neutron_region.nPerAbsorption()) +
                   1; // determine number of neutrons born

    for (auto i = 0; i < n_born; i++) {
      Neutron fission_neutron =
          Neutron(neutron.pos(),
                  _regions); // create neutron at location with isotropic angle
      _new_fission_bank.push_back(fission_neutron); // add to fission bank
    }

    neutron.kill(); // kill current neutron
  } else {
    // this would be where we could calculate the energy lost in scatter.
    // however, because the material properties are energy-indendent, we only
    // need the new angle to proceed.
    neutron.randomIsoAngle();
  }
}

// double MCSlab::MFP(const unsigned int id) {
//   for (auto region : _regions) {
//     if (region.id() == id)
//       return region.mfp();
//   }
//   // add case for MFP in void
// }

void MCSlab::readInput() {
  tinyxml2::XMLDocument input_file;

  // check if the input file exists
  if (input_file.LoadFile(_input_file_name.c_str()) != tinyxml2::XML_SUCCESS)
    throw std::runtime_error("Input file not found!");

  auto *root = input_file.FirstChildElement("simulation");
  if (!root)
    throw std::runtime_error("No <simulation> root element");

  // loop over regions
  auto *regionsElement = root->FirstChildElement("regions");
  auto *region = root->FirstChildElement("region");
  while (region) {
    // auto id = getAttributeOrThrow<unsigned int>(region, "id");
    auto xmin = getAttributeOrThrow<double>(region, "xmin");
    auto xmax = getAttributeOrThrow<double>(region, "xmax");
    auto n_cells = getAttributeOrThrow<unsigned int>(region, "n_cells");
    auto Sigma_a = getAttributeOrThrow<double>(region, "Sigma_a");
    auto Sigma_s = getAttributeOrThrow<double>(region, "Sigma_s");
    auto nuSigma_f = getAttributeOrThrow<double>(region, "Sigma_f");

    Region region_obj(xmin, xmax, n_cells, Sigma_a, Sigma_s,
                      nuSigma_f);   // create region
    _regions.push_back(region_obj); // add region to list of regions

    region = region->NextSiblingElement("region"); // move to next
  }

  // add void regions, check if regions overlap
  std::vector<Region> new_regions;
  for (auto i = 0; i < _regions.size(); i++) {

    new_regions.push_back(_regions[i]);           // add user-specified region
    _regions[i].setIndex(new_regions.size() - 1); // set region index

    if (_regions[i].xMax() < _regions[i + 1].xMin()) {
      // create void region between user-specified regions
      Region void_region =
          Region::voidRegion(_regions[i].xMax(), _regions[i + 1].xMin(), 10);
      void_region.setIndex(new_regions.size() + 1); // set void region index
      new_regions.push_back(void_region); // add void region to list of regions

    } else if (_regions[i].xMax() > _regions[i + 1].xMin()) {
      throw std::runtime_error("Error! Regions overlap."); // check if overlap

    } else if (_regions[i].xMin() > _regions[i + 1].xMin())
      throw std::runtime_error("Error! Regions are not sorted");
  }

  // add cell bounds
  for (auto region : _regions) {
    for (auto i = 0; i < region.cellBounds().size(); i++) {
      // each region holds its own set of bounds, so we need an inner
      // loop over each region's bounds
      _all_cell_bounds.push_back(region.cellBounds()[i]);
    }
  }

  // load settings
  auto *settings = root->FirstChildElement("settings");
  if (!settings)
    throw std::runtime_error("<settings> element not set!");
  // test for existence of data, then set private attributes
  auto *n_part_attrib = settings->FindAttribute("n_particles");
  auto *n_gen_attrib = settings->FindAttribute("n_generations");
  auto *n_inactive_attrib = settings->FindAttribute("n_inactive");

  _n_particles = getAttributeOrThrow<unsigned int>(settings, "n_particles");
  _n_generations = getAttributeOrThrow<unsigned int>(settings, "n_generations");
  _n_inactive = getAttributeOrThrow<unsigned int>(settings, "n_inactive");
}

void MCSlab::fissionRegions() {
  for (auto region : _regions) {
    if (region.nuSigF() > 1e-8)
      _fissionable_regions.push_back(region);
  }
  _n_fissionable_regions = _fissionable_regions.size();
}

void MCSlab::setMinMax() {
  _domainMin = _regions[0].xMin();
  _domainMax = _regions[0].xMax();

  for (auto region : _regions) {
    if (region.xMin() < _domainMin)
      _domainMin = region.xMin();
    if (region.xMax() > _domainMax)
      _domainMax = region.xMax();
  }
}

unsigned int MCSlab::collisionIndex(const Neutron &neutron) {
  double collision_location = neutron.pos();
  for (auto i = 1; i < _n_total_cells + 1; i++) {
    if (collision_location < _all_cell_bounds[i])
      return i - 1; // CHECK THIS LOGIC
  }
  throw std::runtime_error("Collision location not within domain of problem!");
}

double
MCSlab::shannonEntropy(const std::vector<unsigned long int> &collision_bins) {
  double shannon_entropy = 0; // initialize

  // calculate total number of collisions
  unsigned int n_total_collisions;
  for (auto collisions : collision_bins)
    n_total_collisions += collisions;

  // normalize collision_bins
  std::vector<double> normalized_col_bins(collision_bins.size(), 0);
  for (auto i = 0; i < collision_bins.size(); i++)
    normalized_col_bins[i] = static_cast<double>(collision_bins[i]) /
                             static_cast<double>(n_total_collisions);

  // calculate shannon entropy
  for (auto collision_frac : normalized_col_bins)
    shannon_entropy -= collision_frac * std::log2(collision_frac);

  return shannon_entropy;
}