#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "XMLUtils.h"
#include "tinyxml2.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

MCSlab::MCSlab(const std::string input_file_name)
    : _input_file_name(input_file_name) {
  readInput();
  setMinMax();
};

void MCSlab::k_eigenvalue() {
  // this is where the simulation will be run
  for (auto i = 0; _n_generations; i++) {
    // put something in to not count tallies for (i-1)<n_inactive
    for (auto j = 0; _n_particles; j++) {
      // generate neutrons from fission bank positions first
      // use the length of the bank and compare to j
      Neutron neutron(0);
      neutron.setRandomStartPosition(
          _fissionable_regions); // set location in fuel

      // begin random walk
      double mean_free_path = MCSlab::MFP(neutron.regionID());
      double distanceToCollision = neutron.distanceToCollision(mean_free_path);

      // find distance to nearest edge
      double distanceToEdge = neutron.distanceToEdge();
    }
  }
}

double MCSlab::MFP(const unsigned int id) {
  for (auto region : _regions) {
    if (region.id() == id)
      return region.mfp();
  }
  // add case for MFP in void
}

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
    auto id = getAttributeOrThrow<unsigned int>(region, "id");
    auto xmin = getAttributeOrThrow<double>(region, "xmax");
    auto xmax = getAttributeOrThrow<double>(region, "xmin");
    auto n_cells = getAttributeOrThrow<unsigned int>(region, "n_cells");
    auto Sigma_a = getAttributeOrThrow<double>(region, "Sigma_a");
    auto Sigma_s = getAttributeOrThrow<double>(region, "Sigma_s");
    auto nuSigma_f = getAttributeOrThrow<double>(region, "Sigma_f");

    Region region_obj(id, xmin, xmax, n_cells, Sigma_a, Sigma_s,
                      nuSigma_f);   // create region
    _regions.push_back(region_obj); // add region to list of regions
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
