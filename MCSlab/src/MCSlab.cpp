#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "tinyxml2.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

MCSlab::MCSlab(const std::string input_file_name)
    : _input_file_name(input_file_name) {
  readInput();
};

void MCSlab::k_eigenvalue() {
  // this is where the simulation will be run
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
    auto *id = region->FindAttribute("id");
    if (!id)
      throw std::runtime_error("id not set for region!");

    auto *xmin_attrib = region->FindAttribute("xmin");
    if (!xmin_attrib)
      throw std::runtime_error("xmin not set for region!");
    auto xmin = xmin_attrib->DoubleValue();

    auto *xmax_attrib = region->FindAttribute("xmax");
    if (!xmax_attrib)
      throw std::runtime_error("xmax not set for region!");
    auto xmax = xmax_attrib->DoubleValue();

    auto *n_cells_attrib = region->FindAttribute("n_cells");
    if (!n_cells_attrib)
      throw std::runtime_error("n_cells not set for region!");
    auto n_cells = n_cells_attrib->UnsignedValue();

    auto *Sigma_a_attrib = region->FindAttribute("Sigma_a");
    if (!Sigma_a_attrib)
      throw std::runtime_error("Sigma_a not set for region!");
    auto Sigma_a = Sigma_a_attrib->DoubleValue();

    auto *Sigma_s_attrib = region->FindAttribute("Sigma_s");
    if (!Sigma_s_attrib)
      throw std::runtime_error("Sigma_s not set for region!");
    auto Sigma_s = Sigma_s_attrib->DoubleValue();

    auto *nuSigma_f_attrib = region->FindAttribute("nuSigma_f");
    if (!nuSigma_f_attrib)
      throw std::runtime_error("nuSigma_f not set for region!");
    auto nuSigma_f = nuSigma_f_attrib->DoubleValue();

    Region region_obj(xmin, xmax, n_cells, Sigma_a, Sigma_s,
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

  if (!n_part_attrib)
    throw std::runtime_error("n_particles not set!");
  if (!n_gen_attrib)
    throw std::runtime_error("n_generations not set!");
  if (!n_inactive_attrib)
    throw std::runtime_error("n_inactive not set!");

  _n_particles = n_part_attrib->UnsignedValue();
  _n_generations = n_gen_attrib->UnsignedValue();
  _n_inactive = n_inactive_attrib->UnsignedValue();
}
