#include "FW-CADIS.h"
#include "Cell.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "XMLUtils.h"
#include "tinyxml2.h"
#include "utils.h"
#include "SN.h"

#include <vector>
#include <fstream>
#include <iostream>

FWCADIS::FWCADIS(const std::string input_file_name) : _input_file_name(input_file_name)
{
  _n_total_cells = 0;
  readInput();
}

void
FWCADIS::runForwardFlux()
{
  SN simulation{_cells, _quadrature_order};
  simulation.run();
  _forward_flux = simulation.getScalarFlux();
  _forward_k_eff = simulation.k();
  updateForwardFlux(simulation);
}

void
FWCADIS::runAdjointFlux()
{
  SN simulation{_cells, _quadrature_order, true, _forward_k_eff};
  simulation.run();
  updateAdjointFlux(simulation);
}

void
FWCADIS::updateForwardFlux(const SN & simulation)
{
  std::vector<double> forward_flux = simulation.getScalarFlux();
  for (auto i = 0; i < _cells.size(); i++)
    _cells[i].setForwardFlux(forward_flux[i]);
}

void
FWCADIS::updateAdjointFlux(const SN & simulation)
{
  std::vector<double> adjoint_flux = simulation.getScalarFlux();
  for (auto i = 0; i < _cells.size(); i++)
    _cells[i].setAdjointFlux(adjoint_flux[i]);
}

void
FWCADIS::readInput()
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
    auto vol_source = getAttributeOrThrow<double>(region, "source");

    Region region_obj(
        xmin, xmax, n_cells, Sigma_a, Sigma_s, nuSigma_f, vol_source); // create region

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

  // add cells to global cell vector from regions
  for (auto & region : _regions)
  {
    _n_total_cells += region.nCells();
    for (auto & cell : region.cells())
      _cells.push_back(cell);
  }

  // load settings
  auto * settings = root->FirstChildElement("settings");
  if (!settings)
    throw std::runtime_error("<settings> element not set!");
  // test for existence of data, then set private attributes
  auto * n_part_attrib = settings->FindAttribute("n_particles");
  auto * n_gen_attrib = settings->FindAttribute("n_generations");
  auto * n_inactive_attrib = settings->FindAttribute("n_inactive");
  auto * use_vr = settings->FindAttribute("fw-cadis");

  _n_particles = getAttributeOrThrow<unsigned int>(settings, "n_particles");
  _n_generations = getAttributeOrThrow<unsigned int>(settings, "n_generations");
  _n_inactive = getAttributeOrThrow<unsigned int>(settings, "n_inactive");

  if (!use_vr)
    throw std::runtime_error("fw-cadis parameter not set to true or false in input file!");
  _use_vr = use_vr->BoolValue();
  if (use_vr)
  {
    auto * quad_order = settings->FindAttribute("quadrature-order");
    if (!quad_order)
      throw std::runtime_error("quadrature-order not set!");
    auto * www = settings->FindAttribute("window-width");
    if (!www)
      throw std::runtime_error("window-width not set!");
    _quadrature_order = quad_order->UnsignedValue();
    _window_width = www->DoubleValue();
  }
}