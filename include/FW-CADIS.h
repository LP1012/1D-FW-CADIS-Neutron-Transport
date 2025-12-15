#pragma once

#include <string>
#include "Region.h"
#include "Cell.h"

class FWCADIS
{
public:
  FWCADIS(const std::string input_file_name);

  void runForwardFlux();

private:
  const std::string _input_file_name;
  unsigned int _n_total_cells;
  bool _use_vr;

  std::vector<Region> _regions;
  std::vector<Cell> _cells;

  // Monte Carlo settings:
  unsigned int _n_particles;
  unsigned int _n_generations;
  unsigned int _n_inactive;

  // SN settings
  unsigned int _quadrature_order;

  void readInput();
};