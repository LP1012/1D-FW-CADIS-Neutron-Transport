#pragma once

#include <string>
#include "Region.h"
#include "Cell.h"

class FWCADIS
{
public:
  FWCADIS(const std::string input_file_name);

private:
  const std::string _input_file_name;
  unsigned int _n_total_cells;

  std::vector<Region> _regions;
  std::vector<Cell> _cells;

  // Monte Carlo settings:
  unsigned int _n_particles;
  unsigned int _n_generations;
  unsigned int _n_inactive;

  void readInput();
};