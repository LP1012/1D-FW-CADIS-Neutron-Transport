#pragma once

#include <string>
#include "Region.h"
#include "Cell.h"

class FWCADIS
{
  FWCADIS(const std::string input_file_name);

public:
private:
  const std::string _input_file_name;
  std::vector<Cell> _cells;
  unsigned int _n_total_cells;

  std::vector<Region> _regions;
  std::vector<Cell> _cells;

  // Monte Carlo settings:
  unsigned int _n_particles;
  unsigned int _n_generations;
  unsigned int _n_inactive;

  void readInput();
};