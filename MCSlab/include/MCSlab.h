#pragma once

#include "Neutron.h"
#include "Point.h"
#include "tinyxml2.h"
#include <vector>

class MCSlab {
public:
  MCSlab(const tinyxml2::XMLDocument &input_file);
  /// method to run simulation
  void k_eigenvalue();

protected:
  /// number of neutrons per generation
  unsigned int _n_particles;
  /// number of generations
  unsigned int _n_generations;
  /// number of inactive cycles
  unsigned int _n_inactive;

  /// simulation input file
  const tinyxml2::XMLDocument &_input;

  /// bank of source sites
  std::vector<Neutron> _fission_bank;

  /// flux at each point in mesh
  std::vector<double> _scalar_flux;
};
