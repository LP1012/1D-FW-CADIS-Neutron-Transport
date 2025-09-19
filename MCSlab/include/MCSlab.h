#pragma once

#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "tinyxml2.h"
#include <vector>

class MCSlab {
public:
  MCSlab(const std::string input_file_name);
  /// method to run simulation
  void k_eigenvalue();

  // define getter functions
  unsigned int nParticles() { return _n_particles; }
  unsigned int nGenerations() { return _n_generations; }
  unsigned int nInactive() { return _n_inactive; }

protected:
  /// number of neutrons per generation
  unsigned int _n_particles;
  /// number of generations
  unsigned int _n_generations;
  /// number of inactive cycles
  unsigned int _n_inactive;

  /// simulation input file
  // const tinyxml2::XMLDocument &_input;
  const std::string _input_file_name;

  /// vector of regions
  std::vector<Region> _regions;

  /// bank of source sites
  std::vector<Neutron> _fission_bank;

  /// flux at each point in mesh
  std::vector<double> _scalar_flux;

  void readInput();
};
