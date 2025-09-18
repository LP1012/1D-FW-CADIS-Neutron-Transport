#pragma once

#include "Neutron.h"
#include "Point.h"
#include "tinyxml2.h"
#include <vector>

class MCSlab {
public:
  MCSlab(const unsigned int x_cells, const unsigned int n_particles,
         const unsigned int n_generations, const unsigned int n_inactive);
  /// method to run simulation
  void k_eigenvalue();

protected:
  /// mesh size for Shannon entropy calcs
  const unsigned int _x_cells;

  /// number of neutrons per generation
  const unsigned int _n_particles;
  /// number of generations
  const unsigned int _n_generations;
  /// number of inactive cycles
  const unsigned int _n_inactive;

  /// simulation input file
  const tinyxml2::XMLDocument _input;

  /// bank of source sites
  std::vector<Neutron> _fission_bank;

  /// flux at each point in mesh
  std::vector<double> _scalar_flux;
};
