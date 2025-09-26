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

  void collision(Neutron neutron);

  // define getter functions
  unsigned int nParticles() { return _n_particles; }
  unsigned int nGenerations() { return _n_generations; }
  unsigned int nInactive() { return _n_inactive; }

protected:
  unsigned int _n_particles;   // number of neutrons per generation
  unsigned int _n_generations; // number of generations
  unsigned int _n_inactive;    // number of inactive cycles

  /// simulation input file
  const std::string _input_file_name;

  std::vector<Region> _regions;             // vector of regions
  std::vector<Region> _fissionable_regions; // vector of fissile regions
  unsigned int _n_fissionable_regions;      // number of fissionable regions

  // hold  min and max of computation domain
  double _domainMin;
  double _domainMax;

  // initialize RNG
  UniformRNG _rng;

  // banks of source sites
  std::vector<Neutron> _old_fission_bank;
  std::vector<Neutron> _new_fission_bank;

  double _k; // multiplication constant

  std::vector<double> _scalar_flux; // flux at each point in mesh

  void readInput();
  void fissionRegions();

  void setMinMax();
};
