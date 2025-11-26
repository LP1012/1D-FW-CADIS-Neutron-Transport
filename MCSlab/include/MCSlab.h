#pragma once

#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "tinyxml2.h"
#include <vector>
#include <fstream>

class MCSlab
{
  friend class MCSlabTest;

public:
  MCSlab(const std::string input_file_name);

  /// method to run simulation
  void k_eigenvalue();

  /// @brief throws a random number to see if absorption occurred
  bool testAbsorption(const Neutron & neutron);
  /// @brief carries out events of an absorption, assuming one occurs
  void absorption(Neutron & neutron);
  /// @brief carries out events of a scattering, assuming one occurs
  void scatter(Neutron & neutron);

  /// @brief computes Shannon Entropy for source convergence
  double shannonEntropy(const std::vector<unsigned long int> & collision_bins);

  /// define getter functions
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
  const std::string _input_file_name;

  std::ofstream _collision_outfile;  // output file for storing collisions
  std::ofstream _pathlength_outfile; // output file for pathlength locations

  std::vector<Region> _regions;             // vector of regions
  std::vector<Region> _fissionable_regions; // vector of fissile regions
  unsigned int _n_fissionable_regions;      // number of fissionable regions
  std::vector<double> _all_cell_bounds;     // all cell boundaries
  unsigned int _n_total_cells;              // number of cells in all regions
  std::vector<double> _cell_widths;         // width of each cell
  std::vector<double> _all_cell_centers;    // vector of all cell center locations
  std::vector<double> _Sigma_t_vals;

  UniformRNG _rng; // initialize RNG

  unsigned int _n_neutrons_born; // the number of fission neutrons born in generation

  // banks of source sites
  std::vector<Neutron> _old_fission_bank;
  std::vector<Neutron> _new_fission_bank;

  void calculateK();
  double _k_gen;                  // multiplication constant of the generation
  std::vector<double> _k_gen_vec; // vector containing running k-gen values
  double _k_eff;                  // running k-eff value of simulation
  double _k_std;                  // std_dev of k-eff

  double _shannon_entropy;

  void readInput();
  void fissionRegions();

  unsigned int collisionIndex(const Neutron & neutron);

  /// @brief method initializes output files for flux tallies (collision and pathlength)
  void initializeOutput();

  /// @brief exports regions in order from left to right for future postprocessing
  void exportRegionsToCsv();

  /// @brief record the location and region of a single collision
  /// @param current_generation
  void recordCollisionTally(const int current_generation,
                            const double location,
                            const unsigned int region_num,
                            const bool absorbed);

  /// @brief record start and end location of single neutron movement and the region the movement occurred in
  /// @param current_generation
  void recordPathLenTally(const int current_generation,
                          const double start_pos,
                          const double end_pos,
                          const double mu,
                          const unsigned int region_num);
};
