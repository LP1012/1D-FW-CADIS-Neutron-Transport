#pragma once

#include "Cell.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include "tinyxml2.h"
#include <vector>
#include <deque>
#include <fstream>

class MCSlab
{
  friend class MCSlabTest;

public:
  MCSlab(const std::string input_file_name,
         const unsigned int n_particles,
         const unsigned int n_generations,
         const unsigned int n_inactive,
         const std::vector<Region> regions,
         const std::vector<Cell> cells,
         const bool implicit_capture,
         const bool use_weight_windows);

  /// method to run simulation
  void k_eigenvalue();

  /// define getter functions
  unsigned int nParticles() { return _n_particles; }
  unsigned int nGenerations() { return _n_generations; }
  unsigned int nInactive() { return _n_inactive; }

protected:
  /// number of neutrons per generation
  const unsigned int _n_particles;
  /// number of generations
  const unsigned int _n_generations;
  /// number of inactive cycles
  const unsigned int _n_inactive;

  const bool _implicit_capture;
  const bool _use_wws;
  bool _export_raw_tallies;

  /// @brief throws a random number to see if absorption occurred
  bool testAbsorption(const Neutron & neutron);
  /// @brief carries out events of an absorption, assuming one occurs
  void absorption(Neutron & neutron);
  /// @brief carries out events of a scattering, assuming one occurs
  void scatter(Neutron & neutron);

  /// @brief computes Shannon Entropy for source convergence
  double shannonEntropy(const std::vector<unsigned long int> & collision_bins);

  /// all cells in simulation
  std::vector<Cell> _cells;
  std::vector<Cell> _fissionable_cells;

  /// simulation input file
  const std::string _input_file_name;

  std::ofstream _collision_outfile;  // output file for storing collisions
  std::ofstream _pathlength_outfile; // output file for pathlength locations

  std::vector<Region> _regions;        // vector of regions
  unsigned int _n_fissionable_regions; // number of fissionable regions

  unsigned int _n_total_cells; // number of cells in all regions

  UniformRNG _rng; // initialize RNG

  unsigned int _n_neutrons_born; // the number of fission neutrons born in a generation

  // banks of source sites
  std::deque<Neutron> _old_fission_bank;
  std::deque<Neutron> _new_fission_bank;

  std::deque<Neutron> _split_bank; // bank of neutrons created via splitting

  void calculateK();
  double _k_gen;                  // multiplication constant of the generation
  std::vector<double> _k_gen_vec; // vector containing running k-gen values
  double _k_eff;                  // running k-eff value of simulation
  double _k_std;                  // std_dev of k-eff

  double _shannon_entropy;

  unsigned int collisionIndex(const Neutron & neutron);

  /// @brief method initializes output files for flux tallies (collision and pathlength)
  void initializeOutput();

  /// @brief exports regions in order from left to right for future postprocessing
  void exportRegionsToCsv(const std::string & outfile);

  /// @brief record the location and region of a single collision
  /// @param current_generation
  void
  recordCollisionTally(const int current_generation, const double location, const double weight);

  /// @brief record start and end location of single neutron movement and the region the movement occurred in
  /// @param current_generation
  void recordPathLenTally(const int current_generation,
                          const double start_pos,
                          const double end_pos,
                          const double mu,
                          const double weight);

  void neutronEscapesCell(Neutron & neutron, const unsigned int generation);
  void createFissionCells();
  void runHistory(Neutron & neutron,
                  const unsigned int generation_num,
                  std::vector<unsigned long int> & source_bins);
  Cell randomFissionCell();
  void splitOrRoulette(Neutron & neutron);
  void implicitCapture(Neutron & neutron);
  unsigned int nNeutronsBorn(const Neutron & neutron);
  void addFissionsToBank(const unsigned int n_neutrons_born, const Neutron & neutron);

  void updatePathLengths(const double x_start,
                         const double x_end,
                         const double mu,
                         const double weight,
                         const unsigned int cell_index);
  void updateCollisions(const double weight, const double Sigma_t, const unsigned int cell_index);
  void exportBinnedTallies();
};
