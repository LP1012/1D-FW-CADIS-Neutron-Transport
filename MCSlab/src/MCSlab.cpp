#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Slab.h"
#include "tinyxml2.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

MCSlab::MCSlab(const unsigned int x_cells, const unsigned int n_particles,
               const unsigned int n_generations, const unsigned int n_inactive)
    : _x_cells(x_cells), _n_particles(n_particles),
      _n_generations(n_generations), _n_inactive(n_inactive){};

void MCSlab::k_eigenvalue() {
  // this is where the simulation will be run
}