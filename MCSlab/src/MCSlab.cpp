#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Slab.h"
#include "tinyxml2.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

MCSlab::MCSlab(const tinyxml2::XMLDocument &input_file) : _input(input_file){};

void MCSlab::k_eigenvalue() {
  // this is where the simulation will be run
}