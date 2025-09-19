// include files
#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Slab.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

int main() {
  Slab test_slab(1, 2, 4, 1, 1, 0);
  for (auto center : test_slab.cellCenters()) {
    std::cout << "cell center at: " << center << std::endl;
  }

  for (auto loc : test_slab.cellLocs()) {
    std::cout << "cells bounded by: " << loc[0] << ", " << loc[1] << std::endl;
  }
  return 0;
}