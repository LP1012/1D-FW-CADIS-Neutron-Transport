// include files
#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

int main() {
  MCSlab test_mcslab("test_1_region.xml");
  std::cout << "n_inactive = " << test_mcslab.nInactive() << std::endl;
  return 0;
}