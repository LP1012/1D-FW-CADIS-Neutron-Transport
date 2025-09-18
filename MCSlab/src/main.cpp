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
  std::cout << "Hello World :-)" << std::endl;
  Neutron test_neutron(Point(1, 1));
  std::cout << "mu = " << test_neutron.mu() << std::endl;
  std::cout << "x-loc = " << test_neutron.pos().getX() << std::endl;
  std::cout << "y-loc = " << test_neutron.pos().getY() << std::endl;
  return 0;
}