// include files
#include "FW-CADIS.h"
#include "MCSlab.h"
#include "Neutron.h"
#include "Point.h"
#include "Region.h"

#include <cmath>
#include <iostream>
#include <random>
#include <vector>

int
main(int argc, char * argv[])
{
  if (argc == 1)
    throw std::runtime_error("Program must have input file command line argument!");
  else if (argc > 2)
    throw std::runtime_error("Program only takes one argument!");

  std::string input_filename = argv[1];
  FWCADIS simulation{input_filename};
  // simulation.k_eigenvalue();
  return 0;
}