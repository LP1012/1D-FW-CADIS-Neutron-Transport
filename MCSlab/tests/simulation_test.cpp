#include "MCSlab.h"
#include "Neutron.h"
#include "Region.h"
#include <gtest/gtest.h>
#include <iostream> // delete later
#include <optional>
#include <vector>

#include <iomanip> // delete later

class SimTest : public testing::Test {
protected:
  std::vector<Region> regions{Region(-1, 1, 1, 1, 1, 0),
                              Region(1, 2, 1, 1, 1, 1)};
  Neutron neutron{0, regions, 0};
};