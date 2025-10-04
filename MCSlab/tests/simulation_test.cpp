#include "MCSlab.h"
#include "Neutron.h"
#include "Region.h"
#include <gtest/gtest.h>
#include <iostream> // delete later
#include <optional>
#include <vector>

#include <iomanip> // delete later

class MCSlabTest : public testing::Test {
protected:
  MCSlab test_sim{"../tests/input_files/test_2_region_split.xml"};
};

TEST_F(MCSlabTest, InitializeSimulation) {}