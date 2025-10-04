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
  MCSlab test_sim_1{"../tests/input_files/test_1_region.xml"};
  // MCSlab test_sim_2{"../tests/input_files/test_2_region.xml"};
  // MCSlab test_sim_split{"../tests/input_files/test_2_region_split.xml"};

  // define test-only getters
  unsigned int getNTotalCells() { return test_sim_1._n_total_cells; }
};

TEST_F(MCSlabTest, InitializeSimulation) {
  EXPECT_EQ(test_sim_1.nParticles(), 100);
}