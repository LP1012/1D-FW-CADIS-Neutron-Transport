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
  MCSlab test_sim_2{"../tests/input_files/test_2_region.xml"};
  MCSlab test_sim_split{"../tests/input_files/test_2_region_split.xml"};

  std::vector<MCSlab> test_sims = {test_sim_1, test_sim_2, test_sim_split};

  // define test-only getters
  unsigned int getNTotalCells(MCSlab sim) const { return sim._n_total_cells; }
  std::vector<Region> regions(MCSlab sim) const { return sim._regions; }
  std::vector<Region> fissionRegions(MCSlab sim) const {
    return sim._fissionable_regions;
  }
  unsigned int nFissionableRegions(MCSlab sim) const {
    return sim._n_fissionable_regions;
  }
};

TEST_F(MCSlabTest, InitializeSimulation) {
  for (auto sim : test_sims) {
    EXPECT_EQ(sim.nParticles(), 100);
    EXPECT_EQ(sim.nGenerations(), 10);
    EXPECT_EQ(sim.nInactive(), 2);
  }
  EXPECT_EQ(getNTotalCells(test_sim_1), 10);
  EXPECT_EQ(getNTotalCells(test_sim_2), 16);
  EXPECT_EQ(getNTotalCells(test_sim_split), 26);

  for (auto i = 0; i < test_sims.size(); i++) {
    EXPECT_EQ(regions(test_sims[i]).size(), i + 1);
  }

  EXPECT_EQ(nFissionableRegions(test_sim_1), 1);
  EXPECT_EQ(nFissionableRegions(test_sim_2), 1);
  EXPECT_EQ(nFissionableRegions(test_sim_split), 1);
}