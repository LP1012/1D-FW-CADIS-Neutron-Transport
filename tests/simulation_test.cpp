// #include "MCSlab.h"
// #include "Cell.h"
// #include "Neutron.h"
// #include "Region.h"
// #include <gtest/gtest.h>
// #include <iostream> // delete later
// #include <optional>
// #include <vector>

// #include <iomanip> // delete later

// class MCSlabTest : public testing::Test
// {
// protected:
//   std::vector<Cell> test_sim_1_cells = {Cell{-1.0, -0.5, 1.0, 2.0, 4},
//                                         Cell{-0.5, 0, 1.0, 2.0, 4},
//                                         Cell{0, 0.5, 1.0, 2.0, 4},
//                                         Cell{0.5, 1.0, 1.0, 2.0, 4}};
//   MCSlab test_sim_1{"test_1_region.xml",100,10,2,};
//   MCSlab test_sim_2{"../tests/input_files/test_2_region.xml"};
//   MCSlab test_sim_split{"../tests/input_files/test_2_region_split.xml"};

//   // std::vector<MCSlab> test_sims = {test_sim_1, test_sim_2, test_sim_split};
//   std::vector<MCSlab *> test_sims{&test_sim_1, &test_sim_2, &test_sim_split};

//   // define test-only getters
//   unsigned int getNTotalCells(const MCSlab & sim) const { return sim._n_total_cells; }
//   std::vector<Region> regions(const MCSlab & sim) const { return sim._regions; }

//   std::vector<Region> fissionRegions(const MCSlab & sim) const { return sim._fissionable_regions;
//   } unsigned int nFissionableRegions(const MCSlab & sim) const { return
//   sim._n_fissionable_regions; }

//   unsigned int lenFissionBank(const MCSlab & sim) { return sim._new_fission_bank.size(); }
// };

// TEST_F(MCSlabTest, InitializeSimulation)
// {
//   for (auto & sim : test_sims)
//   {
//     EXPECT_EQ(sim->nParticles(), 100);
//     EXPECT_EQ(sim->nGenerations(), 10);
//     EXPECT_EQ(sim->nInactive(), 2);
//   }

//   EXPECT_EQ(getNTotalCells(test_sim_1), 10);
//   EXPECT_EQ(getNTotalCells(test_sim_2), 16);
//   EXPECT_EQ(getNTotalCells(test_sim_split), 26);

//   for (auto i = 0; i < test_sims.size(); i++)
//   {
//     EXPECT_EQ(regions(*test_sims[i]).size(), i + 1);
//   }

//   EXPECT_EQ(nFissionableRegions(test_sim_1), 1);
//   EXPECT_EQ(nFissionableRegions(test_sim_2), 1);
//   EXPECT_EQ(nFissionableRegions(test_sim_split), 1);
// }

// TEST_F(MCSlabTest, ShannonEntropy)
// {
//   ASSERT_EQ(getNTotalCells(test_sim_1), 10);
//   std::vector<unsigned long int> collision_bins(getNTotalCells(test_sim_1), 1);
//   EXPECT_TRUE(test_sim_1.shannonEntropy(collision_bins) > 0);
// }

// TEST_F(MCSlabTest, Absorption)
// {
//   Region test_region{-1, 1, 1, 1, 0, 4}; // hard-code because I'm tired of debugging
//   std::vector<Region> test_regions{test_region};

//   MCSlab test_abs{"../tests/input_files/test_pure_absorb.xml"};
//   Neutron test_neutron{0, test_regions};
//   EXPECT_TRUE(test_abs.testAbsorption(test_neutron));

//   test_abs.absorption(test_neutron);
//   EXPECT_FALSE(test_neutron.isAlive());

//   EXPECT_EQ(lenFissionBank(test_abs), 4);
// }

// TEST_F(MCSlabTest, Scatter)
// {
//   Region test_region{-1, 1, 1, 0, 1, 0}; // hard-code because I'm tired of debugging
//   std::vector<Region> test_regions{test_region};
//   MCSlab test_scatter{"../tests/input_files/test_pure_scatter.xml"};
//   Neutron test_neutron{0, test_regions};
//   double orig_mu = test_neutron.mu();

//   EXPECT_FALSE(test_scatter.testAbsorption(test_neutron));
//   test_scatter.scatter(test_neutron);
//   EXPECT_TRUE(test_neutron.isAlive());
//   EXPECT_FALSE(orig_mu == test_neutron.mu());
// }
