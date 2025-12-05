#include "Neutron.h"
#include "Region.h"
#include <gtest/gtest.h>
#include <iostream> // delete later
#include <optional>
#include <vector>

#include <iomanip> // delete later

class RegionTest : public testing::Test
{
protected:
  Region region{-1.0, 1.0, 4, 3, 2, 4}; // set generic tester region
};

TEST_F(RegionTest, RegionIndexing)
{
  region.setIndex(1);
  EXPECT_TRUE(region.regionIndex() == 1); // basic test for setter
}

TEST_F(RegionTest, CellLocations)
{
  EXPECT_TRUE(region.cells().size() == 4);

  std::vector<std::vector<double>> expected_cell_locs{{-1, -0.5}, {-0.5, 0}, {0, 0.5}, {0.5, 1}};
  std::vector<double> expected_cell_centers = {-0.75, -0.25, 0.25, 0.75};
  for (auto i = 0; i < region.cells().size(); i++)
  {
    EXPECT_FLOAT_EQ(region.cells()[i].cellCenter(), expected_cell_centers[i]);
    EXPECT_FLOAT_EQ(region.cells()[i].xMin(), expected_cell_locs[i][0]);
    EXPECT_FLOAT_EQ(region.cells()[i].xMax(), expected_cell_locs[i][1]);
  }
}

TEST_F(RegionTest, MaterialProperties)
{
  EXPECT_FLOAT_EQ(region.SigmaT(), 5.0);
  EXPECT_FLOAT_EQ(region.absorptionRatio(), 0.6);
  EXPECT_FLOAT_EQ(region.nPerAbsorption(), 4.0 / 3.0);
}