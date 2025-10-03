#include "Neutron.h"
#include "Region.h"
#include <gtest/gtest.h>
#include <optional>
#include <vector>

class NeutronTest : public testing::Test {
protected:
  void SetUp() override {
    regions.clear();
    regions.emplace_back(-1, 1, 1, 1, 1, 0);
    regions.emplace_back(1, 2, 1, 1, 1, 1);

    fissionable_regions.clear();
    fissionable_regions.emplace_back(regions[1]);

    neutron.emplace(0, regions);
  }
  ~NeutronTest() noexcept override =
      default; // fixes the exception specification

  std::vector<Region> regions;
  std::vector<Region> fissionable_regions;
  std::optional<Neutron> neutron; // can be constructed later
};

TEST_F(NeutronTest, InitialState) {
  EXPECT_TRUE(neutron->isAlive());
  EXPECT_NEAR(neutron->pos(), 0.0, 1e-12);
}

TEST_F(NeutronTest, SetRandomPositionInFuel) {
  neutron->setRandomStartPosition(fissionable_regions);
  EXPECT_TRUE(neutron->pos() < fissionable_regions[0].xMax());
  EXPECT_TRUE(neutron->pos() > fissionable_regions[0].xMin());
}
