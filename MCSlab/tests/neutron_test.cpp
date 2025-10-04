#include "Neutron.h"
#include "Region.h"
#include <gtest/gtest.h>
#include <iostream> // delete later
#include <optional>
#include <vector>

#include <iomanip> // delete later

class NeutronTest : public testing::Test {
protected:
  // void SetUp() override {
  //   regions.clear();
  //   regions.emplace_back(-1, 1, 1, 1, 1, 0);
  //   regions.emplace_back(1, 2, 1, 1, 1, 1);

  //   fissionable_regions.clear();
  //   fissionable_regions.emplace_back(regions[1]);

  //   neutron.emplace(0, regions, 0); // set the seed to be 0 for testing
  // }
  // ~NeutronTest() noexcept override =
  //     default; // fixes the exception specification

  // std::vector<Region> regions;
  // std::vector<Region> fissionable_regions;
  // std::optional<Neutron> neutron; // can be constructed later
  std::vector<Region> regions{Region(-1, 1, 1, 1, 1, 0),
                              Region(1, 2, 1, 1, 1, 1)};
  Neutron neutron{0, regions, 0};

  void setMu(double mu) { neutron._mu = mu; }
};

TEST_F(NeutronTest, InitialState) {
  EXPECT_TRUE(neutron.isAlive());
  EXPECT_NEAR(neutron.pos(), 0.0, 1e-12);
  EXPECT_NEAR(neutron.mu(), 0.185689233033, 1e-12);
}

TEST_F(NeutronTest, SetRandomPositionInFuel) {
  std::vector<Region> fissionable_regions = {regions[1]};
  neutron.setRandomStartPosition(fissionable_regions);
  EXPECT_TRUE(neutron.pos() < fissionable_regions[0].xMax());
  EXPECT_TRUE(neutron.pos() > fissionable_regions[0].xMin());
}

TEST_F(NeutronTest, DistanceToEdge) {
  setMu(0.5); // corresponds to theta=pi/3
  EXPECT_NEAR(neutron.distanceToEdge(), 2, 1e-12);
}

TEST_F(NeutronTest, DistanceToCollision) {
  EXPECT_NEAR(neutron.distanceToCollision(), 0.084643985540910294, 1e-12);
}

TEST_F(NeutronTest, MovePositionAndRegion) {
  neutron.movePositionAndRegion(1.5, regions);
  EXPECT_TRUE(neutron.region().regionIndex() == regions[1].regionIndex());
}