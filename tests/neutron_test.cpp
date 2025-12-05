#include "Neutron.h"
#include "Region.h"
#include <gtest/gtest.h>
#include <iostream> // delete later
#include <optional>
#include <vector>

#include <iomanip> // delete later

class NeutronTest : public testing::Test {
protected:
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
  neutron.setRandomStartPosition(fissionable_regions, regions);
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

  // neutron inside new region
  neutron.movePositionAndRegion(1.5, regions);
  EXPECT_TRUE(neutron.region().regionIndex() == regions[1].regionIndex());

  // neutron on boundary (forward)
  setMu(1);
  neutron.movePositionAndRegion(1, regions);
  EXPECT_TRUE(neutron.region().regionIndex() == regions[1].regionIndex());

  // neutron on boundary (backward)
  setMu(-1);
  neutron.movePositionAndRegion(1, regions);
  EXPECT_TRUE(neutron.region().regionIndex() == regions[0].regionIndex());
}

TEST_F(NeutronTest, Kill) {
  neutron.kill();
  EXPECT_FALSE(neutron.isAlive());
}