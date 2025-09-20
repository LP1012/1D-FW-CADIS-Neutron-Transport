#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include <cmath>
#include <random>

Neutron::Neutron(double position, std::vector<Region> &regions)
    : _pos(position), _regions(regions), _rng() {
  Neutron::randomIsoAngle();
  Neutron::setRegionID();
}

double Neutron::distanceToCollision(const double mfp) {
  double rn = _rng.generateRN();
  return -std::log(rn) * mfp;
}

void Neutron::setRandomStartPosition(
    const std::vector<Region> &fissionable_regions) {

  unsigned int n_fissionable_regions = fissionable_regions.size();
  double region_width = 1.0 / static_cast<double>(n_fissionable_regions);

  double rn = _rng.generateRN();
  for (auto i = 0; i < n_fissionable_regions; i++) {
    if (static_cast<double>(i) *
        region_width<rn &&static_cast<double>(i + 1) * region_width> rn) {

      Region selected_region =
          fissionable_regions[i]; // region of particle location
      _pos = _rng.generateRN() *
                 (selected_region.xMax() - selected_region.xMin()) +
             selected_region.xMin(); // neutron position in region
      break;
    }
  }
}

void Neutron::randomIsoAngle() {
  _mu = 2 * _rng.generateRN() - 1; // sample between -1 and 1
  _ang = std::acos(_mu);
};

void Neutron::movePosition(const double new_position) { _pos = new_position; }

void Neutron::setRegionID() {
  _region_id = 0; // set to 0 (void) if not in region
  for (auto region : _regions) {
    if (region.xMin() < _pos && _pos < region.xMax()) {
      _region_id = region.id();
      break;
    }
  }
}
