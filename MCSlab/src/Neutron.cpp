#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include <cmath>
#include <random>

Neutron::Neutron(double position) : _pos(position), _rng() {
  Neutron::randomIsoAngle();
  _is_alive = true;
}

double Neutron::distanceToCollision() {
  double mfp = 1.0 / _region->SigmaT();
  double rn = _rng.generateRN();
  return -std::log(rn) * mfp;
}

double Neutron::distanceToEdge() {
  if (_mu > 0)
    return _region->xMax() - _pos;
  else
    return _pos - _region->xMin();
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

void Neutron::movePositionAndRegion(const double new_position,
                                    const std::vector<Region> &regions) {
  _pos = new_position;
  setRegion(regions);
}

void Neutron::setRegion(const std::vector<Region> &regions) {
  _region = nullptr;
  for (auto region : regions) {
    if (region.xMin() < _pos && _pos < region.xMax()) {
      _region = &region;
      break;
    }
  }
}

void Neutron::kill() { _is_alive = false; }
