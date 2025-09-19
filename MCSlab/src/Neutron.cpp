#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include <cmath>
#include <random>

Neutron::Neutron(double position, std::vector<Region> &regions)
    : _pos(position), _regions(regions), _rng() {
  Neutron::randomIsoAngle();
  Neutron::setRegionID(regions);
};

void Neutron::randomIsoAngle() {
  _mu = 2 * _rng.generateRN() - 1; // sample between -1 and 1
  _ang = std::acos(_mu);
};

void Neutron::movePosition(const double new_position) { _pos = new_position; }

void Neutron::setRegionID(const std::vector<Region> &regions) {
  _region_id = 0; // set to 0 (void) if not in region
  for (auto region : regions) {
    if (region.xMin() < _pos && _pos < region.xMax()) {
      _region_id = region.id();
      break;
    }
  }
}
