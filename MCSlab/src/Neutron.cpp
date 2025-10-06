#include "Neutron.h"
#include "Point.h"
#include "Region.h"
#include <cmath>
#include <memory>
#include <optional>
#include <random>

Neutron::Neutron(double position, std::vector<Region> &regions,
                 std::optional<unsigned int> seed)
    : _rng(seed.has_value() ? UniformRNG(0, 1.0, seed.value())
                            : UniformRNG(0, 1.0)) {
  Neutron::randomIsoAngle();
  Neutron::movePositionAndRegion(position, regions);
  _is_alive = true;
}

double Neutron::distanceToCollision() {
  double mfp = 1.0 / _region->SigmaT();
  double rn = _rng.generateRN();
  return -std::log(rn) * mfp;
}

double Neutron::distanceToEdge() {
  if (_mu > 0)
    return (_region->xMax() - _pos) / _mu;
  else
    return (_pos - _region->xMin()) / (-_mu);
}

void Neutron::setRandomStartPosition(
    const std::vector<Region> &fissionable_regions,
    std::vector<Region> &regions) {

  unsigned int n_fissionable_regions = fissionable_regions.size();
  double region_width = 1.0 / static_cast<double>(n_fissionable_regions);

  double rn = _rng.generateRN();
  for (auto i = 0; i < n_fissionable_regions; i++) {
    if ((static_cast<double>(i) * region_width < rn) &&
        (static_cast<double>(i + 1) * region_width > rn)) {

      Region selected_region =
          fissionable_regions[i]; // region of particle location
      _pos = _rng.generateRN() *
                 (selected_region.xMax() - selected_region.xMin()) +
             selected_region.xMin(); // neutron position in region
      setRegion(regions);
      break;
    }
  }
}

void Neutron::randomIsoAngle() {
  _mu = 2 * _rng.generateRN() - 1; // sample between -1 and 1
};

void Neutron::movePositionAndRegion(const double new_position,
                                    std::vector<Region> &regions) {
  _pos = new_position;
  setRegion(regions);
}

void Neutron::movePositionWithinRegion(const double new_position) {
  _pos = new_position;
}

void Neutron::setRegion(std::vector<Region> &regions) {
  _region = nullptr;
  for (auto &region : regions) {
    if (region.xMin() < _pos && _pos < region.xMax()) {
      _region = &region;
      break;
    } else if (region.xMin() == _pos) {
      if (_mu > 0)
        _region = &region;
    } else if (_pos == region.xMax()) {
      if (_mu < 0)
        _region = &region;
    }
  }
}

void Neutron::kill() { _is_alive = false; }