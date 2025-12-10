#include "Neutron.h"
#include "Point.h"
#include "Cell.h"
#include <cmath>
#include <memory>
#include <optional>
#include <random>

Neutron::Neutron(double position,
                 double mu,
                 Cell * cell,
                 std::optional<double> weight,
                 std::optional<unsigned int> seed)
  : _pos(position),
    _mu(mu),
    _cell(cell),
    _weight(weight.value_or(1.0)),
    _rng(seed.has_value() ? UniformRNG(0, 1.0, seed.value()) : UniformRNG(0, 1.0))
{
  _is_alive = true;
  checkNeutron();
}

double
Neutron::distanceToCollision()
{
  double mfp = 1.0 / _cell->sigmaT();
  double rn = _rng.generateRN();
  return -std::log(rn) * mfp;
}

double
Neutron::distanceToEdge()
{
  if (_mu > 0)
    return (_cell->xMax() - _pos) / _mu;
  else
    return (_pos - _cell->xMin()) / (-_mu);
}

double
Neutron::randomIsoAngle(UniformRNG rng)
{
  return 2 * rng.generateRN() - 1; // sample between -1 and 1
};

void
Neutron::movePositionAndCell(const double new_position, std::vector<Cell> & cells)
{
  _pos = new_position;
  setCell(cells);
  checkNeutron();
}

void
Neutron::movePositionWithinCell(const double new_position)
{
  _pos = new_position;
}

void
Neutron::setCell(std::vector<Cell> & cells)
{
  _cell = nullptr;
  for (auto & cell : cells)
  {
    if (cell.xMin() < _pos && _pos < cell.xMax())
    {
      _cell = &cell;
      return;
    }
    else if (cell.xMin() == _pos)
    {
      if (_mu > 0)
      {
        _cell = &cell;
        return;
      }
    }
    else if (_pos == cell.xMax())
    {
      if (_mu < 0)
      {
        _cell = &cell;
        return;
      }
    }
  }
  throw std::runtime_error("ERROR: Neutron::setRegion() could not assign region!");
}

void
Neutron::changeWeight(const double weight)
{
  _weight = weight;
}

void
Neutron::kill()
{
  _is_alive = false;
}

void
Neutron::checkNeutron()
{
  if (!(_pos <= _cell->xMax() && _pos >= _cell->xMin()))
    throw std::runtime_error("Neutron cell not set correctly! Position "
                             "not located within bounds!");
}