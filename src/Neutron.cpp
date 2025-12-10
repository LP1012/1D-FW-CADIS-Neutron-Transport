#include "Neutron.h"
#include "Point.h"
#include "Cell.h"
#include <cmath>
#include <memory>
#include <optional>
#include <random>

Neutron::Neutron(double position,
                 double mu,
                 Cell & cell,
                 std::optional<double> weight,
                 std::optional<unsigned int> seed)
  : _weight(weight.value_or(1.0)),
    _rng(seed.has_value() ? UniformRNG(0, 1.0, seed.value()) : UniformRNG(0, 1.0))
{
  _is_alive = true;
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

void
Neutron::setRandomStartPosition(const std::vector<Cell> & fissionable_cells,
                                std::vector<Cell> & cells)
{
  const double eps = 1e-12;
  unsigned int n_fissionable_regions = fissionable_cells.size();
  double region_width = 1.0 / static_cast<double>(n_fissionable_regions);
  double rn = _rng.generateRN();
  for (auto i = 0; i < n_fissionable_regions; i++)
  {
    if ((static_cast<double>(i) * region_width < rn) &&
        (static_cast<double>(i + 1) * region_width > rn))
    {
      Cell selected_cell = fissionable_cells[i]; // region of particle location

      while (true)
      {
        _pos = _rng.generateRN() * (selected_cell.xMax() - selected_cell.xMin()) +
               selected_cell.xMin();

        bool bad_left = (std::abs(_pos - selected_cell.xMin()) < eps && _mu < 0);
        bool bad_right = (std::abs(_pos - selected_cell.xMax()) < eps && _mu > 0);

        if (!bad_left && !bad_right)
          break;
      }

      setCell(cells);
      return;
    }
  }
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
}

void
Neutron::movePositionWithinRegion(const double new_position)
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