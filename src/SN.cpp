
#include "SN.h"
#include "Cell.h"
#include "SNCell.h"

#include <stdexcept>
#include <vector>
#include <cmath>

SN::SN(const std::vector<Cell> & cells, const unsigned int GQ_order) : _gq_order(GQ_order)
{
  if (GQ_order % 2 == 1.0)
    throw std::runtime_error("Gauss quadrature order cannot be odd!");
  populateSNCells(cells);
}

void
SN::populateSNCells(const std::vector<Cell> & cells)
{
  for (auto & cell : cells)
    _sn_cells.push_back(SNCell(cell));
}