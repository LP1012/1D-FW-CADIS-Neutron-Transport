#pragma once

#include <vector>
#include "Cell.h"

class Region
{
  friend class RegionTest;

public:
  Region(const double xmin,
         const double xmax,
         const unsigned int n_cells,
         const double Sigma_a,
         const double Sigma_s,
         const double nu_Sigma_f,
         const double source);

  static Region voidRegion(double xmin, double xmax, unsigned n_cells);

  void setIndex(const unsigned int index_number);

  // define getter funcions
  double xMin() const { return _xmin; }
  double xMax() const { return _xmax; }
  double nuSigF() const { return _nu_Sigma_f; }
  double SigmaA() const { return _Sigma_a; }
  double SigmaS() const { return _Sigma_s; }
  double SigmaT() const { return _Sigma_t; }
  double absorptionRatio() const { return _absorption_ratio; }
  double nPerAbsorption() const { return _n_per_abs; }
  unsigned int regionIndex() const { return _region_index; }
  unsigned int nCells() const { return _n_cells; }
  std::vector<Cell> cells() const { return _region_cells; }

private:
  // bounds of slab
  const double _xmin;
  const double _xmax;

  unsigned int _region_index; // location of region in region vector

  const unsigned int _n_cells; // number of cells in mesh

  // material properties
  const double _Sigma_a;    // absorption cross section
  const double _Sigma_s;    // scattering cross section
  const double _nu_Sigma_f; // fission cross section times number of neutrons produced per fission
  const double _source;
  double _Sigma_t;          // total cross section
  double _absorption_ratio; // absorptions per collision
  double _n_per_abs;        // neutrons produced per absorptions

  // store locations of cells in mesh
  std::vector<Cell> _region_cells;

  void populateCells();
};