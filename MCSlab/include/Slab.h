#pragma once

#include <vector>

class Slab {
public:
  Slab(const double xmin, const double xmax, const unsigned int n_cells,
       const double Sigma_a, const double Sigma_s, const double nu_Sigma_f);

  // define getter funcions
  double xMin() { return _xmin; }
  double xMax() { return _xmax; }
  double mfp() { return _mfp; }
  std::vector<std::vector<double>> cellLocs() { return _cell_locs; }
  std::vector<double> cellCenters() { return _cell_centers; }

private:
  // bounds of slab
  const double _xmin;
  const double _xmax;
  // number of cells in mesh
  const unsigned int _n_cells;
  // material properties
  const double _Sigma_a;
  const double _Sigma_s;
  const double _nu_Sigma_f;
  double _Sigma_t; // total cross section
  double _mfp;     // mean free path

  // store locations of cells in mesh
  std::vector<std::vector<double>> _cell_locs;
  std::vector<double> _cell_centers;
  void populateCellLocs();
};