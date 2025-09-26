#pragma once

#include <vector>

class Region {
public:
  Region(const double xmin, const double xmax, const unsigned int n_cells,
         const double Sigma_a, const double Sigma_s, const double nu_Sigma_f);

  static Region voidRegion(double xmin, double xmax, unsigned n_cells);

  // define getter funcions
  double xMin() const { return _xmin; }
  double xMax() const { return _xmax; }
  double nuSigF() const { return _nu_Sigma_f; }
  double SigmaT() const { return _Sigma_t; }
  std::vector<std::vector<double>> cellLocs() const { return _cell_locs; }
  std::vector<double> cellCenters() const { return _cell_centers; }

private:
  // bounds of slab
  const double _xmin;
  const double _xmax;

  const unsigned int _n_cells; // number of cells in mesh

  // material properties
  const double _Sigma_a;
  const double _Sigma_s;
  const double _nu_Sigma_f;
  double _Sigma_t; // total cross section

  // store locations of cells in mesh
  std::vector<std::vector<double>> _cell_locs;
  std::vector<double> _cell_centers;
  void populateCellLocs();
};