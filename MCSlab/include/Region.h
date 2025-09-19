#pragma once

#include <vector>

class Region {
public:
  Region(const unsigned int id, const double xmin, const double xmax,
         const unsigned int n_cells, const double Sigma_a, const double Sigma_s,
         const double nu_Sigma_f);

  // define getter funcions
  double xMin() { return _xmin; }
  double xMax() { return _xmax; }
  double mfp() { return _mfp; }
  unsigned int id() { return _id; }
  std::vector<std::vector<double>> cellLocs() { return _cell_locs; }
  std::vector<double> cellCenters() { return _cell_centers; }

private:
  // slab id
  const unsigned int _id;
  // bounds of slab
  const double _xmin;
  const double _xmax;

  const unsigned int _n_cells; // number of cells in mesh

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