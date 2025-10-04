#include "Region.h"

Region::Region(const double xmin, const double xmax, const unsigned int n_cells,
               const double Sigma_a, const double Sigma_s,
               const double nu_Sigma_f)
    : _xmin(xmin), _xmax(xmax), _n_cells(n_cells), _Sigma_a(Sigma_a),
      _Sigma_s(Sigma_s), _nu_Sigma_f(nu_Sigma_f) {

  // calculate total cross section and mean-free-path
  _Sigma_t = _Sigma_a + _Sigma_s;
  _absorption_ratio = _Sigma_a / _Sigma_t;
  _n_per_abs = _nu_Sigma_f / _Sigma_a;

  // populate cell locations
  Region::populateCellLocs();
};

Region Region::voidRegion(double xmin, double xmax, unsigned n_cells) {
  // Preset cross sections for void
  constexpr double Sigma_a_void = 0.0;
  constexpr double Sigma_s_void = 0.0;
  constexpr double nuSigma_f_void = 0.0;

  return Region(xmin, xmax, n_cells, Sigma_a_void, Sigma_s_void,
                nuSigma_f_void);
};

void Region::setIndex(const unsigned int region_index) {
  _region_index = region_index;
}

void Region::populateCellLocs() {
  double dx = (_xmax - _xmin) / static_cast<double>(_n_cells);
  for (auto i = 0; i < _n_cells; i++) {
    double lower_bound =
        static_cast<double>(i) * dx + _xmin; // determine lower bound
    double upper_bound =
        static_cast<double>(i + 1) * dx + _xmin; // determine upper bound
    std::vector<double> cell_loc = {lower_bound,
                                    upper_bound}; // create point set
    _cell_locs.push_back(cell_loc);               // add to cell locations
    _cell_centers.push_back((upper_bound + lower_bound) /
                            2.0); // add cell center location

    _cell_bounds.push_back(lower_bound);
    if (i == _n_cells - 1)
      _cell_bounds.push_back(upper_bound);
  }
};