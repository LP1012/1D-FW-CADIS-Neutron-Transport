#pragma once

#include <vector>

class Cell
{
public:
  Cell(const double xmin,
       const double xmax,
       const double Sigma_a,
       const double Sigma_s,
       const double nu_Sigma_f);

  void setWeight(const double center_weight, const double upper_weight, const double lower_weight);
  void setForwardFlux(const double forward_flux) { _forward_flux = forward_flux; };
  void setAdjointFlux(const double adjoint_flux) { _adjoint_flux = adjoint_flux; };

  static unsigned int
  cellIndex(const double position, const double mu, const std::vector<Cell> cells);

  bool static isOnBoundary(const double pos, const double boundary);

  // getters
  const double xMin() const { return _xmin; }
  const double xMax() const { return _xmax; }
  const double cellWeight() const { return _weight; }
  const double forwardFlux() const { return _forward_flux; }
  const double adjointFlux() const { return _adjoint_flux; }
  const double cellCenter() const { return _cell_center; }
  const double sigmaT() const { return _Sigma_t; }
  const double sigmaA() const { return _Sigma_a; }
  const double sigmaS() const { return _Sigma_s; }
  const double nuSigmaF() const { return _nu_Sigma_f; }
  const double absorptionRatio() const { return _abs_ratio; }
  const double nPerAbsorption() const { return _n_per_abs; }

private:
  const double _xmin;
  const double _xmax;
  const double _Sigma_a;
  const double _Sigma_s;
  const double _nu_Sigma_f;

  double _n_per_abs;
  double _abs_ratio;

  double _Sigma_t;
  double _cell_width;
  double _cell_center;

  // set weight window values
  double _center_weight; // "target weight"
  double _upper_weight;
  double _lower_weight;

  double _weight;
  double _forward_flux;
  double _adjoint_flux;
};