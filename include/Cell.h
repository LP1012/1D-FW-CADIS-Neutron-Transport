#pragma once

#include "CellBase.h"
#include "RNG.h"
#include <vector>

class Cell : public CellBase
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

  static double randomPositionInCell(const Cell cell, UniformRNG rng);

  // getters

  const double cellWeight() const { return _weight; }
  const double forwardFlux() const { return _forward_flux; }
  const double adjointFlux() const { return _adjoint_flux; }
  const double absorptionRatio() const { return _abs_ratio; }
  const double nPerAbsorption() const { return _n_per_abs; }

private:
  double _n_per_abs;
  double _abs_ratio;

  // set weight window values
  double _center_weight; // "target weight"
  double _upper_weight;
  double _lower_weight;

  double _weight;
  double _forward_flux;
  double _adjoint_flux;
};