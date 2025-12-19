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
       const double nu_Sigma_f,
       const double source);

  void createWeightWindow(const double fwcadis_adjoint_flux, const double window_width);
  void createSimpleWWs();
  void setForwardFlux(const double forward_flux) { _forward_flux = forward_flux; };
  void setAdjointFlux(const double adjoint_flux) { _adjoint_flux = adjoint_flux; };
  void addToPathlength(const double val) { _pathlength_bin += val; }
  void addToCollisions(const double val) { _collision_bin += val; }
  void addPathlengthSS(const double val) { _pathlength_ss += val; }

  static double randomPositionInCell(const Cell * cell, UniformRNG rng);

  // getters

  const double lowerWeight() const { return _lower_weight; }
  const double upperWeight() const { return _upper_weight; }
  const double targetWeight() const { return _center_weight; }
  const double forwardFlux() const { return _forward_flux; }
  const double adjointFlux() const { return _adjoint_flux; }
  const double absorptionRatio() const { return _abs_ratio; }
  const double nPerAbsorption() const { return _n_per_abs; }
  const double absorptionProbability() const { return _abs_prob; }
  const double pathlength() const { return _pathlength_bin; }
  const double collision() const { return _collision_bin; }
  const double pathlengthSS() const { return _pathlength_ss; }

private:
  double _n_per_abs;
  double _abs_ratio;
  double _abs_prob;

  // set weight window values
  double _center_weight; // "target weight"
  double _upper_weight;
  double _lower_weight;

  double _forward_flux;
  double _adjoint_flux;

  double _pathlength_bin;
  double _collision_bin;
  double _pathlength_ss; // sum of squares
};