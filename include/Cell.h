#pragma once

class Cell
{
public:
  Cell(const double xmin,
       const double xmax,
       const double Sigma_a,
       const double Sigma_s,
       const double nu_Sigma_f);

  void setWeight(const double weight);
  void setForwardFlux(const double forward_flux);
  void setAdjointFlux(const double adjoint_flux);

  // getters
  const double xMin() const { return _xmin; }
  const double xMax() const { return _xmax; }
  const double cellWeight() const { return _weight; }
  const double forwardFlux() const { return _forward_flux; }
  const double adjointFlux() const { return _adjoint_flux; }

private:
  const double _xmin;
  const double _xmax;
  const double _Sigma_a;
  const double _Sigma_s;
  const double _nu_Sigma_f;

  double _Sigma_t;
  double _cell_width;
  double _cell_center;
  double _weight;
  double _forward_flux;
  double _adjoint_flux;
};