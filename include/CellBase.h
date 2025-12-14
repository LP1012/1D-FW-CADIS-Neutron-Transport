#pragma once

#include <vector>

class CellBase
{
public:
  CellBase(const double xmin,
           const double xmax,
           const double Sigma_a,
           const double Sigma_s,
           const double nu_Sigma_f);

  template <typename T>
  static unsigned int
  cellIndex(const double position, const double mu, const std::vector<T> & cells);

  static bool isOnBoundary(const double pos, const double boundary);

  const double xMin() const { return _xmin; }
  const double xMax() const { return _xmax; }
  const double cellCenter() const { return _cell_center; }
  const double sigmaT() const { return _Sigma_t; }
  const double sigmaA() const { return _Sigma_a; }
  const double sigmaS() const { return _Sigma_s; }
  const double nuSigmaF() const { return _nu_Sigma_f; }

protected:
  const double _xmin;
  const double _xmax;
  const double _Sigma_a;
  const double _Sigma_s;
  const double _nu_Sigma_f;

  double _Sigma_t;
  double _cell_width;
  double _cell_center;
};