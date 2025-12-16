#pragma once

#include <vector>
#include <stdexcept>

class CellBase
{
public:
  CellBase(const double xmin,
           const double xmax,
           const double Sigma_a,
           const double Sigma_s,
           const double nu_Sigma_f,
           const double volumetric_source);

  template <typename T>
  static unsigned int
  cellIndex(const double position, const double mu, const std::vector<T> & cells)
  {
    for (auto i = 0; i < cells.size(); i++)
    {
      const double cell_max = cells[i].xMax();
      const double cell_min = cells[i].xMin();
      if (position < cell_max && position > cell_min)
        return i;
      else if (isOnBoundary(position, cell_max) && mu < 0)
        return i;
      else if (isOnBoundary(position, cell_min) && mu > 0)
        return i;
    }
    throw std::runtime_error("Position not located within given cells!");
  };

  static bool isOnBoundary(const double pos, const double boundary);

  const double xMin() const { return _xmin; }
  const double xMax() const { return _xmax; }
  const double cellCenter() const { return _cell_center; }
  const double cellWidth() const { return _cell_width; }
  const double sigmaT() const { return _Sigma_t; }
  const double sigmaA() const { return _Sigma_a; }
  const double sigmaS() const { return _Sigma_s; }
  const double nuSigmaF() const { return _nu_Sigma_f; }
  const double volumetricSource() const { return _vol_source; }

protected:
  const double _xmin;
  const double _xmax;
  const double _Sigma_a;
  const double _Sigma_s;
  const double _nu_Sigma_f;
  const double _vol_source;

  double _Sigma_t;
  double _cell_width;
  double _cell_center;
};