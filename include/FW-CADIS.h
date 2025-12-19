#pragma once

#include <string>
#include "Region.h"
#include "Cell.h"
#include "SN.h"

class FWCADIS
{
public:
  FWCADIS(const std::string input_file_name);

  void runForwardFlux();
  void runAdjointFlux();
  void updateForwardFlux(const SN & simulation);
  void updateAdjointFlux(const SN & simulation);
  void setWeightWindows();
  void kEigenvalueMonteCarlo();

  // bool useIC() const { return _implicit_capture; }
  bool useFWCADIS() const { return _fw_cadis; }

private:
  const std::string _input_file_name;
  unsigned int _n_total_cells;
  bool _fw_cadis;
  // bool _implicit_capture;
  double _window_width;

  std::vector<Region> _regions;
  std::vector<Cell> _cells;

  // Monte Carlo settings:
  unsigned int _n_particles;
  unsigned int _n_generations;
  unsigned int _n_inactive;

  // SN settings
  unsigned int _quadrature_order;

  std::vector<double> _forward_flux;
  double _forward_k_eff;

  void readInput();
};