#pragma once

#include <array>

template <typename F>
double
gaussQuadrature(F f, double a, double b, unsigned short quad)
{
  switch (quad)
  {
    case 8:
    {
      std::array<double, 8> weights = {0.3626837833783620,
                                       0.3626837833783620,
                                       0.3137066458778873,
                                       0.3137066458778873,
                                       0.2223810344533745,
                                       0.2223810344533745,
                                       0.1012285362903763,
                                       0.1012285362903763};
      std::array<double, 8> abscissa = {-0.1834346424956498,
                                        0.1834346424956498,
                                        -0.5255324099163290,
                                        0.5255324099163290,
                                        -0.7966664774136267,
                                        0.7966664774136267,
                                        -0.9602898564975363,
                                        0.9602898564975363};

      double running_sum = 0;
      double front_coeff = (b - a) / 2.0;
      for (auto i = 0; i < quad; i++)
      {
        double t = (b + a) / 2.0 + front_coeff * abscissa[i];
        running_sum += weights[i] * f(t);
      }
      return front_coeff * running_sum;
    }

    case 12:
    {
      std::array<double, 12> weights = {0.2491470458134028,
                                        0.2491470458134028,
                                        0.2334925365383548,
                                        0.2334925365383548,
                                        0.2031674267230659,
                                        0.2031674267230659,
                                        0.1600783285433462,
                                        0.1600783285433462,
                                        0.1069393259953184,
                                        0.1069393259953184,
                                        0.0471753363865118,
                                        0.0471753363865118};
      std::array<double, 12> abscissa = {-0.1252334085114689,
                                         0.1252334085114689,
                                         -0.3678314989981802,
                                         0.3678314989981802,
                                         -0.5873179542866175,
                                         0.5873179542866175,
                                         -0.7699026741943047,
                                         0.7699026741943047,
                                         -0.9041172563704749,
                                         0.9041172563704749,
                                         -0.9815606342467192,
                                         0.9815606342467192};

      double running_sum = 0;
      double front_coeff = (b - a) / 2.0;
      for (auto i = 0; i < quad; i++)
      {
        double t = (b + a) / 2.0 + front_coeff * abscissa[i];
        running_sum += weights[i] * f(t);
      }
      return front_coeff * running_sum;
    }
    default:
      throw std::runtime_error("Invalid number of quadrature points selected!");
  }
}

template <typename F>
double
gauss2DQuadrature(F f, double x1, double x2, double y1, double y2, unsigned short quad)
{
  auto inner_integral = [&](double y)
  {
    auto fx = [&](double x) { return f(x, y); };
    return gaussQuadrature(fx, x1, x2, quad);
  };

  return gaussQuadrature(inner_integral, y1, y2, quad);
}
