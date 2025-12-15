#pragma once

#include <vector>
#include <array>

namespace discreteQuadrature
{
template <std::size_t N>
struct GaussLegendre;

template <>
struct GaussLegendre<2>
{
  static constexpr std::array<double, 2> weights = {1.0, 1.0};
  static constexpr std::array<double, 2> abscissa = {-0.5773502691896257, 0.5773502691896257};
};

template <>
struct GaussLegendre<4>
{
  static constexpr std::array<double, 4> weights = {
      0.6521451548625461, 0.6521451548625461, 0.3478548451374538, 0.3478548451374538};
};

template <>
struct GaussLegendre<8>
{
  static constexpr std::array<double, 8> weights = {0.3626837833783620,
                                                    0.3626837833783620,
                                                    0.3137066458778873,
                                                    0.3137066458778873,
                                                    0.2223810344533745,
                                                    0.2223810344533745,
                                                    0.1012285362903763,
                                                    0.1012285362903763};
  static constexpr std::array<double, 8> abscissa = {-0.1834346424956498,
                                                     0.1834346424956498,
                                                     -0.5255324099163290,
                                                     0.5255324099163290,
                                                     -0.7966664774136267,
                                                     0.7966664774136267,
                                                     -0.9602898564975363,
                                                     0.9602898564975363};
};

template <>
struct GaussLegendre<12>
{
  static constexpr std::array<double, 12> weights = {0.2491470458134028,
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
  static constexpr std::array<double, 12> abscissa = {-0.1252334085114689,
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
};

/// @brief Integrates grid function using Gauss-Legendre quadrature assuming the grid function has been evaluated using the roots of Legendre polynomials.
/// @tparam N
/// @param fs Grid function values.
/// @param a Lower integration bound.
/// @param b Upper integration bound.
/// @return
template <std::size_t N>
double
integrate(const std::vector<double> & fs, const double a, const double b)
{
  static_assert(N > 0, "Quadrature order must be positive");

  if (fs.size() != N)
    throw std::runtime_error("Grid function size does not match quadrature order");

  const double front_coeff = (b - a) / 2.0;

  double sum = 0.0;
  for (std::size_t i = 0; i < N; ++i)
    sum += GaussLegendre<N>::weights[i] * fs[i];

  return front_coeff * sum;
};

}