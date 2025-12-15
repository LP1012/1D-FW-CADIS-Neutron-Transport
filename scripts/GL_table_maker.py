import numpy as np
from numpy.polynomial.legendre import leggauss  # computes abscissa and weights

quadratures = np.arange(2, 101)

with open("gl_weights_and_abscissa.txt", "w") as f:
    for quad in quadratures:
        ordinates, weights = leggauss(quad)
        print("template<>", file=f)
        print(f"struct GaussLegendre<{quad}>", file=f)
        print("{", file=f)
        print(f"static constexpr std::array<double, {quad}> weights = ", end="", file=f)
        print("{", end="", file=f)
        for val in weights[:-1]:
            print(f"{val}, ", end="", file=f)
        print(f"{weights[-1]}", end="", file=f)
        print("};", file=f)
        print(
            f"static constexpr std::array<double, {quad}> abscissa = ",
            end="",
            file=f,
        )
        print("{", end="", file=f)
        for val in ordinates[:-1]:
            print(f"{val}, ", end="", file=f)
        print(f"{ordinates[-1]}", end="", file=f)
        print("};", file=f)
        print("};\n", file=f)
