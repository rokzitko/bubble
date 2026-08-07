#!/usr/bin/env python3

"""Generate high-precision references for the built-in clean-limit kernels.

This script is not needed by the normal test suite. It requires mpmath and
prints data suitable for clean_kernel_references.dat.
"""

import mpmath as mp


mp.mp.dps = 100


def phi(kernel, epsilon):
    if kernel == 1:
        return mp.mpf(1)
    if kernel == 2:
        return epsilon
    if kernel == 3:
        return mp.sqrt(max(mp.mpf(0), 1 - epsilon * epsilon))
    if kernel == 4:
        return epsilon * mp.sqrt(max(mp.mpf(0), 1 - epsilon * epsilon))
    if kernel == 5:
        return max(mp.mpf(0), 1 - epsilon * epsilon) ** mp.mpf("1.5")
    if kernel == 6:
        return epsilon * max(mp.mpf(0), 1 - epsilon * epsilon) ** mp.mpf("1.5")
    gaussian = mp.exp(-2 * epsilon * epsilon)
    return gaussian if kernel == 7 else epsilon * gaussian


def transformed_integrand(theta, kernel, power, x, y):
    cosine = mp.cos(theta)
    epsilon = x + y * mp.tan(theta)
    return phi(kernel, epsilon) * cosine ** (2 * power - 2) / (
        mp.pi**power * y ** (power - 1)
    )


def reference(kernel, power, x_text, y_text):
    # Match the exact binary64 inputs parsed by the C++ test executable. This
    # matters when a point is only a few linewidths from a band edge.
    x = mp.mpf(float(x_text))
    y = mp.mpf(float(y_text))

    if kernel <= 6 and x > 1 + 32 * y:
        # There is no Lorentzian peak in the integration interval. Direct
        # energy quadrature avoids subtracting angles close to -pi/2.
        def direct(epsilon):
            spectral = y / (mp.pi * ((x - epsilon) ** 2 + y * y))
            return phi(kernel, epsilon) * spectral**power

        return mp.quad(direct, [-1, -mp.mpf("0.5"), 0, mp.mpf("0.5"), 1])

    if kernel <= 6:
        energies = [-1, -mp.mpf("0.5"), 0, mp.mpf("0.5"), 1]
        points = [mp.atan((energy - x) / y) for energy in energies]
    else:
        energies = [-12, -8, -5, -3, -2, -1, 0, 1, 2, 3, 5, 8, 12]
        points = [-mp.pi / 2]
        points.extend(mp.atan((mp.mpf(energy) - x) / y) for energy in energies)
        points.extend([mp.mpf(0), mp.pi / 2])
        points = sorted(set(points))

    return mp.quad(
        lambda theta: transformed_integrand(theta, kernel, power, x, y),
        points,
    )


def cases():
    result = []

    for kernel in range(1, 9):
        for power in (2, 3):
            result.extend(
                (kernel, power, x, y)
                for x, y in (("0.37", "0.2"), ("0.37", "1e-8"))
            )

    for kernel in range(1, 7):
        for power in (2, 3):
            result.extend(
                (kernel, power, x, y)
                for x, y in (
                    ("1", "1e-8"),
                    ("1.000001", "7.49e-7"),
                    ("2", "1e-8"),
                    ("20", "0.1"),
                )
            )

    for kernel in (1, 2):
        for power in (2, 3):
            result.extend(
                (kernel, power, "1.001", y)
                for y in ("0.000749", "0.000751")
            )

    for kernel in (7, 8):
        for power in (2, 3):
            result.extend(
                (kernel, power, x, y)
                for x, y in (
                    ("3.75", "1e-8"),
                    ("5", "1e-12"),
                    ("8", "1e-12"),
                    ("13.3", "0.01"),
                    ("0.5", "5"),
                    ("1", "8"),
                    ("5.5", "0.3"),
                    ("2", "0.001"),
                    ("2", "1e-18"),
                    ("4.19", "1e-16"),
                    ("1e-12", "1e-16"),
                )
            )

    return result


print("# m n x y expected relative_tolerance absolute_tolerance")
for kernel, power, x, y in cases():
    value = reference(kernel, power, x, y)
    print(
        kernel,
        power,
        x,
        y,
        mp.nstr(value, 25),
        "5e-11",
        "1e-300",
    )
