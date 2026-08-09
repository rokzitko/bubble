#!/usr/bin/env python3

"""Generate high-precision full-domain references for built-in optical kernels."""

import mpmath as mp
import math


mp.mp.dps = 250


def binary64(text):
    return mp.mpf(float(text))


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


def transformed_segment(kernel, anchor, other, lower, upper):
    anchor_x, anchor_y = anchor
    other_x, other_y = other
    theta_lower = -mp.pi / 2 if lower == -mp.inf else mp.atan((lower - anchor_x) / anchor_y)
    theta_upper = mp.pi / 2 if upper == mp.inf else mp.atan((upper - anchor_x) / anchor_y)

    def integrand(theta):
        epsilon = anchor_x + anchor_y * mp.tan(theta)
        other_spectral = other_y / (
            mp.pi * ((other_x - epsilon) ** 2 + other_y**2)
        )
        return phi(kernel, epsilon) * other_spectral / mp.pi

    points = [theta_lower]
    other_theta = mp.atan((other_x - anchor_x) / anchor_y)
    if theta_lower < other_theta < theta_upper:
        points.append(other_theta)
    points.append(theta_upper)
    return mp.quad(integrand, points)


def reference(kernel, x1_text, y1_text, x2_text, y2_text):
    first = (binary64(x1_text), binary64(y1_text))
    second = (binary64(x2_text), binary64(y2_text))
    lower = mp.mpf(-1) if kernel <= 6 else -mp.inf
    upper = mp.mpf(1) if kernel <= 6 else mp.inf
    separation = abs(first[0] - second[0])

    if separation > 8 * max(first[1], second[1]):
        left, right = (first, second) if first[0] <= second[0] else (second, first)
        split = (left[0] + right[0]) / 2
        if lower < split < upper:
            return transformed_segment(kernel, left, right, lower, split) + transformed_segment(
                kernel, right, left, split, upper
            )

    anchor, other = (first, second) if first[1] <= second[1] else (second, first)
    return transformed_segment(kernel, anchor, other, lower, upper)


def cases():
    for kernel in range(1, 9):
        yield kernel, "0.2", "0.05", "0.35", "0.08"

    for kernel in range(1, 7):
        yield kernel, "1e8", "1e-3", "100000000.00001", "1e-3"

    for kernel in (7, 8):
        yield kernel, "8.5", "1e-100", "8.50000000000091", "1e-100"

    for kernel in range(1, 9):
        yield kernel, "0.37", "0.2", "0.37", "0.2"
        yield kernel, "0.37", "0.2", "0.37000000000000005", "0.2"
        yield kernel, "0.2", "0.1", "0.204", "0.1"
        yield kernel, "0.2", "0.1", "0.206", "0.1"

    for kernel in (7, 8):
        yield kernel, "7.999999999999999", "1e-20", "8.000000000000002", "1e-20"
        yield kernel, "6", "6", "6.1", "6"
        pole_linewidth = 0.25/8.5
        for linewidth in (
            math.nextafter(pole_linewidth, 0.0),
            pole_linewidth,
            math.nextafter(pole_linewidth, math.inf),
        ):
            yield kernel, "8.5", repr(linewidth), "8.501", repr(linewidth)


print("# m x1 y1 x2 y2 expected relative_tolerance absolute_tolerance")
for values in cases():
    value = reference(*values)
    print(*values, mp.nstr(value, 25), "5e-8", "1e-300")
