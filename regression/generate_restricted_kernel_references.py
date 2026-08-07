#!/usr/bin/env python3

"""Generate high-precision references for explicitly bounded epsilon kernels."""

import mpmath as mp


mp.mp.dps = 180


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


def phi_at_offset(kernel, center, offset):
    epsilon = center + offset
    if kernel < 3 or kernel > 6:
        return phi(kernel, epsilon)
    one_minus_square = max(mp.mpf(0), (1 - center - offset) * (1 + center + offset))
    root = mp.sqrt(one_minus_square)
    if kernel == 3:
        return root
    if kernel == 4:
        return epsilon * root
    if kernel == 5:
        return one_minus_square * root
    return epsilon * one_minus_square * root


def distance_reference(kernel, power, x, y, lower, upper):
    scale = y ** (1 - power) / mp.pi**power

    def integrate_distances(direction, start, stop):
        if not start < stop:
            return mp.mpf(0)
        total = mp.mpf(0)
        if start < 1:
            direct_stop = min(stop, mp.mpf(1))
            total += mp.quad(
                lambda distance: phi_at_offset(kernel, x, direction * y * distance)
                / (1 + distance**2) ** power,
                [start, direct_stop],
                maxdegree=10,
            )
            start = direct_stop
        if start < stop:
            log_start = mp.log(start)
            log_stop = mp.log(stop)
            step = mp.log(8)
            points = [log_start]
            point = log_start + step
            while point < log_stop:
                points.append(point)
                point += step
            points.append(log_stop)

            def log_integrand(argument):
                distance = mp.exp(argument)
                return (
                    phi_at_offset(kernel, x, direction * y * distance)
                    * distance
                    / (1 + distance**2) ** power
                )

            total += mp.fsum(
                mp.quad(log_integrand, [a, b], maxdegree=10)
                for a, b in zip(points, points[1:])
            )
        return total

    if x < lower:
        return scale * integrate_distances(1, (lower - x) / y, (upper - x) / y)
    if x > upper:
        return scale * integrate_distances(-1, (x - upper) / y, (x - lower) / y)
    return scale * (
        integrate_distances(-1, 0, (x - lower) / y)
        + integrate_distances(1, 0, (upper - x) / y)
    )


def reference(kernel, power, x_text, y_text, lower_text, upper_text):
    x = binary64(x_text)
    y = binary64(y_text)
    lower = binary64(lower_text)
    upper = binary64(upper_text)

    if kernel % 2 == 0 and lower == -upper and (power == 0 or x == 0):
        return mp.mpf(0)

    if power == 0:
        points = [lower]
        if kernel >= 7:
            points.extend(point for point in map(mp.mpf, (-10, -5, -2, 0, 2, 5, 10))
                          if lower < point < upper)
        else:
            points.append((lower + upper) / 2)
        points.append(upper)
        return mp.quad(lambda epsilon: phi(kernel, epsilon), points)

    if abs(y) < mp.mpf("1e-80") and power == 2:
        if kernel == 5 and ((x == upper == 1) or (x == lower == -1)):
            return mp.sqrt(abs(y)) / (2 * mp.pi)
        if lower < x < upper and phi(kernel, x) != 0:
            return phi(kernel, x) / (2 * mp.pi * abs(y))

    if abs(y) < mp.mpf("1e-80"):
        return distance_reference(kernel, power, x, abs(y), lower, upper)

    theta_lower = mp.atan((lower - x) / y)
    theta_upper = mp.atan((upper - x) / y)
    theta_middle = min(theta_upper, max(theta_lower, mp.mpf(0)))

    def integrand(theta):
        cosine = mp.cos(theta)
        epsilon = x + y * mp.tan(theta)
        return phi(kernel, epsilon) * cosine ** (2 * power - 2) / (
            mp.pi**power * y ** (power - 1)
        )

    points = [theta_lower]
    if theta_lower < theta_middle < theta_upper:
        points.append(theta_middle)
    points.append(theta_upper)
    return mp.quad(integrand, points)


def cases():
    result = []
    for kernel in range(1, 9):
        for power in range(4):
            result.append((kernel, power, "0.17", "0.08", "-0.45", "0.6"))

    for kernel in range(1, 9):
        for power in (1, 2, 3):
            result.append((kernel, power, "0.21", "1e-12", "-0.4", "0.55"))

    for kernel in range(1, 7):
        for power in (1, 2, 3):
            result.append((kernel, power, "0.99", "1e-8", "0.8", "1"))
            result.append((kernel, power, "2", "1e-8", "-0.5", "0.5"))

    for kernel in (7, 8):
        for power in (1, 2, 3):
            result.append((kernel, power, "4.19", "1e-16", "-0.75", "0.9"))
            result.append((kernel, power, "4.19", "1e-16", "4", "4.4"))

    for kernel in (7, 8):
        result.append((kernel, 0, "0", "0.1", "-100000", "100000"))

    for kernel in (2, 4, 6, 8):
        for power in (1, 2, 3):
            result.append((kernel, power, "0", "1e-8", "-0.5", "0.5"))

    result.extend((
        (1, 1, "-0.5", "1", "0.5", "0.5000000000000001"),
        (1, 2, "0.500000000017", "1e-12", "-0.5", "0.5"),
        (3, 1, "1", "1e-20", "0", "1"),
        (7, 0, "0", "0.1", "-1e308", "1e308"),
        (5, 2, "1", "1e-300", "0", "1"),
        (7, 2, "20", "1e-300", "19", "21"),
        (2, 3, "1e-120", "1e-100", "-0.5", "0.5"),
        (1, 2, "-4.097e-197", "1e-200", "0", "1"),
        (1, 2, "0.500000000001", "1e-16", "-0.5", "0.5"),
        (2, 1, "1.7e21", "1e20", "-0.5", "0.5"),
    ))

    return result


print("# m n x y lower upper expected relative_tolerance absolute_tolerance")
for values in cases():
    value = reference(*values)
    print(*values, mp.nstr(value, 25), "5e-9", "1e-290")
