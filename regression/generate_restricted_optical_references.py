#!/usr/bin/env python3

"""Generate high-precision references for bounded two-peak optical kernels."""

import mpmath as mp


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


def segment(kernel, anchor, other, lower, upper):
    anchor_x, anchor_y = anchor
    other_x, other_y = other

    if anchor_x < lower or anchor_x > upper:
        direction = 1 if anchor_x < lower else -1
        near = (lower - anchor_x) / anchor_y if direction > 0 else (anchor_x - upper) / anchor_y
        far = (upper - anchor_x) / anchor_y if direction > 0 else (anchor_x - lower) / anchor_y
        log_span = mp.log(far / near)

        def log_integrand(argument):
            distance = near * mp.exp(argument)
            epsilon = anchor_x + direction * anchor_y * distance
            other_spectral = other_y / (
                mp.pi * ((other_x - epsilon) ** 2 + other_y**2)
            )
            return (
                phi(kernel, epsilon)
                * other_spectral
                * distance
                / (mp.pi * (1 + distance**2))
            )

        return mp.quad(log_integrand, [0, log_span])

    theta_lower = mp.atan((lower - anchor_x) / anchor_y)
    theta_upper = mp.atan((upper - anchor_x) / anchor_y)

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


def reference(kernel, x1_text, y1_text, x2_text, y2_text, lower_text, upper_text):
    first = (binary64(x1_text), binary64(y1_text))
    second = (binary64(x2_text), binary64(y2_text))
    lower = binary64(lower_text)
    upper = binary64(upper_text)
    separation = abs(first[0] - second[0])

    if (
        kernel % 2 == 0
        and lower == -upper
        and first[0] == -second[0]
        and first[1] == second[1]
    ):
        return mp.mpf(0)

    def relevant(peak):
        center, width = peak
        distance = max(mp.mpf(0), lower - center, center - upper)
        return distance <= 8192 * width

    first_relevant = relevant(first)
    second_relevant = relevant(second)
    if first_relevant and second_relevant and separation > 8 * max(first[1], second[1]):
        if first[0] <= second[0]:
            left, right = first, second
        else:
            left, right = second, first
        split = (left[0] + right[0]) / 2
        if lower < split < upper:
            return segment(kernel, left, right, lower, split) + segment(kernel, right, left, split, upper)

    def interval_distance(peak):
        center, _ = peak
        return max(mp.mpf(0), lower - center, center - upper)

    first_much_narrower = first[1] < second[1] / 16
    second_much_narrower = second[1] < first[1] / 16
    if first_much_narrower:
        anchor, other = first, second
    elif second_much_narrower:
        anchor, other = second, first
    elif first_relevant and not second_relevant:
        anchor, other = first, second
    elif second_relevant and not first_relevant:
        anchor, other = second, first
    elif interval_distance(first) < interval_distance(second):
        anchor, other = first, second
    elif interval_distance(second) < interval_distance(first):
        anchor, other = second, first
    elif first[1] <= second[1]:
        anchor, other = first, second
    else:
        anchor, other = second, first
    return segment(kernel, anchor, other, lower, upper)


def cases():
    geometries = (
        ("0.1", "0.08", "0.35", "0.12", "-0.45", "0.6"),
        ("0.2", "1e-8", "0.20000003", "2e-8", "-0.4", "0.55"),
        ("-0.2", "1e-8", "0.3", "2e-8", "-0.4", "0.55"),
        ("-0.2", "1e-16", "0.3", "2e-16", "-0.4", "0.55"),
        ("0.2", "0.05", "-2e-8", "1e-8", "0", "1"),
        ("0.2", "0.05", "-2e-12", "1e-12", "0", "1"),
        ("0.2", "0.05", "-2e-16", "1e-16", "0", "1"),
        ("0.5000000000002", "1e-13", "0.5000000000015", "1e-13", "-0.5", "0.5"),
        ("0.500000000000002", "1e-15", "0.500000000000015", "1e-15", "-0.5", "0.5"),
        ("1.1", "1e-8", "0.2", "0.05", "-0.4", "0.55"),
        ("-0.2", "0.05", "0.2", "0.05", "-0.5", "0.5"),
        ("0.2", "0.05", "-5e-17", "1e-20", "0", "1"),
        ("0.2", "0.05", "-0.1", "1e-200", "0", "1"),
        ("0.5", "1e-100", "0.5000000000000001", "1e-100",
         "0.49999999999999994", "0.5000000000000002"),
        ("0.5", "1e-200", "0.5000000000000001", "1e-200",
         "0.49999999999999994", "0.5000000000000002"),
        ("-8.193e-197", "1e-200", "-9e-197", "1e-200", "0", "1"),
        ("0.2", "1e-10", "0.200001", "1e-10", "-0.5", "0.5"),
    )
    for kernel in range(1, 9):
        for geometry in geometries:
            yield (kernel, *geometry)


print("# m x1 y1 x2 y2 lower upper expected relative_tolerance absolute_tolerance")
for values in cases():
    value = reference(*values)
    print(*values, mp.nstr(value, 25), "2e-8", "1e-290")
