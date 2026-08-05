#!/usr/bin/env python3
"""Measure the LBFAST capillary-wave frequency from zero crossings."""

from __future__ import print_function

import argparse
import math
from pathlib import Path


def zero_crossings(times, signal):
    roots = []
    for index in range(1, len(signal)):
        y0, y1 = signal[index - 1], signal[index]
        if y0 == 0.0:
            roots.append(times[index - 1])
        elif y0 * y1 < 0.0:
            fraction = abs(y0) / (abs(y0) + abs(y1))
            roots.append(times[index - 1] + fraction * (times[index] - times[index - 1]))
    return roots


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("data", nargs="?", type=Path, default=Path("capillary_wave.dat"))
    args = parser.parse_args()

    rows = []
    for line in args.data.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 7:
            rows.append(tuple(float(value) for value in fields))
    if len(rows) < 10:
        raise SystemExit("Not enough capillary-wave samples")

    times = [row[0] for row in rows]
    signal = [row[1] for row in rows]
    roots = zero_crossings(times, signal)
    if len(roots) < 3:
        raise SystemExit("At least three zero crossings are required")
    half_periods = [roots[i] - roots[i - 1] for i in range(1, len(roots))]
    mean_half_period = sum(half_periods) / len(half_periods)
    omega = math.pi / mean_half_period
    omega_theory = rows[-1][6]
    error = omega / omega_theory - 1.0
    max_mass_error = max(abs(row[4]) for row in rows)

    print("samples={}".format(len(rows)))
    print("zero_crossings={}".format(len(roots)))
    print("omega_theory={:.12e}".format(omega_theory))
    print("omega_numerical={:.12e}".format(omega))
    print("frequency_error_percent={:.6f}".format(100.0 * error))
    print("period_numerical={:.12e}".format(2.0 * mean_half_period))
    print("max_mass_error={:.12e}".format(max_mass_error))
    if len(rows[0]) >= 9:
        eta_rates = [row[7] for row in rows[1:]]
        velocity_modes = [row[8] for row in rows[1:]]
        velocity_norm = sum(value * value for value in velocity_modes)
        if velocity_norm > 0.0:
            gain = sum(rate * velocity for rate, velocity in zip(eta_rates, velocity_modes)) / velocity_norm
            residual_norm = math.sqrt(sum(
                (rate - velocity) ** 2 for rate, velocity in zip(eta_rates, velocity_modes)
            ) / velocity_norm)
            print("kinematic_gain_deta_dt_over_vn={:.12e}".format(gain))
            print("kinematic_relative_rms_error={:.12e}".format(residual_norm))


if __name__ == "__main__":
    main()
