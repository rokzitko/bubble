# DMFT Bubble Integrals

`bubble` is a compact C++ program for evaluating real-frequency bubble
integrals from a local DMFT self-energy. For a set of commonly used band and
transport kernels, the band-energy integral is evaluated analytically, leaving
only a one-dimensional adaptive quadrature over frequency.

This is particularly useful in the low-temperature Fermi-liquid regime, where
the spectral function becomes narrow and a direct two-dimensional integration
over band energy and frequency can become slow or inaccurate.

The program is intended for single-band calculations with a local scalar
self-energy and no vertex corrections. It computes dimensionless transport
moments, not conductivities in SI units. Charge, spin/orbital degeneracy,
velocity, volume, lattice-spacing, and other model-dependent prefactors must be
supplied separately.

## Features

- Analytic band-energy integrals for flat, semicircular, Bethe lattice, and
  Gaussian kernels.
- Numerical support for arbitrary tabulated kernels \(\Phi(\epsilon)\).
- Spectral-function powers \(n=0,1,2,3\) and frequency moments \(\omega^o\).
- Linear, cubic-spline, or Akima interpolation of a tabulated self-energy.
- Adaptive Gauss-Kronrod integration through GSL.
- Fermi-derivative transport moments, occupied spectral moments, and
  frequency-resolved DOS output.
- Regression tests

## Quick Start

### Requirements

- A C++ compiler
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/)
- GNU Make and `pkg-config`
- Perl, only for `make check`

The complex error-function implementation needed by the Gaussian kernels is
included in `Faddeeva.cc` and `Faddeeva.hh`.

### Build

Build with the default C++ compiler and GSL flags supplied by `pkg-config`:

```bash
make
```

The compiler and flags can be overridden in the usual way, for example with
`make CXX=clang++`. Run `make clean` before switching compilers. Run the stable
24-case regression suite with:

```bash
make check
```

Remove generated build artifacts with `make clean`.

On systems providing the `intel/2022a` environment module (e.g. at IJS), the
supplied build script can be used instead:

```bash
./compile
```

The script purges the current module environment, loads `intel/2022a`, and
builds the program with `icpc`. It is therefore site-specific.

### Bethe-Lattice Example

The repository contains a real-frequency self-energy from a DMFT calculation
and corresponding independent reference results:

```bash
./bubble 5 2 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```

The expected raw transport moment is

```text
1.556196516510379
```

## Physical Quantity

For a retarded self-energy, the momentum-resolved Green's function and spectral
function are

$$
G^R_\epsilon(\omega) =
\frac{1}{\omega+\mu-\epsilon-\Sigma^R(\omega)},
\qquad
\rho_\epsilon(\omega) =
-\frac{1}{\pi}\mathrm{Im}\,G^R_\epsilon(\omega).
$$

The general dimensionless moments (described in the accompanying notes) are

$$
I_{mno} = \frac{1}{D}\int d\epsilon\,\Phi_m(\epsilon)
\int d\omega\,
\left(-\frac{\partial f}{\partial\omega}\right)
[D\rho_\epsilon(\omega)]^n
\left(\frac{\omega}{D}\right)^o,
$$

where

$$
f(\omega)=\frac{1}{1+e^{\omega/T}},
\qquad
-\frac{\partial f}{\partial\omega}
=\frac{1}{T[2+2\cosh(\omega/T)]}.
$$

The code first evaluates

$$
J_{mn}(\Omega)=\frac{1}{D}\int d\epsilon\,\Phi_m(\epsilon)
\left[-\frac{1}{\pi}\mathrm{Im}\,
\frac{D}{\Omega-\epsilon}\right]^n,
$$

with \(\Omega=\omega+\mu-\Sigma^R(\omega)\), and then performs the remaining
frequency integral numerically.

The current analytic kernels use \(D=1\) and \(\sigma=1\). Consequently,
temperature, chemical potential, frequency, band energy, and self-energy must
all be expressed in the same energy unit, normally the half-bandwidth. The
conventions \(k_B=\hbar=1\) are used.

## Transport Moments

The most common index combinations are:

| Application | `n` | `o` |
|---|---:|---:|
| DC conductivity moment | 2 | 0 |
| Thermoelectric moment | 2 | 1 |
| Thermal-transport moment | 2 | 2 |
| Hall-conductivity moment | 3 | 1 |
| Hall thermoelectric moment | 3 | 2 |

Physical observables are formed from these moments and appropriate prefactors.
For example, a Seebeck coefficient involves a ratio of the `o=1` and `o=0`
moments rather than the `o=1` result alone. The exact conversion depends on the
normalization of \(\Phi_m\), units, degeneracies, and the conventions used for
charge and current operators.

This program evaluates DC moments. Here \(\omega\) is the internal fermionic
frequency, not an external optical probe frequency.

## Built-in Kernels

For the analytic kernels, \(D=\sigma=1\):

| `m` | Kernel \(\Phi_m(\epsilon)\) | Domain | Description |
|---:|---|---|---|
| 0 | User-supplied table | Table interval | Generic numerical kernel |
| 1 | \(1\) | \([-1,1]\) | Flat, even |
| 2 | \(\epsilon\) | \([-1,1]\) | Flat, odd |
| 3 | \(\sqrt{1-\epsilon^2}\) | \([-1,1]\) | Semicircular, even |
| 4 | \(\epsilon\sqrt{1-\epsilon^2}\) | \([-1,1]\) | Semicircular, odd |
| 5 | \((1-\epsilon^2)^{3/2}\) | \([-1,1]\) | Bethe transport, even |
| 6 | \(\epsilon(1-\epsilon^2)^{3/2}\) | \([-1,1]\) | Bethe transport, odd |
| 7 | \(e^{-2\epsilon^2}\) | \(( -\infty,\infty )\) | Gaussian, even |
| 8 | \(\epsilon e^{-2\epsilon^2}\) | \(( -\infty,\infty )\) | Gaussian, odd |

These kernels do not include lattice normalization factors. For example, the
Bethe-lattice benchmark applies the factor \(2/\pi\) externally. Odd kernels
are useful in Hall and occupied-energy moments.

For `m=0`, the kernel is read from a whitespace-separated two-column file:

```text
epsilon_0  Phi(epsilon_0)
epsilon_1  Phi(epsilon_1)
...
```

The energy grid must be strictly increasing and should contain at least five
points for GSL Akima interpolation. The energy integration is restricted to the
table interval. Use `-p FILE` to select the table.

## Command Line

```text
bubble [options] m n o T mu resigma.dat imsigma.dat
```

The positional arguments are:

| Argument | Meaning |
|---|---|
| `m` | Kernel index, `0` through `8` |
| `n` | Power of the spectral function |
| `o` | Power of frequency in the outer integral |
| `T` | Temperature in the chosen energy unit |
| `mu` | Chemical potential |
| `resigma.dat` | Table of \(\omega,\mathrm{Re}\,\Sigma(\omega)\) |
| `imsigma.dat` | Table of \(\omega,\mathrm{Im}\,\Sigma(\omega)\) |

Options:

| Option | Meaning |
|---|---|
| `-v` | Print parameters, numerical result, and outer integration error estimate |
| `-i I` | Self-energy interpolation: `1` linear, `2` cubic spline, `3` Akima |
| `-k K` | GSL rule: `1` through `6` select 15, 21, 31, 41, 51, or 61 points |
| `-a A` | Absolute integration tolerance, default `1e-7` |
| `-r R` | Relative integration tolerance, default `1e-8` |
| `-c C` | Frequency cutoff in units of temperature, default `15` |
| `-p FILE` | Kernel table for `m=0`, default `Phi.dat` |
| `-e E` | Multiply a tabulated kernel by \(\epsilon^E\), default `0` |
| `-f` | Use the Fermi function instead of its derivative |
| `-d` | Skip the frequency integral and write `dos.dat`; requires `m=0` |

Without `-v`, normal integration mode prints one number to standard output.
This makes the executable convenient to call from shell, Python, Julia, or
other DMFT post-processing workflows.

## Self-Energy Input

The real and imaginary parts are supplied as separate two-column tables:

```text
# resigma.dat              # imsigma.dat
omega_0  ReSigma_0         omega_0  ImSigma_0
omega_1  ReSigma_1         omega_1  ImSigma_1
...                        ...
```

Both files should have the same number of rows and identical, strictly
increasing frequency grids. The retarded convention
\(\mathrm{Im}\,\Sigma^R\le 0\) is required.

The following numerical policies are important when interpreting results:

- Input values with \(\mathrm{Im}\,\Sigma>-10^{-8}\) are replaced by
  \(-10^{-8}\), including positive noncausal values.
- Outside the input interval, `ReSigma` is held at its nearest endpoint and
  `ImSigma` is set to \(-10^{-10}\).
- The default Fermi-derivative integration interval is
  \([-15T,15T]\). Increase `-c` if broader thermal tails matter.

Interpolation mode `-i 1` is the conservative default for noisy DMFT data.
Higher-order interpolation can be useful for smooth, well-resolved
self-energies but may introduce artifacts when the grid is sparse.

## Custom Kernels and Special Modes

### Additional Energy Powers

For `m=0`, `-e E` replaces the supplied kernel by

$$
\Phi(\epsilon)\longrightarrow\epsilon^E\Phi(\epsilon).
$$

For example:

```bash
./bubble -e 1 -p DOS.dat \
  0 1 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```

### Occupied Spectral Moments

`-f` replaces \(-\partial f/\partial\omega\) by the Fermi function itself:

$$
\int_{\omega_{\min}}^{CT}d\omega\,
f(\omega)\,\omega^oJ_{mn}(\omega).
$$

This is useful for quantities such as occupied spectral moments and
kinetic-energy checks. For example:

```bash
./bubble -f 4 1 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```


### Frequency-Resolved DOS

`-d` bypasses the outer frequency integration and writes 10,001 pairs

$$
\{\omega_i,J_{0n}(\omega_i)\}
$$

over the self-energy input interval to lowercase `dos.dat`:

```bash
./bubble -d -p DOS.dat \
  0 1 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```

The output is a physical interacting DOS only for `n=1` when the supplied
kernel is itself a noninteracting DOS. The positional `T` and `o` arguments are
still syntactically required but are not used in this mode. `-d` overwrites
`dos.dat` and cannot be combined with `-f`; it can be combined with `-e`.

## Numerical Method

- The analytic \(J_{mn}\) expressions were generated from Mathematica
  calculations documented under `notes/`.
- The remaining frequency integral uses adaptive GSL `qag` quadrature.
- The default rule is the 15-point Gauss-Kronrod rule with a workspace of
  1,000 intervals.
- A custom `m=0` kernel introduces a nested adaptive energy integral.
- Self-energy interpolation is selectable, while custom-kernel interpolation
  is always Akima.
- In verbose mode, the reported error is GSL's estimate for the outer
  frequency quadrature. The custom-kernel inner error estimate is not reported.

## Validation

The repository contains two complementary sets of checks.

### Legacy Regression Data

`regression/` contains constant and linearly varying self-energies together
with historical reference values over a range of temperatures and kernel
indices. A representative suite is run with:

```bash
cd regression
./run_tests 1
```

The current implementation passes 372 of the 384 historical comparisons. The
12 remaining cases are confined to the clean-limit `n=3` and suite `22` cases
marked as problematic in `regression/README`.

### Bethe-Lattice DMFT Benchmark

`Bethe_lattice_test/` contains a real-frequency DMFT self-energy and independent
Mathematica reference values for conductivity, resistivity, thermopower,
thermal conductivity, Lorenz ratio, and \(ZT\). With the documented
normalization, the current executable agrees with these references to about
\(10^{-5}\) relative accuracy. The independent occupied-energy check agrees to
better than \(10^{-7}\) relative accuracy.

Some scripts in this directory belong to the original larger DMFT/NRG workflow
and require tools that are not distributed here. The `bubble` invocations and
the checked-in input/reference data are sufficient for the examples above.

## Precision

For precision-sensitive work, vary `-a`, `-r`, `-c`, `-k`, and the input grid,
and confirm that the physical result is stable.

## Method Notes

The original derivation, kernel definitions, and implementation discussion are
available in:

- [`notes/bubble.pdf`](notes/bubble.pdf)
- [`notes/bubble.tex`](notes/bubble.tex)

The invocation section of the 2017 notes describes an earlier
interface. The command line documented in this README and printed by the
current executable is authoritative for version 1.4.

## References

1. A. Khurana, "Electrical conductivity in the infinite-dimensional Hubbard
   model," *Physical Review Letters* **64**, 1990 (1990),
   [doi:10.1103/PhysRevLett.64.1990](https://doi.org/10.1103/PhysRevLett.64.1990).
2. L.-F. Arsenault and A.-M. S. Tremblay, "Transport functions for hypercubic
   and Bethe lattices," *Physical Review B* **88**, 205109 (2013),
   [doi:10.1103/PhysRevB.88.205109](https://doi.org/10.1103/PhysRevB.88.205109).
3. G. Bevilacqua, "Some integrals related to the Fermi function,"
   [arXiv:1303.6206](https://arxiv.org/abs/1303.6206).

If this software contributes to published work, please cite the repository
release together with the method reference most relevant to the calculation.
A dedicated software DOI is not currently available.

## Authors

- **Amina Alic**: original implementation
- **Rok Žitko**: method notes and subsequent development

Issues and pull requests that improve numerical robustness, validation, or
support for additional lattice transport functions are welcome.

## License

The original project code and data are released under the
[GNU General Public License, version 3 or later](https://www.gnu.org/licenses/gpl-3.0.html).

`Faddeeva.cc` and `Faddeeva.hh` are by Steven G. Johnson, Massachusetts
Institute of Technology, and are distributed under the MIT license included in
their source headers.
