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
- Numerical support for arbitrary tabulated kernels $\Phi(\epsilon)$.
- Built-in spectral-function powers $n=0,1,2,3$, arbitrary nonnegative
  powers for tabulated kernels, and frequency moments $\omega^o$.
- Linear, cubic-spline, or Akima interpolation of a tabulated self-energy.
- Adaptive Gauss-Kronrod integration through GSL.
- Fermi-derivative transport moments, finite-frequency optical conductivity,
  occupied spectral moments, and frequency-resolved DOS output.
- DC, optical-limit, and independent Bethe-lattice regression tests.
- Clean-limit kernel checks against independent high-precision quadrature.

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
`make CXX=clang++`. Run `make clean` before switching compilers. Run the default
DC regression and optical validation suite with:

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
-\frac{1}{\pi}\mathrm{Im}\ G^R_\epsilon(\omega).
$$

The general dimensionless moments (described in the accompanying notes) are

$$
I_{mno} = \frac{1}{D}\int d\epsilon\ \Phi_m(\epsilon)
\int d\omega
\ \left(-\frac{\partial f}{\partial\omega}\right)
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
J_{mn}(z)=\frac{1}{D}\int d\epsilon\ \Phi_m(\epsilon)
\left[-\frac{1}{\pi}\mathrm{Im}
\ \frac{D}{z-\epsilon}\right]^n,
$$

with $z=\omega+\mu-\Sigma^R(\omega)$, and then performs the remaining
frequency integral numerically.

The current analytic kernels use $D=1$. Consequently,
temperature, chemical potential, frequency, band energy, and self-energy must
all be expressed in the same energy unit, normally the half-bandwidth. The
conventions $k_B=\hbar=1$ are used.

### Clean-Limit Kernels

The phrase *clean-limit kernel* refers to $J_{mn}(z)$ evaluated when the
single-particle scattering rate is very small. It is an evaluation regime, not
an additional transport function $\Phi_m$. Writing

$$
z=x+iy,
\qquad
x=\omega+\mu-\mathrm{Re}\ \Sigma^R(\omega),
\qquad
y=-\mathrm{Im}\ \Sigma^R(\omega)>0,
$$

the spectral factor for $D=1$ is the Lorentzian

$$
A_y(x-\epsilon)
=\frac{y}{\pi[(x-\epsilon)^2+y^2]}.
$$

The clean limit is $y\to0^+$, where this Lorentzian becomes increasingly
narrow and tends to $\delta(x-\epsilon)$. If $x$ is strictly inside a band
and $\Phi_m$ is smooth there, the leading behavior is

$$
J_{m1}(x+iy)\longrightarrow\Phi_m(x),
\qquad
J_{m2}(x+iy)\sim\frac{\Phi_m(x)}{2\pi y},
\qquad
J_{m3}(x+iy)\sim\frac{3\Phi_m(x)}{8\pi^2y^2}.
$$

Thus the `n=2` and `n=3` kernels become large even though their energy width
shrinks. Outside a finite band they instead vanish as $O(y^n)$, while at a
band edge the behavior depends on how $\Phi_m$ approaches zero. These very
different scales make direct quadrature and cancellation-prone closed formulas
unreliable. The built-in clean-limit implementation uses factored finite-band
expressions, exterior moment expansions, and conditioned Gaussian/Faddeeva
identities to retain accuracy in all three regions. The `-s` option sets the
minimum linewidth used when the supplied self-energy is even cleaner.

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
normalization of $\Phi_m$, units, degeneracies, and the conventions used for
charge and current operators.

Without `-O`, the program evaluates these DC moments. Here $\omega$ is the
internal fermionic frequency.

## Finite-Frequency Optical Conductivity

With `-O OMEGA`, `bubble` evaluates the two-spectral-function moment

$$
I_{mo}(\Omega)=
\int_{-CT-\Omega}^{CT}d\omega
\ \frac{f(\omega)-f(\omega+\Omega)}{\Omega}\ \omega^o
K_m[z(\omega),z(\omega+\Omega)],
$$

where

$$
K_m(z_1,z_2)=\int d\epsilon\ \Phi_m(\epsilon)
\left[-\frac{1}{\pi}\mathrm{Im}\frac{1}{z_1-\epsilon}\right]
\left[-\frac{1}{\pi}\mathrm{Im}\frac{1}{z_2-\epsilon}\right],
\qquad
z(\omega)=\omega+\mu-\Sigma^R(\omega).
$$

The optical frequency $\Omega$ is external and distinct from the integration
frequency $\omega$. The self-energy is evaluated independently at
$\omega$ and $\omega+\Omega$. The frequency quadrature is split at
$-\Omega$ and zero. It therefore requires self-energy data over approximately
$[-CT-\Omega,CT+\Omega]$; unless `-q` is used, the program warns before
using its normal out-of-table extrapolation policy.

Optical mode requires `n=2` and a finite `OMEGA >= 0`. It cannot be combined
with `-f` or `-d`. At `OMEGA=0`, the existing Fermi derivative and DC kernel are
used exactly. For small positive `OMEGA`, the finite-frequency result therefore
approaches the DC value without a subtractive Fermi-function cancellation.
The physical longitudinal optical conductivity uses `o=0`; other powers remain
available as generalized moments.

For example:

```bash
./bubble -O 0.1 5 2 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```

This is the regular absorptive bubble contribution. It does not add a separate
Drude delta function or diamagnetic term, vertex corrections, or model-specific
normalization factors.

## Built-in Kernels

For the analytic kernels, $D=\sigma=1$:

| `m` | Kernel $\Phi_m(\epsilon)$ | Domain | Description |
|---:|---|---|---|
| 0 | User-supplied table | Table interval | Generic numerical kernel |
| 1 | $1$ | $[-1,1]$ | Flat, even |
| 2 | $\epsilon$ | $[-1,1]$ | Flat, odd |
| 3 | $\sqrt{1-\epsilon^2}$ | $[-1,1]$ | Semicircular, even |
| 4 | $\epsilon\sqrt{1-\epsilon^2}$ | $[-1,1]$ | Semicircular, odd |
| 5 | $(1-\epsilon^2)^{3/2}$ | $[-1,1]$ | Bethe transport, even |
| 6 | $\epsilon(1-\epsilon^2)^{3/2}$ | $[-1,1]$ | Bethe transport, odd |
| 7 | $e^{-2\epsilon^2}$ | $( -\infty,\infty )$ | Gaussian, even |
| 8 | $\epsilon e^{-2\epsilon^2}$ | $( -\infty,\infty )$ | Gaussian, odd |

These kernels do not include lattice normalization factors. For example, the
Bethe-lattice benchmark applies the factor $2/\pi$ externally. Odd kernels
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
| `resigma.dat` | Table of $\omega,\mathrm{Re}\ \Sigma(\omega)$ |
| `imsigma.dat` | Table of $\omega,\mathrm{Im}\ \Sigma(\omega)$ |

Options:

| Option | Meaning |
|---|---|
| `-v` | Print parameters, numerical result, and outer integration error estimate |
| `-q` | Suppress non-fatal warnings without suppressing results or errors |
| `-i I` | Self-energy interpolation: `1` linear, `2` cubic spline, `3` Akima |
| `-k K` | GSL `qag` rule: `1` through `6` select 15, 21, 31, 41, 51, or 61 points |
| `-a A` | Absolute integration tolerance, default `1e-7` |
| `-r R` | Relative integration tolerance, default `1e-8` |
| `-c C` | Frequency cutoff in units of temperature, default `15` |
| `-s S` | Positive clipping floor for $-\mathrm{Im}\ \Sigma$, default `1e-8` |
| `-M M` | Restrict epsilon to a half-width `M` around the interacting Fermi level; default `0` is unrestricted |
| `-p FILE` | Kernel table for `m=0`, default `Phi.dat` |
| `-e E` | Multiply a tabulated kernel by $\epsilon^E$, default `0` |
| `-f` | Use the Fermi function instead of its derivative |
| `-d` | Skip the frequency integral and write `dos.dat`; requires `m=0` |
| `-O OMEGA` | External optical frequency; requires `OMEGA >= 0` and `n=2` |

Without `-v`, normal integration mode prints one number to standard output.
This makes the executable convenient to call from shell, Python, Julia, or
other DMFT post-processing workflows.

`-q` suppresses only non-fatal warning messages. It does not hide the numerical
result, fatal diagnostics, or the additional output requested by `-v`.

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
$\mathrm{Im}\ \Sigma^R\le 0$ is required.

The following numerical policies are important when interpreting results:

- Input values with $\mathrm{Im}\ \Sigma>-S$ are replaced by $-S$,
  including positive noncausal values. The `-s S` option selects this floor and
  defaults to $S=10^{-8}$.
- Outside the input interval, `ReSigma` is held at its nearest endpoint and
  `ImSigma` is set to $-10^{-10}$.
- The default Fermi-derivative integration interval is
  $[-15T,15T]$. Increase `-c` if broader thermal tails matter.
- Optical mode extends the lower endpoint by `OMEGA` and evaluates the shifted
  propagator up to `cutoff*T+OMEGA`.

Interpolation mode `-i 1` is the conservative default for noisy DMFT data.
Higher-order interpolation can be useful for smooth, well-resolved
self-energies but may introduce artifacts when the grid is sparse.

## Restricted Epsilon Window

With `-M M`, for `M > 0`, every band-energy integral is restricted by

$$
-M < \mu-\epsilon-\mathrm{Re}\ \Sigma^R(0) < M.
$$

Equivalently, the fixed epsilon interval is

$$
\epsilon_F-M < \epsilon < \epsilon_F+M,
\qquad
\epsilon_F=\mu-\mathrm{Re}\ \Sigma^R(0).
$$

The value of $\mathrm{Re}\ \Sigma^R(0)$ is obtained with the interpolation
selected by `-i`. Therefore a positive `M` requires the self-energy input grid
to contain $\omega=0$; the program reports an error rather than extrapolating
the Fermi-level center. `M` uses the same energy units as the bandwidth,
temperature, chemical potential, and self-energy.

The restriction applies to DC and optical calculations, occupied moments with
`-f`, tabulated `m=0` kernels, and DOS output with `-d`. The requested interval
is intersected with the kernel's natural domain: `[-1,1]` for built-in finite
bands and the tabulated interval for `m=0`. An empty intersection gives zero.
If the interval completely covers a finite natural domain, the existing
full-band evaluator is retained. Every positive finite `M` truncates the
Gaussian kernels because their natural domain is the full real line.

The strict endpoint inequalities do not affect these ordinary integrals, so
the numerical quadrature includes the endpoints. The value `M=0` is a sentinel
for an infinite window, not a zero-width interval. Both an omitted `-M` and an
explicit `-M 0` follow the previous unrestricted code paths exactly.
If a positive `M` is smaller than floating-point resolution at the resolved
center, so that the stored bounds do not lie on opposite sides of
$\epsilon_F$, the program reports an error instead of using a one-sided
window. It likewise rejects bounds whose rounding erases a nonzero resolved
center and would make the stored interval spuriously symmetric.

For example, this restricts the Bethe transport integral to a half-width of
`0.25` around the interacting Fermi level:

```bash
./bubble -M 0.25 5 2 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```

Restricted built-in kernels require nested numerical epsilon integration and
can therefore be slower than their unrestricted analytic counterparts. The
implementation removes narrow Lorentzian peaks with a tangent transformation
or a logarithmic distance transformation for exterior peaks. Finite-band edge
charts and shifted transformed coordinates avoid cancellation at sub-ULP
linewidths and intervals. The outer frequency integral is split at
linewidth-aware crossings and nearby extrema of the window boundaries.
Smooth tangencies use a curvature-scaled neighborhood, while slope reversals
at interpolation knots receive explicit one-sided neighborhoods.
Restricted inner integrals use adaptive 61-point quadrature;
breakpoint-based `qagp` integrations do not use the `-k` selection. A
roundoff-limited result is accepted only within the requested error bound or
after agreement with independent 64- and 128-point rules, and emits a warning
unless `-q` is active.

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

`-f` replaces $-\partial f/\partial\omega$ by the Fermi function itself:

$$
\int_{\omega_{\min}}^{CT}d\omega
\ f(\omega)\ \omega^oJ_{mn}(\omega).
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

The output is a physical interacting DOS (lattice spectral function) only for `n=1` when the supplied
kernel is itself a noninteracting DOS. The positional `T` and `o` arguments are
still syntactically required but are not used in this mode. `-d` overwrites
`dos.dat` and cannot be combined with `-f`; it can be combined with `-e`.

## Numerical Method

- The $n=0,1$ analytic $J_{mn}$ expressions originate from Mathematica
  calculations documented under `notes/`.
- Built-in $n=2,3$ kernels use cancellation-free finite-band formulas,
  exterior moment expansions, and conditioned Faddeeva identities. Gaussian
  transition cases have a direct positive-quadrature fallback.
- Built-in optical kernels use analytic Hilbert transforms and stable divided
  differences rather than generated two-frequency expressions.
- Restricted kernels use bounded, peak-removing epsilon quadrature; separated
  optical peaks are integrated with independent local regions and exterior
  peaks use logarithmic distance coordinates.
- The remaining frequency integral uses adaptive GSL `qag` quadrature.
- The default rule is the 15-point Gauss-Kronrod rule with a workspace of
  1,000 intervals.
- A custom `m=0` kernel introduces a nested adaptive energy integral with
  linewidth-aware breakpoints around both spectral peaks.
- Self-energy interpolation is selectable, while custom-kernel interpolation
  is always Akima.
- In verbose mode, the reported error is GSL's estimate for the outer
  frequency quadrature. Nested restricted and custom-kernel inner uncertainty
  is validated separately but is not included in that printed estimate.

## Validation

The repository contains complementary DC, clean-kernel, restricted-window, and
optical checks.
`make check` runs legacy suite 1, the audited clean-limit legacy suites, the
pointwise kernel suites, restricted-window command tests, and both optical
checks without invoking Mathematica.

### Legacy Regression Data

`regression/` contains constant and linearly varying self-energies together
with historical reference values over a range of temperatures and kernel
indices. A representative suite is run with:

```bash
cd regression
./run_tests 1
```

The current implementation passes all 384 legacy suite rows. References in
clean-limit suites `5` and `22*` were independently regenerated for the finite
integration interval used by `run_tests`; parity-forced zeros were preserved
exactly.

### Clean-Kernel Regression

`regression/clean_kernel_references.dat` contains pointwise values for all
built-in kernels with `n=2,3`. It covers interior points, exact and nearby band
edges, clean exterior points, Gaussian pole/tail transitions, broad
linewidths, and exact parity. Most frozen values were generated at 100 decimal
digits after the substitution

$$
\epsilon=x+y\tan\theta,
$$

which removes the narrow Lorentzian peak from the numerical quadrature. Finite
exterior points without a peak in the band use direct high-precision energy
quadrature. To print the regenerated table, install the optional Python package
`mpmath` and run:

```bash
python3 regression/generate_clean_kernel_references.py
```

The bounded single- and two-peak values are stored in
`restricted_kernel_references.dat` and `restricted_optical_references.dat`.
Their corresponding `generate_restricted_*_references.py` scripts use
180-250-decimal-digit quadrature and controlled clean-limit asymptotics. They
include exact odd symmetry, sub-ULP
transformed spans, finite-band edges, adjacent-double optical centers, huge
Gaussian windows, and exterior linewidths down to `1e-200`.
`run_window_tests` additionally checks command parsing, exact `M=0` behavior,
full and empty intersections, tabulated kernels, optical and occupied moments,
`n=0` and `n=1` DOS output, and narrow outer-window crossings.

### Optical Regression

`regression/run_optical_tests` checks the exact `-O 0` DC path, small-frequency
convergence, independent finite-frequency values for the built-in kernel
families, a large-frequency shifted-range case, a tabulated `m=0` kernel, and
invalid option combinations. It can be run directly with:

```bash
./regression/run_optical_tests
```

### Bethe-Lattice DMFT Benchmark

`Bethe_lattice_test/` contains a real-frequency DMFT self-energy and independent
Mathematica reference values for conductivity, resistivity, thermopower,
thermal conductivity, Lorenz ratio, and $ZT$. With the documented
normalization, the current executable agrees with these references to about
$10^{-5}$ relative accuracy. The independent occupied-energy check agrees to
better than $10^{-7}$ relative accuracy.

The finite-frequency reference contains 50 optical frequencies. Generate
`cond.opt.geo.dat` on that exact mesh and compare it with the frozen Mathematica
output using:

```bash
./Bethe_lattice_test/cond.opt.geo
```

The script applies the benchmark's $2\pi$ normalization, `20T` cutoff, and
$10^{-16}$ self-energy clipping floor. It accepts
$|\Delta|\le10^{-12}+3\times10^{-5}|\mathrm{reference}|$, reports the largest
errors, and exits unsuccessfully on any mismatch. `cond.opt.geo.mma` is retained
only as provenance; validation reads `cond.opt.geo-mma.dat` and never launches
Mathematica.

Some scripts in this directory belong to the original larger DMFT/NRG workflow
and require tools that are not distributed here. The `bubble` invocations and
the checked-in input/reference data are sufficient for the examples above.

## Precision

For precision-sensitive work, vary `-a`, `-r`, `-c`, `-k`, `-s`, `-M`, and the
input grid, and confirm that the physical result is stable.

## Method Notes

The original derivation, kernel definitions, and implementation discussion are
available in:

- [`notes/bubble.pdf`](notes/bubble.pdf)
- [`notes/bubble.tex`](notes/bubble.tex)

The invocation section of the 2017 notes describes an earlier
interface. The command line documented in this README and printed by the
current executable is authoritative for version 1.6.

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
