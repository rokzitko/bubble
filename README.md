# DMFT Bubble Integrals

`bubble` is a compact C++ program for evaluating real-frequency bubble
integrals from a local DMFT self-energy. For a set of commonly used band and
transport kernels, the band-energy integral is evaluated with
conditioned analytic formulas and numerical fallbacks, leaving a
one-dimensional adaptive quadrature over frequency. Restricted and tabulated
kernels use a nested band-energy quadrature.

This code is particularly useful in the low-temperature Fermi-liquid regime, where
the spectral function becomes narrow and a direct two-dimensional integration
over band energy and frequency can become slow or inaccurate.

The program is intended for single-band calculations with a local scalar
self-energy and no vertex corrections. It computes dimensionless transport
moments. Charge, spin/orbital degeneracy, velocity, volume, lattice-spacing, 
and other model-dependent prefactors must be supplied separately.

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

## Applications

See [DMFT_NRG_KLM](https://github.com/rokzitko/DMFT_NRG_KLM) for an example of a dynamical mean
field theory (DMFT) calculation using numerical renormalization group (NRG) as a solver for
the Kondo lattice model (KLM).

## Quick Start

### Requirements

- A C++ compiler
- [GNU Scientific Library (GSL)](https://www.gnu.org/software/gsl/) 2.0 or newer
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
`make CXX=clang++`. Run `make clean` before switching compilers. Run the
regression, clean-limit, optical, and malformed-input validation suites with:

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
Note, with this compiler it is recommended to use the following coimpiler flags:

```text
CXXFLAGS='-O2 -std=c++11 -fp-model precise -no-ftz'
```

### Bethe-Lattice Example

The repository contains a real-frequency self-energy from a DMFT calculation
and corresponding independent reference results:

```bash
./bubble -E -q 5 2 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```

The expected raw transport moment is

```text
1.556196516510379
```

## Physical Quantity

For a retarded self-energy, define

$$
z(\omega)=x(\omega)+i\gamma(\omega)
=\omega+\mu-\Sigma^R_{\mathrm{eff}}(\omega),
\qquad
x(\omega)=\omega+\mu-\mathrm{Re}\thinspace{}\Sigma^R(\omega),
\qquad
\gamma(\omega)=-\mathrm{Im}\thinspace{}\Sigma^R_{\mathrm{eff}}(\omega)>0.
$$

The Green's function and spectral function at band energy $\epsilon$ are

$$
G^R_\epsilon(\omega)=\frac{1}{z(\omega)-\epsilon},
\qquad
A_\epsilon(z)=-\frac{1}{\pi}\mathrm{Im}\frac{1}{z-\epsilon}
=\frac{\gamma}{\pi[(x-\epsilon)^2+\gamma^2]}.
$$

For the selected energy domain $\mathcal E$, `bubble` first evaluates

$$
J_{mn}(z;\mathcal E)=
\int_{\mathcal E}d\epsilon\thinspace{}\Phi_m(\epsilon)[A_\epsilon(z)]^n.
$$

For $n=0$, this notation means the explicitly defined kernel mass
$J_{m0}=\int_{\mathcal E}\Phi_m(\epsilon)d\epsilon$, independent of $z$; it
does not rely on an $A^0$ convention. The ordinary transport mode then computes

$$
I^{\mathrm{DC}}_{mno}=\int_{-CT}^{CT}d\omega\thinspace{}
[-f'(\omega)]\thinspace{}\omega^oJ_{mn}[z(\omega);\mathcal E],
$$

where

$$
f(\omega)=\frac{1}{1+e^{\omega/T}},
\qquad
-f'(\omega)=\frac{1}{T[2+2\cosh(\omega/T)]},
$$

and $C$ is the cutoff selected by `-c`. With `-f`, the outer integral instead is

$$
I^{\mathrm{occ}}_{mno}=\int_{\omega_{\min}}^{CT}d\omega\thinspace{}
f(\omega)\thinspace{}\omega^oJ_{mn}[z(\omega);\mathcal E].
$$

These are the unmasked ranges. A positive spectral-frequency window intersects
them as described under [Spectral-Frequency Mask](#spectral-frequency-mask).

The integration frequency $\omega$ is measured relative to the Fermi level;
$\mu$ enters $z(\omega)$, not the Fermi function. The implementation has no
bandwidth parameter and uses $D=1$ throughout, including for a tabulated
`m=0` kernel. Temperature, chemical potential, frequency, band energy, and
self-energy must therefore be nondimensionalized consistently before use,
normally by the half-bandwidth. A custom table is not rescaled automatically.
The conventions $k_B=\hbar=1$ are used.

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
identities. Tabulated kernels remove the same narrow peaks with tangent or
logarithmic distance coordinates instead of integrating directly in epsilon.
The `-s` option sets the minimum linewidth used when the supplied self-energy
is even cleaner.

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
K_m(z_1,z_2;\mathcal E)=\int_{\mathcal E}d\epsilon\ \Phi_m(\epsilon)
\left[-\frac{1}{\pi}\mathrm{Im}\frac{1}{z_1-\epsilon}\right]
\left[-\frac{1}{\pi}\mathrm{Im}\frac{1}{z_2-\epsilon}\right],
\qquad
z(\omega)=\omega+\mu-\Sigma^R_{\mathrm{eff}}(\omega).
$$

The optical frequency $\Omega$ is external and distinct from the integration
frequency $\omega$. The power $\omega^o$ uses the unshifted incoming frequency,
not $\omega+\Omega$ or $\omega+\Omega/2$. The self-energy is evaluated
independently at $\omega$ and $\omega+\Omega$. The frequency quadrature is split
at $-\Omega$ and zero. Its two propagator arguments require the union

$$
[-CT-\Omega,CT]\cup[-CT,CT+\Omega]
=[-CT-\Omega,CT+\Omega].
$$

This is the unmasked range. A positive spectral-frequency window restricts both
propagator arguments as described under
[Spectral-Frequency Mask](#spectral-frequency-mask).

Optical mode requires `n=2` and a finite `OMEGA >= 0`. It cannot be combined
with `-f` or `-d`. At `OMEGA=0`, the exact DC `n=2` path is used, with
$K_m(z,z)=J_{m2}(z)$. For small positive `OMEGA`, the finite-frequency result
therefore approaches the DC value without a subtractive Fermi-function
cancellation.
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

The finite-band half-width and Gaussian scale are fixed by the displayed
dimensionless formulas:

| `m` | Kernel $\Phi_m(\epsilon)$ | Natural domain | $J_{m0}$ | Description |
|---:|---|---|---:|---|
| 0 | User-supplied table | Table interval | Spline integral | Generic numerical kernel |
| 1 | $1$ | $[-1,1]$ | $2$ | Flat, even |
| 2 | $\epsilon$ | $[-1,1]$ | $0$ | Flat, odd |
| 3 | $\sqrt{1-\epsilon^2}$ | $[-1,1]$ | $\pi/2$ | Semicircular, even |
| 4 | $\epsilon\sqrt{1-\epsilon^2}$ | $[-1,1]$ | $0$ | Semicircular, odd |
| 5 | $(1-\epsilon^2)^{3/2}$ | $[-1,1]$ | $3\pi/8$ | Bethe transport, even |
| 6 | $\epsilon(1-\epsilon^2)^{3/2}$ | $[-1,1]$ | $0$ | Bethe transport, odd |
| 7 | $e^{-2\epsilon^2}$ | $(-\infty,\infty)$ | $\sqrt{\pi/2}$ | Gaussian, even |
| 8 | $\epsilon e^{-2\epsilon^2}$ | $(-\infty,\infty)$ | $0$ | Gaussian, odd |

Kernels `m=1` through `6` vanish outside `[-1,1]`; Gaussian kernels `m=7,8`
have support on the full real line. Built-ins support `n=0,1,2,3`. A tabulated
`m=0` kernel supports any nonnegative `n`, although optical mode always requires
`n=2`.

No hidden density-of-states or transport normalization is included. For
example, the normalized Bethe benchmark transport function is
$\Phi_{\mathrm{benchmark}}=(2/\pi)(1-\epsilon^2)^{3/2}$. If a reference formula
uses $[\mathrm{Im}\thinspace{}G]^2$ rather than
$A^2=[-\mathrm{Im}\thinspace{}G/\pi]^2$, the combined conversion from the built-in `m=5`
result is $(2/\pi)\pi^2=2\pi$. Odd kernels are useful in Hall and
occupied-energy moments.

All tabulated inputs use the same strict line format. Each nonempty data line
must contain exactly two complete finite floating-point values separated by
whitespace. A `#` starts a comment, either on its own line or after the two
values; blank and comment-only lines are ignored. Malformed values, extra or
missing columns, overflow, and underflow are rejected with a filename and line
number.

For `m=0`, the kernel is read from a two-column file:

```text
epsilon_0  Phi(epsilon_0)
epsilon_1  Phi(epsilon_1)
...
```

The energy grid must be strictly increasing. The `-e E` transformation is
applied to the knots first, producing $\epsilon_i^E\Phi_i$, and that effective
table is then interpolated. Every effective value must remain finite. Tables
whose effective values are all
nonnegative use GSL Steffen interpolation and require at least three rows. This
monotonic cubic method cannot undershoot below the supplied values. Signed
effective kernels retain Akima interpolation and require at least five rows.
The self-energy option `-i 2` uses GSL's natural cubic spline. The energy
integration is restricted to the table interval. Use `-p FILE` to select the
table.

## Command Line

```text
bubble [options] m n o T mu resigma.dat imsigma.dat
```

Options must precede the seven positional arguments. The final seven tokens are
always interpreted positionally, so a negative finite chemical potential such
as `mu=-0.25` does not require a preceding `--`. Numeric arguments must be
complete finite values; malformed suffixes, overflow, and underflow are
rejected.

The positional arguments are:

| Argument | Meaning |
|---|---|
| `m` | Kernel index, `0` through `8` |
| `n` | Nonnegative power of the spectral function; built-in kernels support `0` through `3` |
| `o` | Nonnegative power of frequency in an outer integral; unused by `-d` |
| `T` | Positive temperature in the chosen energy unit; unused by `-d` |
| `mu` | Any finite chemical potential, including negative values |
| `resigma.dat` | Table of $\omega,\mathrm{Re}\ \Sigma(\omega)$ |
| `imsigma.dat` | Table of $\omega,\mathrm{Im}\ \Sigma(\omega)$ |

Options:

| Option | Meaning |
|---|---|
| `-v` | Print parameters, numerical result, and outer integration error estimate |
| `-q` | Suppress non-fatal warnings without suppressing results or errors |
| `-E`, `--ignore-integration-errors` | Continue after an integration failure when a finite partial result is available |
| `--allow-legacy-self-energy-extrapolation` | Explicitly permit the historical out-of-table self-energy model |
| `-i I` | Self-energy interpolation: `1` linear, `2` cubic spline, `3` Akima |
| `-k K` | GSL `qag` rule: `1` through `6` select 15, 21, 31, 41, 51, or 61 points; also used by transformed tabulated-kernel quadrature |
| `-a A` | Nonnegative absolute integration tolerance, default `1e-7` |
| `-r R` | Nonnegative relative integration tolerance, default `1e-8` |
| `-c C` | Positive frequency cutoff in units of temperature, default `30`; unused by `-d` |
| `-s S` | Positive in-table clipping floor for $-\mathrm{Im}\ \Sigma$, default `1e-8` |
| `-M M` | Restrict epsilon to a half-width `M` around the interacting Fermi level; default `0` is unrestricted |
| `-Mp MP` | Mask spectral functions outside the frequency half-width $M'$; default `0` is unrestricted |
| `-p FILE` | Kernel table for `m=0`, default `Phi.dat` |
| `-e E` | Multiply a tabulated kernel by $\epsilon^E$ for nonnegative `E`, default `0` |
| `-f` | Use the Fermi function instead of its derivative |
| `-d` | Skip the frequency integral and write `dos.dat`; requires `m=0` |
| `-O OMEGA` | External optical frequency; requires `OMEGA >= 0` and `n=2` |

### Output and Exit Status

Without `-v`, a normal DC or optical integration prints one number with 16
significant digits to standard output. With `-v`, the program prints the parsed
parameters followed by `Result` and the outer quadrature's `Error` estimate.
This makes the executable convenient to call from shell, Python, Julia, or
other DMFT post-processing workflows.

DOS mode `-d` writes `dos.dat` and normally prints no scalar result. If `n!=1`,
its nonfatal warning is currently written to standard output unless `-q` is
used. The executable has no `-h` or `--help` option; invoking it without all
seven positional arguments prints its usage and exits unsuccessfully.

Fatal diagnostics are written to standard error and return a nonzero status.
`-q` suppresses only nonfatal warnings; it does not hide numerical results,
fatal diagnostics, or output requested by `-v`.

By default, an adaptive integration failure is fatal and no numerical result is
printed. A restricted evaluator may accept a roundoff or singularity status
only after an independent quadrature check agrees within the requested bound.
`-E` downgrades other integration failures or invalid error estimates to a
warning when GSL supplied a finite partial result. Combining `-E -q` reproduces
the historical silent best-effort behavior. This is not an accuracy guarantee:
allocation failures and non-finite results remain fatal, and `-E` cannot detect
an inaccurate calculation for which GSL returned success.

## Self-Energy Input

The real and imaginary parts are supplied as separate two-column tables:

```text
# resigma.dat              # imsigma.dat
omega_0  ReSigma_0         omega_0  ImSigma_0
omega_1  ReSigma_1         omega_1  ImSigma_1
...                        ...
```

Both files must have the same number of rows and exactly identical, strictly
increasing frequency grids. Generate or copy the frequency column once rather
than independently rounding it in the two files. The selected interpolation
requires at least two rows for linear mode, three for cubic-spline mode, and
five for Akima mode. Physical retarded input should satisfy
$\mathrm{Im}\thinspace{}\Sigma^R\le0$, but the program corrects rather than rejects
positive values. For clipping floor $S$, it constructs

$$
\widetilde{\Sigma''}_i=\min(\Sigma''_i,-S),
\qquad
\Sigma''_{\mathrm{eff}}(\omega)=
\min\left(\mathop{\mathrm{interp}}\lbrace\widetilde{\Sigma''}_i\rbrace,-S\right),
$$

so the spectral linewidth obeys $\gamma=-\Sigma''_{\mathrm{eff}}\ge S$ even
if a higher-order spline overshoots between causal knots. Because `-s` changes
the linewidth, it can strongly affect clean-limit `n=2,3` moments.

For `n>0`, the self-energy tables must cover every frequency used by the
selected calculation. With the spectral-frequency mask disabled, coverage is
checked as the following exact closed interval before quadrature starts:

| Mode | Required self-energy interval |
|---|---|
| DC and `-O 0` | $[-CT,CT]$ |
| Optical `OMEGA > 0` | $[-(CT+\mathrm{OMEGA}),CT+\mathrm{OMEGA}]$ |
| Occupied `-f` | $[\omega_{\min},CT]$ |
| DOS `-d` | The table interval used for the output mesh |

An incomplete interval is fatal even with `-q` or `-E`. Extend both tables or
reduce `-c`, `-O`, or `-Mp` as appropriate. Calculations with `n=0` do not
evaluate the frequency-dependent self-energy and are exempt from this coverage
check. The test is exact and includes both endpoints, so a table endpoint even one
binary64 step short is insufficient. Both files remain mandatory for `n=0` and
are still parsed, grid-matched, and used to construct interpolants. A positive
`-M` also evaluates $\mathrm{Re}\thinspace{}\Sigma(0)$ and therefore requires the table to
contain zero; if the resulting epsilon intersection is empty, no additional
frequency coverage is required.
For positive $M'$, only the active frequencies described under
[Spectral-Frequency Mask](#spectral-frequency-mask) require coverage.

The following numerical policies are important when interpreting results:

The `-s S` option selects the in-table clipping floor and defaults to
$S=10^{-8}$. `--allow-legacy-self-energy-extrapolation` explicitly restores
the historical exterior model:

```math
\Sigma^R_{\mathrm{ext}}(\omega)=
\mathrm{Re}\thinspace{}\Sigma^R(\omega_{\min})-i10^{-10}
\quad(\omega<\omega_{\min}),
```

```math
\Sigma^R_{\mathrm{ext}}(\omega)=
\mathrm{Re}\thinspace{}\Sigma^R(\omega_{\max})-i10^{-10}
\quad(\omega>\omega_{\max}).
```

The fixed exterior width is not controlled by `-s`. The program warns once and
splits quadrature at every unshifted and optical-shifted table boundary; `-q`
suppresses this warning. The default Fermi-derivative integration interval is
$[-30T,30T]$; increase `-c` if broader thermal tails matter. Optical mode extends
the lower endpoint by `OMEGA` and evaluates the shifted propagator up to
`cutoff*T+OMEGA`.

Interpolation mode `-i 1` is the conservative default for noisy DMFT data.
Higher-order interpolation can be useful for smooth, well-resolved
self-energies but may introduce artifacts when the grid is sparse. The
post-interpolation imaginary-part floor preserves causality but can create
continuous floor plateaus with derivative kinks.

## Restricted Epsilon Window

The `-M M` option controls $M=M_\epsilon$, the epsilon-window half-width. It is
distinct from the spectral-frequency half-width $M'=M_\omega$ selected by
`-Mp MP`.

Define the interacting Fermi-level band energy

$$
\epsilon_F=\mu-\mathrm{Re}\thinspace{}\Sigma^R(0).
$$

The energy domain selected by `-M M` is

$$
\mathcal E_m(M)=
\begin{cases}
\mathcal E_m,&M=0,\\
\mathcal E_m\cap[\epsilon_F-M,\epsilon_F+M],&M>0,
\end{cases}
$$

with

$$
\mathcal E_0=[\epsilon_{\min},\epsilon_{\max}],\qquad
\mathcal E_{1\ldots6}=[-1,1],\qquad
\mathcal E_{7,8}=\mathbb R.
$$

The value of $\mathrm{Re}\ \Sigma^R(0)$ is obtained with the interpolation
selected by `-i`. Therefore a positive `M` requires the self-energy input grid
to contain $\omega=0$; the program reports an error rather than extrapolating
the Fermi-level center. `M` uses the same energy units as the bandwidth,
temperature, chemical potential, and self-energy.

The restriction applies to DC and optical calculations, occupied moments with
`-f`, tabulated `m=0` kernels, and DOS output with `-d`. It truncates only the
epsilon integral: it does not shorten the frequency interval, follow the
frequency-dependent center $x(\omega)$, or renormalize the retained part of the
kernel. An empty intersection gives zero. If the interval completely covers a
finite natural domain, the full-band evaluator is retained. Every positive
finite `M` truncates a Gaussian kernel because its natural domain is the full
real line.

The strict endpoint inequalities do not affect these ordinary integrals, so
the numerical quadrature includes the endpoints. The value `M=0` is a sentinel
for an infinite window, not a zero-width interval. Both an omitted `-M` and an
explicit `-M 0` select the full natural kernel domain and are exactly
equivalent.
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
Restricted built-in inner integrals use adaptive 61-point quadrature;
their breakpoint-based `qagp` integrations do not use the `-k` selection. A
roundoff-limited result is accepted only within the requested error bound or
after agreement with independent 64- and 128-point rules, and emits a warning
unless `-q` is active.

## Spectral-Frequency Mask

The `-Mp MP` option controls $M'=M_\omega$, a window in physical frequency.
This is different from $M=M_\epsilon$: `-M` changes the epsilon integration
domain, while `-Mp` leaves that domain unchanged and masks every spectral
function according to its own frequency argument. Both quantities use the same
energy units as temperature and the self-energy input.

Define

$$
\chi_{M'}(\nu)=
\begin{cases}
1,& |\nu|<M',\\
0,& |\nu|\ge M',
\end{cases}
$$

and replace each spectral function by

$$
A^{(M')}_\epsilon(\nu)=\chi_{M'}(\nu)A_\epsilon(\nu).
$$

Omitted `-Mp` and `-Mp 0` are equivalent and leave the
calculation unrestricted. The active-range
formulas below apply only when $M'>0$; for $M'=0$, the unmasked ranges given in
the preceding sections apply.

For positive $M'$, the mask intersects the existing thermal cutoff; it does not
replace `-c`. For DC moments with `n>0`, the active integration interval is

$$
\mathcal W_{\mathrm{DC}}=[-CT,CT]\cap(-M',M'),
$$

or, in terms of its limiting endpoints,

$$
L_{\mathrm{DC}}=\max(-CT,-M'),\qquad
U_{\mathrm{DC}}=\min(CT,M').
$$

For positive $C$, $T$, and $M'$, this interval is nonempty because it contains
zero. At `-O 0`, the same DC rule applies.

For an occupied moment selected by `-f`,

$$
\mathcal W_{\mathrm{occ}}
=[\omega_{\min},CT]\cap(-M',M'),
$$

with

$$
L_{\mathrm{occ}}=\max(\omega_{\min},-M'),\qquad
U_{\mathrm{occ}}=\min(CT,M').
$$

The frequency-window condition for a finite occupied contribution is

$$
L_{\mathrm{occ}} < U_{\mathrm{occ}}.
$$

At positive external optical frequency, the two spectral factors have
arguments $\omega$ and $\omega+\Omega$. Both must be active, so

$$
\begin{aligned}
\mathcal W_{\mathrm{opt}}
={}&[-CT-\Omega,CT]
\cap(-M',M')\\
&\cap(-\Omega-M',-\Omega+M'),
\end{aligned}
$$

The lower endpoint of an intersection is the maximum of its lower bounds, and
the upper endpoint is the minimum of its upper bounds. Here the limiting
endpoints are

$$
L_{\mathrm{opt}}
=\max(-CT-\Omega,-M',-\Omega-M')
=\max(-CT-\Omega,-M'),
$$

$$
U_{\mathrm{opt}}
=\min(CT,M',-\Omega+M')
=\min(CT,-\Omega+M').
$$

Here the simplification uses the supported $\Omega\ge0$, for which the overlap
of the two mask intervals is

$$
(-M',M')\cap(-\Omega-M',-\Omega+M')
=(-M',-\Omega+M').
$$

The frequency-window condition for a finite optical contribution is

$$
L_{\mathrm{opt}} < U_{\mathrm{opt}}.
$$

For positive $C$, $T$, and $M'$ and the supported $\Omega\ge0$, this is
equivalent to

$$
\Omega<2M'.
$$

The mask applies to every positive spectral power, including tabulated `m=0`
kernels. It does not affect `n=0`, whose kernel mass contains no spectral
function. In DOS mode, the output mesh and its 10,001 rows are unchanged, but
values with `n>0` and $|\omega|\ge M'$ are written as exact zeros.

When `-M` and `-Mp` are used together, their restrictions apply independently:
$M=M_\epsilon$ first selects the epsilon domain and $M'=M_\omega$ selects the
active frequencies. For `n>0`, a finite contribution requires both retained
epsilon support and frequency limits satisfying $L<U$; for `n=0`, $M'$ has no
effect.

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

For `n>0`, tabulated DC and optical kernels use the same bounded tangent and
logarithmic peak-removing coordinates as restricted built-in kernels. Their
adaptive pieces honor `-a`, `-r`, and `-k`, with the absolute tolerance divided
among pieces. The outer frequency quadrature is also split where either
effective-frequency argument crosses a table endpoint. For `n=0`, the result
does not depend on the propagator and is evaluated directly as the exact
integral of the selected interpolation spline, then cached for repeated DOS or
outer-integrand calls.

### Occupied Spectral Moments

`-f` replaces $-\partial f/\partial\omega$ by the Fermi function itself:

$$
\int_{\omega_{\min}}^{CT}d\omega
\ f(\omega)\ \omega^oJ_{mn}[z(\omega)].
$$

This is useful for quantities such as occupied spectral moments and
kinetic-energy checks. Without a positive `-Mp`, the upper endpoint `CT` must
lie in the self-energy table and the lower endpoint is the table's own
`omega_min`. A positive spectral-frequency window applies the range and
coverage rules stated above. For example:

```bash
./bubble -f 4 1 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```


### Frequency-Resolved DOS

`-d` bypasses the outer frequency integration and writes 10,001 pairs

$$
\lbrace\omega_i,J_{0n}(\omega_i)\rbrace
$$

over the self-energy input interval to lowercase `dos.dat`:

```bash
./bubble -d -p DOS.dat \
  0 1 0 0.27 0 \
  Bethe_lattice_test/resigma.dat \
  Bethe_lattice_test/imsigma.dat
```

The output is a physical interacting DOS (lattice spectral function) only for
`n=1` when the supplied kernel is itself a noninteracting DOS. The positional
`T` and `o` arguments are still syntactically required but are not used in this
mode, and neither is `-c`. They must be complete finite numeric values but may
have either sign. A positive `-Mp` still masks the spectral values as described
above. `-d` overwrites `dos.dat` and cannot be combined with `-f`; it can be
combined with `-e`. Rows are written with enough significant digits to
round-trip `double` values. The complete result is first written to a checked
temporary file in the same directory and then atomically replaces `dos.dat`, so
an existing output is preserved if calculation or output fails. For `n=0`, the
kernel mass is independent of Sigma and the output loop does not evaluate the
self-energy spline.

## Numerical Method

The analytic starting point is the Hilbert transform

$$
H_m(z)=\int_{\mathcal E_m}
\frac{\Phi_m(\epsilon)}{z-\epsilon}\thinspace{}d\epsilon,
\qquad
J_{m1}(z)=-\frac{\mathrm{Im}\thinspace{}H_m(z)}{\pi}.
$$

Full finite-band `n=2,3` kernels use factored long-double expressions near the
band edges and moment expansions outside the band. Gaussian kernels use
conditioned Faddeeva/Hilbert identities together with the exponentially small
on-shell pole that an algebraic asymptotic series omits. If cancellation or a
sign diagnostic rejects those identities, the code uses direct positive
quadrature.

For an unrestricted optical kernel, define the continued divided difference

$$
\mathcal D_m(a,b)=
\begin{cases}
[H_m(a)-H_m(b)]/(b-a),&a\ne b,\\
-H'_m(a),&a=b.
\end{cases}
$$

Then

$$
K_m(z_1,z_2)=
\frac{\mathrm{Re}\thinspace{}\mathcal D_m(z_1^*,z_2)
-\mathrm{Re}\thinspace{}\mathcal D_m(z_1,z_2)}{2\pi^2}.
$$

Near-coincident arguments use a derivative limit rather than subtracting two
nearby Hilbert transforms. Non-finite values, a forbidden sign, a Gaussian pole
outside the safe asymptotic sector, or an estimated cancellation loss trigger
direct epsilon quadrature. Finite-band odd kernels are symmetry-paired in this
fallback so a tiny residual is not lost by final subtraction. Gaussian fallback
quadrature separates the central region, narrow peaks, and tails.

Restricted and tabulated kernels use tangent coordinates for peaks inside the
energy interval and logarithmic distance coordinates for peaks outside it.
Well-separated optical peaks receive separate local regions; finite-band edge
charts, scaled integrands, and compensated sums protect very narrow or
sub-ULP geometries. A custom `m=0`, `n=0` kernel instead integrates its
interpolation spline once and caches the mass.

The remaining frequency integral uses adaptive GSL `qag` quadrature, with a
15-point Gauss-Kronrod rule and a workspace of 1,000 intervals by default.
Restricted outer integrals add linewidth-aware crossings of energy-window
boundaries and nearby extrema. Self-energy table boundaries, including shifted
optical boundaries, are also explicit quadrature partitions when legacy
extrapolation is enabled. A positive spectral-frequency window adds its
unshifted and optical-shifted mask boundaries.

The user options `-a`, `-r`, and `-k` govern the outer quadrature and applicable
custom-kernel pieces. Built-in fallback paths use tighter internal policies, and
restricted breakpoint integrations use fixed internal rules. Every adaptive
status and result is checked; route-specific roundoff or singularity failures
are accepted only after a stated error bound or independent fixed-rule check.
Unverified failures are fatal unless `-E` explicitly permits a finite partial
result. In verbose mode, `Error` is only the outer frequency-quadrature estimate;
nested epsilon uncertainty is checked separately and is not folded into it.

## Validation

The repository contains complementary DC, clean-kernel, restricted-window,
optical, Faddeeva, input, and DMFT benchmark checks. `make check` runs the
embedded Faddeeva suite, strict runner self-tests, all 384 legacy rows, all
pointwise kernel tables, restricted-window and input command tests, optical
checks, and the Bethe optical comparison without invoking Mathematica.

### Legacy Regression Data

`regression/` contains constant and linearly varying self-energies together
with historical reference values over a range of temperatures and kernel
indices. Run the complete manifest with:

```bash
cd regression
./run_tests --all
```

To run one table, pass its suite number and optional suffix, for example
`./run_tests 22 b`. Suite `11` uses `ReSigma=0`, `ImSigma=-0.01` and includes
the broader `n=0,1,2,3` set. The runner validates child status and output and
exits nonzero on any failed row. It accepts
$|\Delta|<10^{-4}$ when $|R|<10^{-6}$ and
$|\Delta|/|R|<10^{-4}$ otherwise.

### Clean-Kernel Regression

The pointwise reference reader accepts a computed value $V$ for reference $R$
when

$$
|V-R|\le\max(a,r|R|).
$$

The frozen tables are:

| Reference table | Rows | `mpmath` working precision | $r$ | $a$ |
|---|---:|---:|---:|---:|
| `clean_kernel_references.dat` | 370 | 160 digits | `5e-11` | `1e-300` |
| `restricted_kernel_references.dat` | 128 | 180 digits | `5e-9` | `1e-290` |
| `restricted_optical_references.dat` | 136 | 250 digits | `2e-8` | `1e-290` |
| `clean_optical_references.dat` | 58 | 250 digits | `5e-8` | `1e-300` |

Each generator first converts textual coordinates through binary64, matching
the C++ inputs, then performs high-precision quadrature and prints 25
significant digits. The generated files are frozen regression inputs; normal
testing and source refactoring do not regenerate them.

`regression/clean_kernel_references.dat` covers every built-in kernel with
`n=1,2,3`: interior points, exact and nearby band edges, clean exterior points,
Gaussian pole/tail transitions, broad linewidths, and exact spatial and
linewidth parity. Its generator uses the substitution

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
controlled clean-limit asymptotics where direct quadrature becomes singular.
They include exact odd symmetry, sub-ULP transformed spans, finite-band edges,
adjacent-double optical centers, huge Gaussian windows, and exterior linewidths
down to `1e-200`. Full-domain two-frequency values are generated separately by
`generate_clean_optical_references.py` and stored in
`clean_optical_references.dat`.
`run_window_tests` additionally checks command parsing, exact `M=0` behavior,
full and empty intersections, clean and narrow-domain tabulated kernels,
optical and occupied moments, exact `n=0` spline mass, `n=1` DOS output, and
narrow outer-window crossings.

The `regression/faddeeva_tests` target builds the bundled implementation with
its upstream `TEST_FADDEEVA` suite. It checks the complex error function,
`erf`, `erfc`, `erfi`, `erfcx`, and Dawson values against embedded
WolframAlpha/Maple references, including extreme and non-finite arguments, at a
maximum componentwise relative error of `1e-13`. Real/complex consistency is
also sampled over 10,000 logarithmically spaced magnitudes. The clean Gaussian
$J$ and $K$ tables then test how `bubble` uses these functions.

### Optical Regression

`regression/run_optical_tests` checks the exact `-O 0` DC path, small-frequency
convergence, frozen finite-frequency values for the built-in kernel families, a
large-frequency shifted-range case, a tabulated `m=0` kernel, shifted
self-energy boundaries, and invalid option combinations. It also retains
regressions for a far-exterior finite-band cancellation and an exponentially
small Gaussian exterior pole. It can be run directly with:

```bash
./regression/run_optical_tests
```

### Bethe-Lattice DMFT Benchmark

`Bethe_lattice_test/` contains a real-frequency DMFT self-energy and independent
Mathematica reference values for conductivity, resistivity, thermopower,
thermal conductivity, Lorenz ratio, and $ZT$. With the documented
normalization, the current executable agrees with these references to about
$10^{-5}$ relative accuracy. The occupied-energy check agrees to about
$1.3\times10^{-7}$ relative accuracy.

The finite-frequency reference contains 50 optical frequencies. Compare that
mesh with the frozen Mathematica output without changing files using:

```bash
./Bethe_lattice_test/cond.opt.geo --check
```

The script applies the benchmark's $2\pi$ normalization, `20T` cutoff, and
$10^{-16}$ self-energy clipping floor. It accepts
$|\Delta|\le3\times10^{-11}+3\times10^{-5}|\mathrm{reference}|$, reports the
largest errors, and exits unsuccessfully on any mismatch. The absolute allowance
covers one tiny optical tail where the frozen Mathematica interpolant clips
input knots but not noncausal inter-knot overshoot. `cond.opt.geo.mma` is retained
only as provenance; validation reads `cond.opt.geo-mma.dat` and never launches
Mathematica.

Running `./Bethe_lattice_test/cond.opt.geo` without `--check` performs the same
comparison and also rewrites `cond.opt.geo.dat`.

Some scripts in this directory belong to the original larger DMFT/NRG workflow
and require tools that are not distributed here. The `bubble` invocations and
the checked-in input/reference data are sufficient for the examples above.

## Precision

For precision-sensitive work, vary `-a`, `-r`, `-c`, `-k`, `-s`, `-M`, `-Mp`,
and the input grid, and confirm that the physical result is stable.

## Method Notes

The original derivation, kernel definitions, and implementation discussion are
available in:

- [`notes/bubble.pdf`](notes/bubble.pdf)
- [`notes/bubble.tex`](notes/bubble.tex)

The invocation section of the 2017 notes describes an earlier
interface. The command line documented in this README and printed by the
current executable is authoritative for version 1.11.

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

- **Amina Alić**: original implementation
- **Rok Žitko**: method notes and subsequent development

Issues and pull requests that improve numerical robustness, validation, or
support for additional lattice transport functions are welcome.

## License

The original project code and data are released under the
[GNU General Public License, version 3 or later](https://www.gnu.org/licenses/gpl-3.0.html).

`Faddeeva.cc` and `Faddeeva.hh` are by Steven G. Johnson, Massachusetts
Institute of Technology, and are distributed under the MIT license included in
their source headers.
