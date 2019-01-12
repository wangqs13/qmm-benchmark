# qmm-benchmark

This is a benchmark for equivalence checking of quantum circuits using quantum Mealy machine [WLY19].

## Basic Definitions

A quantum Mealy machine (QMM) is a 5-tuple $\mathcal{M} = (\Sigma, \Gamma, \mathcal{H}, U, M)$, where:
* $\Sigma$ is a finite input alphabet;
* $\Gamma$ is a finite output alphabet;
* $\mathcal{H}$ is a finite-dimensional Hilbert space;
* $U = \{ U_\sigma: \sigma \in \Sigma \}$ is a set of unitary operators. For each $\sigma \in \Sigma$, $U_\sigma$ is a unitary operator on $\mathcal{H}$; and
* $M = \{ M_m: m \in \Gamma \}$ is a quantum measurement in $\mathcal{H}$, that is, $M_m$ is a linear operator on $\mathcal{H}$ for each $m \in \Gamma$ and $\sum_m M_m^\dag M_m = I$.

See [WLY19] for more details.

## Data Format

testXXX.data is a binary file, whose contents are listed below:
> <dim, int>
> <sigma, int>
> <gamma, int>
> <limit, int>
> <eps, double>
> for a in 0..sigma-1:
>     for i in 0..dim-1:
>         for j in 0..dim-1:
>             <unitary[a][i][j].real, double>
>             <unitary[a][i][j].imag, double>
> for b in 0..sigma-1:
>     for i in 0..dim-1:
>         for j in 0..dim-1:
>             <measure[b][i][j].real, double>
>             <measure[b][i][j].imag, double>
> for i in 0..dim-1:
>     for j in 0..dim-1:
>         <rho1[i][j].real, double>
>         <rho1[i][j].imag, double>
> for i in 0..dim-1:
>     for j in 0..dim-1:
>         <rho2[i][j].real, double>
>         <rho2[i][j].imag, double>

The offsets of the first few data are shown below for calibration:
> offset        <content, type>
>  0000	        <dim, int>
>  0004	        <sigma, int>
>  0008	        <gamma, int>
>  000c	        <limit, int>
>  0010	        <eps, double>
>  ....

## Examples

We select test002 as an illustrative example, in which case:
* dim = 4.
* sigma = 2.
* gamma = 2.
* limit = -1.
* eps = 1e-8.
* $\text{unitary[0]} = \frac 1 {\sqrt2} \begin{bmatrix} 
  1 & 0 & 1 & 0 \\
  0 & 1 & 0 & 1 \\
  1 & 0 & -1 & 0 \\
  0 & 1 & 0 & -1
\end{bmatrix}$,
$\text{unitary[1]} = \begin{bmatrix} 
  1 & 0 & 0 & 0 \\
  0 & 1 & 0 & 0 \\
  0 & 0 & 0 & 1 \\
  0 & 0 & 1 & 0
\end{bmatrix}$.
* $\text{measure[0]} = \begin{bmatrix} 
  1 & 0 & 0 & 0 \\
  0 & 1 & 0 & 0 \\
  0 & 0 & 0 & 0 \\
  0 & 0 & 0 & 0
\end{bmatrix}$,
$\text{measure[1]} = \begin{bmatrix} 
  0 & 0 & 0 & 0 \\
  0 & 0 & 0 & 0 \\
  0 & 0 & 1 & 0 \\
  0 & 0 & 0 & 1
\end{bmatrix}$.
* $\text{rho1} = \begin{bmatrix} 
  1 & 0 & 0 & 0 \\
  0 & 0 & 0 & 0 \\
  0 & 0 & 0 & 0 \\
  0 & 0 & 0 & 0
\end{bmatrix}$,
$\text{rho2} = \begin{bmatrix} 
  0 & 0 & 0 & 0 \\
  0 & 1 & 0 & 0 \\
  0 & 0 & 0 & 0 \\
  0 & 0 & 0 & 0
\end{bmatrix}$.

## References

[WLY19] Q. S. Wang, J. Y. Liu and M. S. Ying. Equivalence Checking of Quantum Finite-State
Machines. In: [arXiv:1901.02173](https://arxiv.org/pdf/1901.02173.pdf).
