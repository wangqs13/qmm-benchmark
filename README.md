# qmm-benchmark

This is a benchmark for equivalence checking of quantum circuits using quantum Mealy machine [WLY19].

## Basic Definitions

A quantum Mealy machine (QMM) is a 5-tuple <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{M}&space;=&space;(\Sigma,&space;\Gamma,&space;\mathcal{H},&space;U,&space;M)" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{M}&space;=&space;(\Sigma,&space;\Gamma,&space;\mathcal{H},&space;U,&space;M)" title="\mathcal{M} = (\Sigma, \Gamma, \mathcal{H}, U, M)" /></a>, where:
* <a href="https://www.codecogs.com/eqnedit.php?latex=\Sigma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Sigma" title="\Sigma" /></a> is a finite input alphabet;
* <a href="https://www.codecogs.com/eqnedit.php?latex=\Gamma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\Gamma" title="\Gamma" /></a> is a finite output alphabet;
* <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{H}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{H}" title="\mathcal{H}" /></a> is a finite-dimensional Hilbert space;
* <a href="https://www.codecogs.com/eqnedit.php?latex=U&space;=&space;\{&space;U_\sigma:&space;\sigma&space;\in&space;\Sigma&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U&space;=&space;\{&space;U_\sigma:&space;\sigma&space;\in&space;\Sigma&space;\}" title="U = \{ U_\sigma: \sigma \in \Sigma \}" /></a> is a set of unitary operators. For each <a href="https://www.codecogs.com/eqnedit.php?latex=\sigma&space;\in&space;\Sigma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sigma&space;\in&space;\Sigma" title="\sigma \in \Sigma" /></a>, <a href="https://www.codecogs.com/eqnedit.php?latex=U_\sigma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?U_\sigma" title="U_\sigma" /></a> is a unitary operator on <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{H}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{H}" title="\mathcal{H}" /></a>; and
* <a href="https://www.codecogs.com/eqnedit.php?latex=M&space;=&space;\{&space;M_m:&space;m&space;\in&space;\Gamma&space;\}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?M&space;=&space;\{&space;M_m:&space;m&space;\in&space;\Gamma&space;\}" title="M = \{ M_m: m \in \Gamma \}" /></a> is a quantum measurement in <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{H}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{H}" title="\mathcal{H}" /></a>, that is, <a href="https://www.codecogs.com/eqnedit.php?latex=M_m" target="_blank"><img src="https://latex.codecogs.com/gif.latex?M_m" title="M_m" /></a> is a linear operator on <a href="https://www.codecogs.com/eqnedit.php?latex=\mathcal{H}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathcal{H}" title="\mathcal{H}" /></a> for each <a href="https://www.codecogs.com/eqnedit.php?latex=m&space;\in&space;\Gamma" target="_blank"><img src="https://latex.codecogs.com/gif.latex?m&space;\in&space;\Gamma" title="m \in \Gamma" /></a> and <a href="https://www.codecogs.com/eqnedit.php?latex=\sum_m&space;M_m^\dag&space;M_m&space;=&space;I" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\sum_m&space;M_m^\dag&space;M_m&space;=&space;I" title="\sum_m M_m^\dag M_m = I" /></a>.

See [WLY19] for more details.

## Data Format

testXXX.data is a binary file, whose contents are listed below:

```
<dim, int>
<sigma, int>
<gamma, int>
<limit, int>
<eps, double>
for a in 0..sigma-1:
    for i in 0..dim-1:
        for j in 0..dim-1:
            <unitary[a][i][j].real, double>
            <unitary[a][i][j].imag, double>
for b in 0..sigma-1:
    for i in 0..dim-1:
        for j in 0..dim-1:
            <measure[b][i][j].real, double>
            <measure[b][i][j].imag, double>
for i in 0..dim-1:
    for j in 0..dim-1:
        <rho1[i][j].real, double>
        <rho1[i][j].imag, double>
for i in 0..dim-1:
    for j in 0..dim-1:
        <rho2[i][j].real, double>
        <rho2[i][j].imag, double>
```

The offsets of the first few data are shown below for calibration:

```
offset        <content, type>
 0000	        <dim, int>
 0004	        <sigma, int>
 0008	        <gamma, int>
 000c	        <limit, int>
 0010	        <eps, double>
 ....
```

## Examples

We select test002 as an illustrative example, in which case:
* dim = 4.
* sigma = 2.
* gamma = 2.
* limit = -1.
* eps = 1e-8.
* <a href="https://www.codecogs.com/eqnedit.php?latex=\text{unitary[0]}&space;=&space;\frac&space;1&space;{\sqrt2}&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;1&space;\\&space;1&space;&&space;0&space;&&space;-1&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;-1&space;\end{bmatrix},\text{unitary[1]}&space;=&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\end{bmatrix}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\text{unitary[0]}&space;=&space;\frac&space;1&space;{\sqrt2}&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;1&space;\\&space;1&space;&&space;0&space;&&space;-1&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;-1&space;\end{bmatrix},\text{unitary[1]}&space;=&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\end{bmatrix}." title="\text{unitary[0]} = \frac 1 {\sqrt2} \begin{bmatrix} 1 & 0 & 1 & 0 \\ 0 & 1 & 0 & 1 \\ 1 & 0 & -1 & 0 \\ 0 & 1 & 0 & -1 \end{bmatrix},\text{unitary[1]} = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 1 \\ 0 & 0 & 1 & 0 \end{bmatrix}." /></a>
* <a href="https://www.codecogs.com/eqnedit.php?latex=\text{measure[0]}&space;=&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix},&space;\text{measure[1]}&space;=&space;\begin{bmatrix}&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\end{bmatrix}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\text{measure[0]}&space;=&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix},&space;\text{measure[1]}&space;=&space;\begin{bmatrix}&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;1&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;1&space;\end{bmatrix}." title="\text{measure[0]} = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix}, \text{measure[1]} = \begin{bmatrix} 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 1 & 0 \\ 0 & 0 & 0 & 1 \end{bmatrix}." /></a>
* <a href="https://www.codecogs.com/eqnedit.php?latex=\text{rho1}&space;=&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix},&space;\text{rho2}&space;=&space;\begin{bmatrix}&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\text{rho1}&space;=&space;\begin{bmatrix}&space;1&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix},&space;\text{rho2}&space;=&space;\begin{bmatrix}&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;1&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\\&space;0&space;&&space;0&space;&&space;0&space;&&space;0&space;\end{bmatrix}." title="\text{rho1} = \begin{bmatrix} 1 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix}, \text{rho2} = \begin{bmatrix} 0 & 0 & 0 & 0 \\ 0 & 1 & 0 & 0 \\ 0 & 0 & 0 & 0 \\ 0 & 0 & 0 & 0 \end{bmatrix}." /></a>

## References

[WLY19] Q. S. Wang, J. Y. Liu and M. S. Ying. Equivalence Checking of Quantum Finite-State
Machines. In: [arXiv:1901.02173](https://arxiv.org/pdf/1901.02173.pdf).
