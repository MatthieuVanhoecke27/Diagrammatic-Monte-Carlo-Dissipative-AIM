# Real-Time Diagrammatic Monte Carlo Solver (DiagMC)

![Language](https://img.shields.io/badge/C++-17-blue.svg) ![Library](https://img.shields.io/badge/Armadillo-Linear%20Algebra-red.svg) ![Status](https://img.shields.io/badge/Status-Research%20Grade-success.svg)

## Project Overview

This repository contains a high-performance C++ implementation of the **Continuous-Time Diagrammatic Monte Carlo (CT-QMC)** algorithm. It is designed to solve **Open Quantum Systems** and **Non-Hermitian Hamiltonians** in real-time.

Unlike standard Exact Diagonalization techniques which scale exponentially with system size, this stochastic solver samples the perturbation series of the impurity model, allowing for the simulation of complex quantum impurity problems with high accuracy.

**Key Application Areas:** Quantum Transport, Kondo Physics, Non-Equilibrium Dynamics.

---

## Technical Highlights

This project demonstrates proficiency in **High-Performance Computing (HPC)** and **Algorithmic Optimization**:

* **Fast Rank-1 Matrix Updates:** Implemented the **Sherman-Morrison formula** to update inverse matrices during Monte Carlo steps. This reduces the complexity of determining acceptance probabilities from $O(N^3)$ (standard inversion) to **$O(N^2)$**, significantly speeding up the sampling loop.
* **Custom Data Structures:** Built a specialized **Doubly Linked List** (`clist`) to manage the continuous-time configuration of operators. This ensures $O(1)$ insertion/deletion time for diagrammatic updates, avoiding the memory overhead of `std::vector` resizing.
* **Memory Management:** Rigorous management of dynamic memory allocation for Worldline nodes to prevent memory leaks during millions of stochastic steps.
* **Numerical Stability:** Handled complex number arithmetic and fermionic sign problems inherent in Quantum Monte Carlo simulations.

---

## The Algorithm: Metropolis-Hastings

The simulation explores the configuration space of the perturbation expansion using the Metropolis-Hastings algorithm. The core engine relies on three ergodic updates:

1.  **Insertion :** Adds a pair of creation/annihilation operators at random times.
2.  **Removal :** Removes an existing pair of operators.
3.  **Shifting :** Moves an operator in time, updating the hybridization matrix continuously via row/column shifts.

### Code Structure

| File | Description |
| :--- | :--- |
| `main.cpp` | Entry point. Parses CLI args, initializes the bath, and drives the time-evolution loop. |
| `MonteCarlo_Adding.cpp` | Implements the **add** move. Handles matrix expansion and determinant ratios. |
| `MonteCarlo_Removing.cpp` | Implements the **remove** move. Handles matrix contraction. |
| `MonteCarlo_Shifting.cpp` | Implements the **Shift** move. Handles row/column updates in the hybridization matrix. |
| `clist.cpp` | Custom Linked-List implementation for managing the time-ordered operator sequence. |

---

## Dependencies & Build Instructions

### Prerequisites
* **C++ Compiler** (GCC 7+ or Clang) supporting C++17.
* **Armadillo**: C++ linear algebra library for high-speed matrix operations.
* **BLAS/LAPACK**: Underlying linear algebra kernels.

### Compilation

    ```bash
    g++ main.cpp MonteCarlo_Adding.cpp MonteCarlo_Removing.cpp MonteCarlo_Shifting.cpp clist.cpp -o DiagMC -O3 -larmadillo -std=c++17
    ```
*(Note: The `-O3` flag is critical for enabling compiler optimizations for Monte Carlo loops.)*

---

## Usage

The executable takes physical parameters and simulation constraints as command-line arguments:

```bash
./DiagMC [SavePath] [MC_Steps] [EffHamil] [T_Max_Steps] [U] [GammaUp] [GammaDown] [FinalTime]
```

### Parameters
* **MC_Steps**: Number of Metropolis steps (e.g., $10^6$ to $10^8$).
* **U**: Coulomb interaction strength.
* **GammaUp/Down**: Dephasing rates.
* **FinalTime**: The maximum physical time to simulate evolution.

---

## Sample Output

The simulation outputs time-resolved observables, including:
1. **Density Matrix components:** $\rho(t)$
2. **Number Operator:** $\langle n(t) \rangle$
3. **Average Perturbation Order:** To monitor convergence.

---
## Citation

If you use this solver in your research or project, please cite it as follows:

**BibTeX:**
```bibtex
@article{Vanhoecke_2024,
   title={Diagrammatic Monte Carlo for dissipative quantum impurity models},
   volume={109},
   ISSN={2469-9969},
   url={http://dx.doi.org/10.1103/PhysRevB.109.125125},
   DOI={10.1103/physrevb.109.125125},
   number={12},
   journal={Physical Review B},
   publisher={American Physical Society (APS)},
   author={Vanhoecke, Matthieu and Schirò, Marco},
   year={2024},
   month=mar }

```
---
**Author:** Matthieu Vanhoecke