# DMD Introductory Exercises  
**Hall Group – New Graduate Student Onboarding**

This document outlines the recommended training path for new graduate students who will be working with **Discontinuous Molecular Dynamics (DMD)** and related simulation codes in the group. The goal is to build a strong foundation in Linux, programming, and statistical mechanics before you start your research projects.

---

## Recommended Software

The following tools are **suggested** to make it easier to work with lab computing clusters and the NCSU HPC Hazel.  
Students may use equivalent tools based on personal preference or operating system.

### Commonly Used Options (Windows)
- **Notepad++** – Lightweight text editor for viewing and editing source code  

- **WinSCP** – Secure file transfer between local machines and lab computers  

- **PuTTY** – SSH client for remote access to lab computers  
  *(Alternatives: Windows Terminal, OpenSSH)*

### macOS / Linux
- Mac/Linux users can use native terminal and SCP tools

> **Note:**  
> Regardless of tool choice, all students should be comfortable with:
> - Editing plain-text source files  
> - Transferring files to and from remote machines  
> - Accessing HPC and lab clusters via SSH  

---

## Spring Semester – First Year  
**Foundations**

During your first semester, focus on building computational fundamentals.

### Core Skills
- Become comfortable working in a **Linux environment**
- Choose and begin learning a **programming language** if you are not familiar with programing

> **Language guidance:**  
> - **Fortran** – highly recommended for future work with PRIME20  
> - **C/C++** – useful for performance-critical DMD development  
> - **Python** – useful for data analysis, automation, and visualization  

### Recommended Beginner Resources
- **Linux**:  
  https://www.coursera.org/learn/hands-on-introduction-to-linux-commands-and-shell-scripting  
- **Git/GitHub**:  
  https://www.coursera.org/learn/introduction-git-github  
- **C Programming**:  
  https://www.coursera.org/specializations/c-programming  
- **Fortran**:  
  https://fortran-lang.org/learn/  
- **Python**:  
  https://www.coursera.org/learn/python-programming-intro  

---

## Summer / Fall – First Year  
**DMD Programming Exercises**

The following exercises introduce core DMD concepts using progressively more complex systems. Students are encouraged to write clean, well-documented code and validate results against published literature.

---

## 1. Hard Sphere (HS) Simulations

### Objectives
Write a code to simulate hard spheres using DMD. Two existing Fortran codes are provided as references. You may rewrite the code in **Fortran, C/C++, or another suitable language**.

### Required Components
Your code should be able to:

1. **Initialize particle positions and velocities**
   - Place spheres in a simulation box on an **FCC lattice** (`coordinate_and_vel`)
   - Optionally generate **random initial positions** (recommended for very low or very high densities)
   - Assign initial velocities from a **Gaussian distribution**

2. **Ensure valid configurations**
   - Verify that particles do not overlap (`check_position`)

3. **Set up event scheduling**
   - Calculate initial collision times (`uplist`)
   - Advance particles and handle collisions (Start of *Dynamics* section in `HS_Haoyu.F90`)
   - Required subroutines:
     - `find_tmin_and_update`
     - `uplist`
     - `dnlist`

4. **Compute physical observables**
   - **Compressibility factor** using the virial accumulator  
     (lines 443–445 in `find_tmin_and_update`; see Allen & Tildesley, pp. 46–48)
   - **Radial distribution function (RDF)**  
     (`grsort` in `hardsphererd.f90`; Allen & Tildesley, pp. 54–55)

### Analysis Tasks
- Compute the compressibility factor over a range of densities and compare with literature values  
- Generate RDF plots at different densities and compare with published results  

### Notes and Practical Tips
- `HS_Haoyu.F90` uses **reduced number density**, while most literature reports **volume fraction (η)**  
  → Convert between the two  
  → Use **volume fraction** in plots and reports (Dr. Hall’s preference)
- Visualization:
  - Use **VMD** to inspect particle configurations
  - XYZ output can be enabled by uncommenting lines **95, 96, and 109** in `HS_Haoyu.F90`
  - Ask for help if you need example output formats

**References/Resources:**  
`Literature References/Hard Spheres`

---

## 2. Square-Well (SW) Spheres

### Objectives
Extend the hard-sphere code to include attractive interactions.

### Required Modifications
- Add a **square-well potential**
  - Modify `bump`, `uplist`, and `dnlist`
- Implement an **Andersen thermostat** (`ghostcoll` in `SW_Haoyu.f90`)
- Compute:
  - **Kinetic energy** and system temperature (`kecal`)
  - **Potential energy** (`pecal`)
  - **Total energy** to confirm equilibration

### Analysis Tasks
- Compute compressibility factors at multiple **temperatures and densities**
- Generate RDFs for different state points
- Compare all results with literature values

**References/Resources:**  
`Literature References/Square-well Spheres`

---

## 3. Square-Well Chains

### Objectives
Extend the square-well model to simulate **polymeric chains**.

### System Requirements
- Chains of **4 or 16 spheres**
- Bonded neighbors must satisfy:  
  `σ(1 − δ) ≤ r ≤ σ(1 + δ)`

### Initial Configurations
- Random initial configurations (recommended starting point)
- Lattice-based configurations if random initialization fails

### Analysis Tasks
- Compute compressibility factors vs. temperature and density
- Generate RDFs and compare with literature

**References/Resources:**  
`Literature References/Square-well Chains`

---

## 4. Efficiency Techniques (Advanced)

### Objectives
Improve performance using standard DMD efficiency techniques.

### Required Reading
Smith, Hall, Freeman  
*Molecular Dynamics for Polymeric Fluids Using Discontinuous Potentials*  
**Journal of Computational Physics**, 1996, 134, 16–30

### Implementation Path
1. **False positioning**
   - Easiest to implement
   - Limited benefit without neighbor lists
2. **Neighbor lists / link lists**
   - Implement together for best performance gains

### Additional Resources
See the directory **Former Student Notes**, which includes:
- Complete DMD assignments in C++ (Dr. Ryan Malony)
- Results and notes from Corey Febo and Haoyu Wang
- Van Nguyen’s personal notes on learning DMD
