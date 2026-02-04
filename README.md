# DMD-Assigments
This is the DMD assignment for new student entering the Hall group.

Programs to download

Windows:

Notepad++: Text editor that will display code neatly

WinSCP: File transfer to lab computers

PuTTY: Remote access to lab computers

Spring First Year

Choose a programming language (Fortran is useful for working with PRIME20 in the future, C++ and Python are currently useful for peptide design projects).

Summer/Fall First Year

- Hard Sphere Code:
  - Create a code to run simulations of hard spheres. To do this you can use the two provided Fortran codes as a guide, and work to rewrite these codes into C++. Your code will need to
    - Place spheres in a simulation box on an FCC lattice (subroutine createcoord in swnvt4.f90)
    - Assign spheres an initial velocity based on a Gaussian distribution (subroutines createvel, gauss in swnvt4.f90)
    - Ensure that none of the particles are overlapping (subroutine check in hardspehrerd.f90)
    - Calculate initial collision times (subroutine uplist in hardsphererd.f90)
    - Begin moving the particles and letting them collide (Start of Dynamics section in hardsphererd.f90, requires subroutines bump, uplist, dnlist in hardsphererd.f90)
    - Compute compressibility factor (done using virial accumulator and the calculations after END OF DYNAMICS in hardsphererd.f90)
    - Compute radial distribution function (subroutine grsort in hardsphererd.f90)
  - compute the compressibility factor at various densities and compare to literature values
  - create a radial distribution diagram for systems at different densities and compare to literature

OTHER THINGS TO CONSIDER:

- hardsphererd.f90 takes reduced (number) density as an input, most of the literature reports values in terms of η (volume fraction), so you will need to convert between the two--also Dr. Hall likes to see volume fractions rather than reduced density.
- If you want to see what is happening in your simulations, you can use VMD to visualize how the spheres are positioned in the box. Once you get the particles placed on the initial lattice, you can use VMD to check and see that it is right. Ask for code on writing out files that VMD can read, as the program requires specific input files.
- Valgrind can be a useful program to finding errors in your code, please see the separate document for instructions on using this.

References/Resources

- Fortran hard sphere code (note that this code does not have information on creating the initial configuration and assigning initial velocity, use the square well code swnvt.f90 code for this)
- Allen, Tildesley. Computer Simulation of Liquids, 1987.
- Carnahan, N.F., Starling, K.E., Equation of State for Nonattracting Rigid Spheres.
- Alder, B.J., Hoover, W.G., Young, D.A., Studies in Molecular Dynamics. V.
- Alder, B.J., Wainwright, T.E., Studies in Molecular Dynamics. I.
- HardSphereMD_2 document (included in literature reference folder--explains the basics of pairwise collisions, the integration strategy, uplists, and downlists with figures)

- Square-well Sphere Code:
  - modify your hard sphere code to include an attractive interaction based on a square well potential (will require changes in bump, uplist, dnlist)
  - use an Andersen thermostat to maintain a constant temperature (subroutine ghostcoll in swnvt4.f90)
  - Write code to calculate kinetic and potential energy of the system, use kinetic energy to calculate the system temperature and ensure your thermostat is keeping the temperature constant (subroutine kecalc, pecalc in swnvt4.f90)
  - Calculate the total energy to ensure the system has reached equilibrium
  - compute the compressibility factor at various temperatures and densities and compare to literature values
  - create a radial distribution diagram for systems at different densities

References/Resources

- Fortran square-well sphere code
- Allen, Tildesley. Computer Simulation of Liquids, 1987.
- Henderson D. Madden W.G., Fitts D.D. Monte Carlo and hypernetted chain equation of state for the square-well fluid. 1976.
- Lee, R.J., Chao, K.C. Equation of state for square-well fluids. 1988.
- Tang, Y., Lu, B.C.-Y. An analytical analysis of the square-well fluid behaviors. 1994

- Square-well Chains:
  - modify your square-well code to model chains of either 4 or 16 spheres
    - spheres next to each other on a chain are allowed to move between σ(1±δ)
    - There are two ways to create the initial configuration: random configuration (done in the Fortran Code) or along a "lattice". Try the random configuration first and if it doesn't work out, move on to trying a lattice.
  - compute the compressibility factor at various temperatures and densities and compare to literature values
  - create a radial distribution diagram for systems at different densities

References/Resources

- Fortran code swcnopdb.f90
- Yeom, M.S., Chang, J., Kim, H. Development of the semi-empirical equation of state for square-well chain fluid based on statistical associating fluid theory (SAFT). 1999
- Tavares, F., Chang, J., Sandler, S. Equation of state for the square-well chain fluid based on the dimer version of Wertheim's perturbation theory. 1995

- Efficiency Techniques
  - Implement the efficiency techniques discussed in: Smith, Hall, Freeman. Molecular Dynamics for Polymeric Fluids Using Discontinuous Potentials. Journal of Computational Physics, 1996, **134**, 16-30
  - Start with False Positioning (the easiest to implement, but without neighbor lists will only slow down your code)
  - Move on to Neighbor lists/Link Lists
    - These work together so it is a good idea to work on them at the same time.
