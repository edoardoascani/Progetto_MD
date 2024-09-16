! Author: Edoardo Ascani
!
! Description: 
! This is a molecular dynamics program developed to simulate atom particles in a gas phase. Verlet integration
! was employed to numerically integrate Newton’s equation of motion, allowing for
! the calculation of particle trajectories. The interaction between gas particles was simulated
! through a Lennard-Jones potential, where every atom experiences a cumulative Lennard-Jones force generated
! by other particles.
! 
!
! Run:
! 1. Create an input file named "input.TXT" with a structure similar to the provided samples 
! 2. Compile and run the program using the following commands:
!    gfortran -c utilities.f90
!    gfortran utilities.o MD.f90 -o ./simulation
! 3. Two outputs will be generated:
!    - "md.xyz": a file containing the coordinates of the system at every step (trajectory). 
!    - "energy.dat": a file containing the energy data of the system. 


program noblegas_simulation

  use functions
  
  implicit none

! Parameters for conversions to atomic units

  double precision, parameter  :: kcal_to_Hartree = 0.0015936015d0
  double precision, parameter  :: ev_to_Hartree = 0.0367492929d0
  double precision, parameter  :: ang_to_bohr = 1.8897261d0
  double precision, parameter  :: fs_to_timeau = 41.34137314d0
  double precision, parameter  :: amu_to_au = 1822.8884850d0
  double precision, parameter  :: timeau_to_fs = 0.02418884254d0

! Variable declarations

  integer                                        :: i, j, k, l, step, nstep, natom
  double precision                               :: timestep, mass, eps, sig, kin_ener, pot_ener, energy, time
  double precision, allocatable, dimension(:,:)  :: positions, velocities, distances, forces
  character(len=100)                             :: filename = "input.TXT"
  ! Data taken from article "Neutral Gases Adsorption With Hydrogen on Silicon Nanotubes: A Fuel Cell Investigation"
  character(len=2), allocatable, dimension(:)    :: atom


! Opening output files

  open(10, file="trajectory.xyz")
  open(11, file="energy.dat")

! Reading input variables by calling the subroutine from the module

  call lecture(filename, positions, atom, velocities, nstep, mass, timestep, sig, eps, natom)

! Allocate memory

  allocate(forces(natom,3), distances(natom,natom))

! Write input parameters
  
  write(6,"(A)") "Initial parameters"
  write(6,*)
  write(6,"(A,I4)") "number of steps : ", nstep
  write(6,"(A,F6.3)") "time step (fs)  : ", timestep
  write(6,"(A,I2)") "number of atoms : ", natom
  write(6,"(A,F9.6)") "masses (amu)    : ", mass
  write(6,"(A,F9.6)") "epsilon (eV)   : ", eps 
  write(6,"(A,F9.6)") "sigma (ang)     : ", sig
  write(6,*)

! Write geometry information

  write(6,"(A)") "Geometry (ang)"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), positions(i,:)
  enddo
  write(6,*)

! Write velocities information

  write(6,"(A)") "Velocities (ang/fs)"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), velocities(i,:)
  enddo

! Conversion of quantities to atomic units

  mass = mass * amu_to_au
  timestep = timestep * fs_to_timeau
  sig = sig * ang_to_bohr
  eps = eps * ev_to_Hartree
  positions = positions * ang_to_bohr
  velocities = velocities * ang_to_bohr / fs_to_timeau

! Write parameters at step 0 in atomic units

  write(6,*)
  write(6,"(A)") "-------------------------------------"
  write(6,"(A)") "Initial parameters in atomic units"
  write(6,*) " "
  write(6,"(A,F10.6)") "time step  : ", timestep
  write(6,"(A,F12.6)") "masses     : ", mass
  write(6,"(A,E13.6)") "epsilon    : ", eps 
  write(6,"(A,F9.6)") "sigma      : ", sig
  write(6,*)

! Write geometry information

  write(6,"(A)") "Geometry"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), positions(i,:)
  enddo
  write(6,*)

! Write velocities 

  write(6,"(A)") "Velocities"
  write(6,"(8X,3(A1,10X))") "x", "y", "z"
  do i=1,natom
    write(6,"(A2,1X,3(F10.6,1X))") atom(i), velocities(i,:)
  enddo

! Compute distances for LJ potential

  call get_distances(positions, distances, natom)

! Compute LJ forces (need for Verlet)

  call get_forces(forces, distances, positions, sig, eps, natom)

! Calculate total energy

  call get_kinetic(velocities, natom, kin_ener, mass)
  call get_potential(distances, natom, sig, eps, pot_ener)
  energy = kin_ener + pot_ener
  
! Initialization of step

  step = 0

! Write geometry and energy to the output

  time = step * timestep * timeau_to_fs
  write(10,"(I3)") natom
  write(10,"(A,I4,A,F6.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, & 
                                     & ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
  do i = 1,natom
     write(10,*) atom(i), positions(i,:)
  enddo
  write(11,*) "# step, time, energy, kin_ener, pot_ener"
  write(11,*) step, time, energy, kin_ener, pot_ener


! *****************************************************

! Trajectories are computed in a loop that repeats for the number of steps we set

  do step = 1, nstep

! Calculate velocities using Verlet integration

  do i=1,natom
    do j=1,3
      velocities(i,j)=velocities(i,j)+0.5d0*timestep*forces(i,j)/mass
    enddo
  enddo

  do i=1,natom
    do j=1,3
      positions(i,j)=positions(i,j)+timestep*velocities(i,j)
    enddo
  enddo

  call get_distances(positions,distances,natom)
  call get_forces(forces,distances,positions,sig,eps,natom)

  do i=1,natom
    do j=1,3
      velocities(i,j)=velocities(i,j)+0.5d0*timestep*forces(i,j)/mass
    enddo
  enddo

! Calculate energies by calling the subroutines get_kinetic and get_potential from the module

  call get_kinetic(velocities, natom, kin_ener, mass)
  call get_potential(distances, natom, sig, eps, pot_ener)
  energy = kin_ener + pot_ener

! Debug
! write(6,*) step, ") ", ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener

! Write geometry and energy to the output

  time = step * timestep * timeau_to_fs
  write(10,"(I3)") natom                                 
  write(10,"(A,I4,A,F9.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, & 
                                     & ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
  do i = 1, natom
     write(10,*) atom(i), positions(i,:)
  enddo

  write(11,*) step, time, energy, kin_ener, pot_ener
  
  enddo

  close(10)
  close(11)

  deallocate(forces, distances)

end program noblegas_simulation
