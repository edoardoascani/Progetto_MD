! Author: Edoardo Ascani
!
! Brief description:
! This program is written in Fortran 90 and simulates molecular dynamics of atoms in a gas phase
! inside a box.
! It uses:w
! Verlet integration to numerically solve Newton’s equations of motion, calculating the
! trajectories of particles. The interaction between gas particles is modeled using the
! Lennard-Jones potential, where each atom experiences the cumulative force generated by all
! other atoms.
!
! How to Run:
! 1. Create an input file named "input.TXT" with the structure previously provided.
! 2. Compile the module by running: 
!    gfortran -c functions.F90
!    This will generate two files: functions.o and functions.mod.
! 3. Generate the executable by running:
!    gfortran functions.o nbsimulation.F90 -o ./simulation
!    This creates "simulation.exe" in your working directory.
! 4. Run the simulation with:
!    ./simulation
!
! Outputs:
! - "trajectory.xyz": Contains the system's coordinates at each simulation step.
! - "energy.dat": Contains the system's energy data at each simulation step.
! - "statistics.dat": Contains the system's statistic data at each simulation step.



program noblegas_simulation
!!
!!  
use functions
!  
implicit none

! Parameters for conversions to atomic units

  double precision, parameter  :: kcal_to_Hartree = 0.0015936015d0
  double precision, parameter  :: ev_to_Hartree = 0.0367492929d0
  double precision, parameter  :: ang_to_bohr = 1.8897261d0
  double precision, parameter  :: fs_to_timeau = 41.34137314d0
  double precision, parameter  :: amu_to_au = 1822.8884850d0
  double precision, parameter  :: timeau_to_fs = 0.02418884254d0

! Variable declarations

  integer                                        :: i, j, step, nstep, natom
  double precision                               :: timestep, mass, eps, sig, kin_ener, pot_ener, energy,prev_energy
  double precision                               :: pressure, temperature, virial, volume, temp, time, density
  double precision                               :: box_size, tol
  character(len=2), allocatable, dimension(:)    :: atom
  double precision, allocatable, dimension(:,:)  :: positions, velocities, distances, forces
  character(len=100)                             :: filename = "input.TXT"

  ! Data taken from article "Neutral Gases Adsorption With Hydrogen on Silicon Nanotubes: A Fuel Cell Investigation"

! Opening output files

  open(10, file="trajectory.xyz")
  open(11, file="energy.dat")
  open(12, file="statistics.dat")

! Reading input variables by calling the subroutine from the module

  call lecture(filename, positions, atom, nstep, box_size, temp, mass, timestep, sig, eps, natom)

! Allocate memory

  allocate(forces(natom,3), distances(natom,natom))
  allocate(velocities(natom,3))

! Write input parameters
  
  write(6,"(A)") "Initial parameters"
  write(6,*)
  write(6,"(A,I4)") "number of steps : ", nstep
  write(6,"(A,F6.3)") "time step (fs)  : ", timestep
  write(6,"(A,I2)") "number of atoms : ", natom
  write(6,"(A,F9.6)") "masses (amu)    : ", mass
  write(6,"(A,F9.6)") "epsilon (eV)   : ", eps 
  write(6,"(A,F9.6)") "sigma (ang)     : ", sig
  write(6,"(A,F10.6)") "temperature (K)     : ", temp
  write(6,"(A,F10.6)") "lato box     : ", box_size
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


! Write new velocities information

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
  box_size = box_size* ang_to_bohr

! Write parameters at step 0 in atomic units

  write(6,*)
  write(6,"(A)") "-------------------------------------"
  write(6,"(A)") "Initial parameters in atomic units"
  write(6,*) " "
  write(6,"(A,F10.6)") "time step  : ", timestep
  write(6,"(A,F12.6)") "masses     : ", mass
  write(6,"(A,E13.6)") "epsilon    : ", eps 
  write(6,"(A,F9.6)") "sigma      : ", sig
  write(6,"(A,F10.6)") "temperature (K)     : ", temp
  write(6,"(A,F10.6)") "lato box     : ", box_size
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

  call initialize_velocities(velocities, natom, temp, mass)

  volume = box_size **3


  density = natom* mass / volume

  tol = 1.0d-7

  call get_distances(positions, distances,box_size, natom)
  call get_forces(forces, virial, distances, positions, sig, eps, natom)

  call get_kinetic(velocities, natom, kin_ener, mass)
  call get_potential(distances, natom, sig, eps, pot_ener)

  energy = kin_ener + pot_ener
  
  call get_statistic(kin_ener, natom, volume, density, virial, temperature, pressure)
  
  step = 0

  time = step * timestep * timeau_to_fs

  write(10,"(I3)") natom
  write(10,"(A,I4,A,F6.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, & 
                                     & ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
  do i = 1,natom
     write(10,*) atom(i), positions(i,:)
  enddo

  write(11,*) "# step, time, energy, kin_ener, pot_ener"
  write(11,*) step, time, energy, kin_ener, pot_ener

  write(12,*) "# step, time, energy, temperature, pressure, virial, density"
  write(12,*) step, time, energy, temp, pressure, virial, density

  do step = 1, nstep

    do i=1,natom
      do j=1,3

        velocities(i,j)=velocities(i,j)+0.5d0*timestep*forces(i,j)/mass

      enddo
    enddo

    do i=1,natom
      do j=1,3

        positions(i,j) = positions(i,j) + timestep * velocities(i,j)
        
      enddo
    enddo

    do i=1,natom
      do j=1,3

        if (positions(i,j) > box_size) then

          positions(i,j) = positions(i,j) - box_size 

        else if (positions(i,j) < 0.0d0) then

          positions(i,j) = positions(i,j) + box_size 
          
        end if
      enddo
    enddo
    call get_distances(positions,distances,box_size,natom)
    call get_forces(forces,virial,distances,positions,sig,eps,natom)
    do i=1,natom
      do j=1,3

        velocities(i,j)=velocities(i,j)+0.5d0*timestep*forces(i,j)/mass

      enddo
    enddo

    call get_kinetic(velocities, natom, kin_ener, mass)
    call get_potential(distances, natom, sig, eps, pot_ener)

    energy = kin_ener + pot_ener

    if (step > 1 .and. abs(energy - prev_energy) >= tol) then

      print *, "Energy difference exceeds tolerance. Stopping simulation at step:", step

      exit  

  endif

    prev_energy = energy

    call get_statistic(kin_ener, natom, volume, density, virial, temperature, pressure)

    write(6,*) step, ") ", ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener 


    time = step * timestep * timeau_to_fs
    write(10,"(I3)") natom                                 
    write(10,"(A,I4,A,F9.3,3(A,E15.6))") "# step: ", step, ",  t(fs): ", time, & 
                                       &  ",   E =", energy, ",   Ek =", kin_ener, ",   Ep =", pot_ener
    do i = 1, natom
       write(10,*) atom(i), positions(i,:)
    enddo

    write(11,*) step, time, energy, kin_ener, pot_ener
    write(12,*) step, time, energy, temperature, pressure, virial, density

  enddo
  
  close(10)
  close(11)
  close(12)

  deallocate(forces, distances)

end program noblegas_simulation
