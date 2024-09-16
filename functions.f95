module functions

  public::lecture                 ! Parse and store atomic positions from an .xyz file
  public::get_distances           ! Compute interatomic distances
  public::get_forces              ! Derive forces based on atom positions
  public::get_kinetic             ! Compute kinetic energy of the system
  public::get_potential           ! Compute potential energy of the system

  contains

! **********************************************************

  subroutine lecture (filename, positions, atom, velocities, nstep, mass, timestep, sig, eps, natom)

! Load input parameters from file

  implicit none

  integer                                        :: i , nstep, natom
  double precision                               :: timestep, mass, eps, sig
  double precision, allocatable, dimension(:,:)  :: positions, velocities
  character(len=100)                             :: filename
  character(len=2), allocatable, dimension(:)    :: atom

  open(1,file=filename,status="old")

  read(1,*) 
  read(1,*) nstep
  
  read(1,*)
  read(1,*) timestep

  read(1,*)
  read(1,*) natom

  ! Reserve memory for arrays

  allocate(positions(natom,3))
  allocate(atom(natom))
  allocate(velocities(natom,3))

  read(1,*)
  read(1,*) mass
  
  read(1,*)
  do i=1,natom
     read(1,*) atom(i), positions(i,:)
  enddo

  read(1,*)
  read(1,*) sig

  read(1,*)
  read(1,*) eps

  read(1,*)
  do i=1,natom
     read(1,*) velocities(i,:)
  enddo

  ! Close the input file
  close(1)

  end subroutine lecture

! **********************************************************

  subroutine get_distances (positions, distances, natom)

! Compute the distances between all pairs of atoms

  implicit none

  integer, intent(in)           :: natom
  double precision, intent(in)  :: positions(natom,3)
  double precision, intent(out) :: distances(natom,natom)

  integer                       :: i, j

!  Set initial distances to zero
  distances=0.0d0

  do i = 1, natom
      do j = 1, natom
        distances(i,j) = dsqrt((positions(i,1) - positions(j,1))**2 + & 
                             & (positions(i,2) - positions(j,2))**2 + & 
                             & (positions(i,3) - positions(j,3))**2)
      enddo
  enddo

  end subroutine get_distances

! **********************************************************

  subroutine get_forces(forces, distances, positions, sig, eps, natom)

! Compute Lennard-Jones forces between atoms

  implicit none

  integer, intent(in)             :: natom
  double precision, intent(in)    :: eps, sig
  double precision, intent(inout) :: forces(natom,3)
  double precision, intent(in)    :: positions(natom,3), distances(natom,natom)

  integer                         :: i, j, k


! Initialize forces to zero
  forces = 0.d0

  do i = 1, natom  ! Iterate over the atom being acted upon
    do k = 1,3  ! Iterate over the spatial coordinates (x, y, z)
        do j = 1, natom  ! Iterate over the atoms applying the force
          if (i == j) cycle  ! Skip self-interaction
          forces(i,k) = forces(i,k) + 48*eps*(positions(i,k)-positions(j,k)) * & 
                                    & (sig**12/distances(i,j)**14 - 0.5d0*sig**6/distances(i,j)**8) 
        enddo
    enddo
  enddo

  end subroutine get_forces

! **********************************************************

  subroutine get_kinetic(velocities,natom,kin_ener,mass)

! Calculate the system's kinetic energy

  implicit none

  double precision, intent(in)     :: velocities(natom,3), mass
  double precision, intent(out)    :: kin_ener
  integer, intent(in)              :: natom
  integer                          :: i,k

! Initialize the kinetic energy accumulator
  kin_ener=0.0d0

  do i=1,natom  ! Loop over all atoms
    do k=1,3  ! Loop over x, y, z components of velocity
      kin_ener = kin_ener + 0.5d0 * mass * velocities(i,k)*velocities(i,k)
    enddo
  enddo

  end subroutine get_kinetic

! **********************************************************

  subroutine get_potential (distances,natom,sig,eps,pot_ener)

! Calculate the system's potential energy

  implicit none

  double precision, intent(in)  :: distances(natom,natom), sig, eps
  double precision, intent(out) :: pot_ener
  integer, intent(in)           :: natom
  integer                       :: i,j

! Initialize potential energy to zero
  pot_ener = 0.0d0

! Loop over distinct atom pairs (i, j), skipping double counting
  do i = 1,natom-1
     do j = i+1,natom
        pot_ener = pot_ener + 4.d0 * eps * ((sig/distances(i,j))**12-(sig/distances(i,j))**6)
     enddo
  enddo

  end subroutine get_potential

end module functions
