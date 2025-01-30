module functions
!!
!!
!!
public ::  lecture                 ! Parse and store atomic positions from an .xyz file
public ::  get_distances           ! Compute interatomic distances
public ::  get_forces              ! Derive forces based on atom positions
public ::  get_kinetic             ! Compute kinetic energy of the system
public ::  get_potential           ! Compute potential energy of the system
public ::  get_statistic           ! Compute new T,P,vir
public ::  initialize_velocities 

contains

  subroutine lecture(filename, positions, atom, nstep, box_size, temp, mass, timestep, sig, eps, natom)
  !!
  !!
  !!
  implicit none

  integer                                         :: i, nstep, natom
  double precision                                :: temp, mass, timestep, sig, eps
  double precision                                :: box_size
  double precision, allocatable, dimension(:,:)   :: positions
  character(len=2), allocatable, dimension(:)     :: atom
  character(len=100), intent(in)                  :: filename 

  open(1,file=filename,status="old")

  read(1,*) 
  read(1,*) nstep
 
  read(1,*)
  read(1,*) timestep

  read(1,*)
  read(1,*) natom

  ! Allocate arrays

  allocate(positions(natom,3))
  allocate(atom(natom))

  read(1,*)
  read(1,*) mass
  
  read(1,*)
  read(1,*) sig

  read(1,*)
  read(1,*) eps

  read(1,*)
  read(1,*) temp


  read(1,*)
  read(1,*) box_size

  read(1,*)
  do i=1,natom
     read(1,*) atom(i), positions(i,:)
     print*, positions(i,:)
  enddo

  close(1)

  end subroutine lecture


  subroutine get_distances(positions, distances, box_size, natom)
  !!
  !! Compute the distances between all pairs of atoms with PBC.
  !!
  implicit none

  integer                       :: i, j
  double precision              :: dx, dy, dz, distance
  integer, intent(in)           :: natom
  double precision, intent(in)  :: positions(natom, 3), box_size
  double precision, intent(out) :: distances(natom, natom)

  do i = 1, natom
      do j = 1, natom

        distances(i,j) = 0.0d0
                  
        dx = positions(i, 1) - positions(j, 1)
        dy = positions(i, 2) - positions(j, 2)
        dz = positions(i, 3) - positions(j, 3)

        if (abs(dx) > (box_size / 2.0d0)) then
            dx = dx - box_size * anint(dx/box_size)
        endif

        if (abs(dy) > (box_size / 2.0d0)) then
            dy = dy - box_size * anint(dy/box_size)
        endif

        if (abs(dz) > (box_size / 2.0d0)) then
            dz = dz - box_size * anint(dz/box_size)
        endif

        distances(i,j) = sqrt(dx**2 + dy**2 + dz**2)

    enddo
  enddo

  end subroutine get_distances


  subroutine get_forces(forces, virial, distances, positions, sig, eps, natom)
  !!
  !! Compute Lennard-Jones forces between atoms using precomputed distances.
  !!
  implicit none 

  integer                         :: i, j, k
  double precision                :: cutoff, force_mag, delta(3)
  integer, intent(in)             :: natom
  double precision, intent(in)    :: sig, eps, distances(natom, natom), positions(natom, 3)
  double precision, intent(out)   :: forces(natom, 3), virial
  
  cutoff = 2.5d0 * sig
  
  forces = 0.0d0
  virial = 0.0d0

  do i = 1, natom
    do j = i + 1, natom 
      if (distances(i, j) <= cutoff .and. distances(i, j) > 0.0d0) then
  
        force_mag = 48.0d0 * eps * (sig**12 / distances(i, j)**13 - 0.5d0 * sig**6 / distances(i, j)**7)
  
        do k = 1, 3

          delta(k) = positions(i, k) - positions(j, k)
  
          forces(i, k) = forces(i, k) + force_mag * delta(k) / distances(i, j)
          forces(j, k) = forces(j, k) - force_mag * delta(k) / distances(i, j)

        enddo
  
        virial = virial + force_mag * distances(i, j)

      endif
    enddo
  enddo

  endsubroutine get_forces
  

  subroutine get_kinetic(velocities,natom,kin_ener,mass)
  !!
  !! Calculate the system's kinetic energy
  !!
  implicit none

  integer                          :: i,k
  integer, intent(in)              :: natom
  double precision, intent(in)     :: velocities(natom,3), mass
  double precision, intent(out)    :: kin_ener


  kin_ener=0.0d0

  do i=1,natom                      
    do k=1,3                       

      kin_ener = kin_ener + 0.5d0 * mass * velocities(i,k)*velocities(i,k)

    enddo
  enddo

  end subroutine get_kinetic


  subroutine get_potential(distances, natom, sig, eps, pot_ener)
  !!
  !! Compute Lennard-Jones forces between atoms using precomputed distances.
  !!
  implicit none

  integer                         :: i, j
  double precision                :: cutoff, lj_potential, cutoff_potential
  integer, intent(in)             :: natom
  double precision, intent(in)    :: sig, eps, distances(natom, natom)
  double precision, intent(out)   :: pot_ener
  
  cutoff = 2.5d0 * sig
  
  pot_ener = 0.0d0
  
  cutoff_potential = 4.d0 * eps * ((sig / cutoff)**12 - (sig / cutoff)**6)
  
  do i = 1, natom - 1
    do j = i + 1, natom
      if (distances(i, j) <= cutoff .and. distances(i, j) > 0.0d0) then
        lj_potential = 4.d0 * eps * ((sig / distances(i, j))**12 - (sig / distances(i, j))**6) - cutoff_potential
  
        pot_ener = pot_ener + lj_potential

      endif
    enddo
  enddo

  end subroutine get_potential


  subroutine get_statistic(kin_ener, natom, volume, density, virial, temperature, pressure)
  !!
  !! Calculate the system's temperature and pressure.
  !!
  implicit none
  
  double precision                :: dof, kb_Hartree, scale
  integer, intent(in)             :: natom
  double precision, intent(in)    :: kin_ener, volume, density, virial
  double precision, intent(out)   :: temperature, pressure

  kb_Hartree = 3.166811563d-6
  dof = 3.0d0 * natom

  temperature = (2.d0 * kin_ener) / (dof * kb_Hartree)

  pressure = density * kb_Hartree *temperature - (virial / (3.0d0* volume))

  end subroutine get_statistic


  subroutine initialize_velocities(velocities, natom, temp, mass)
  !!
  !!
  !!
  implicit none
  
  integer                         :: i,k
  double precision                :: sig, pi2, x1, x2, x3 ,x4, kb_Hartree, sumv(3)
  double precision                :: ek, scale, velocity
  integer, intent(in)             :: natom  
  double precision, intent(in)    :: temp, mass   
  double precision, intent(out)   :: velocities(natom, 3)  
  
  ! Costanti
  
  kb_Hartree = 3.166811563d-6          
  pi2 = 2.0d0 * 3.141592653589793d0 
  
  sig = sqrt(3.0d0 * temp * kb_Hartree / mass) 
  
  call random_seed()
  sumv = 0.0d0
  ek = 0.0d0
  scale = 0.0d0
  do i = 1, natom
  
    call random_number(x1)
    call random_number(x2)
    call random_number(x3)
    call random_number(x4)
  
    velocities(i,1) = sig * sqrt(-2.0d0 * log(1.0d0 - x1)) * cos(pi2 * x2)
    velocities(i,2) = sig * sqrt(-2.0d0 * log(1.0d0 - x1)) * sin(pi2 * x2) 
    velocities(i,3) = sig * sqrt(-2.0d0 * log(1.0d0 - x3)) * cos(pi2 * x4)
  
    sumv(1) = sumv(1) + velocities(i,1)  
    sumv(2) = sumv(2) + velocities(i,2)  
    sumv(3) = sumv(3) + velocities(i,3)  
  
  enddo
 
  do i = 1, natom
    do k = 1, 3
    velocities(i,k) = velocities(i,k) - sumv(k)
    enddo
  enddo
  
  do i = 1, natom
    do k = 1, 3
    ek = ek + 0.5d0 * mass * velocities(i,k)*velocities(i,k)
    enddo
  enddo
  
  scale = sqrt((3.0d0 * temp * natom * kb_Hartree)/(2.0d0 * ek))

  do i = 1, natom
    do k = 1, 3
  
    velocities(i,k) = velocities(i,k) * (scale)
  
    velocity = sqrt(velocities(i,1)**2 + velocities(i,2)**2 + velocities(i,3)**2)
  
    print*, "velocity", velocity
    enddo
  enddo
  
  end subroutine initialize_velocities

end module functions
