!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code to generate an initial configuration of spherical particles    !
! randomly distributed with initial velocities distributed according  !
! to Maxwell-Boltzmann. Cubic simulation box, considering slit walls  !
! in z direction.                                                     !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! University of Campinas                                              !
! School of Chemical Engineering                                      !
! Prof. Luis Fernando Mercier Franco                                  !
! Updated: Sep 27th, 2026                                             !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disclaimer:                                                         !
! The author does not accept any liability for the use of this code   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program init_conf
        implicit none
        ! Counters 
        integer*8           :: i
        integer*8           :: j
        integer*8           :: k
        integer*8           :: cont
        integer*8           :: trial
        integer*8           :: particle
        ! Number of particles
        integer*8           :: n_particles
        ! Trial positions in x coordinate
        real*8              :: rxtrial
        ! Trial positions in y coordinate
        real*8              :: rytrial
        ! Trial positions in z coordinate
        real*8              :: rztrial
        ! Distances between i and j in x coordinate
        real*8              :: rxij
        ! Distances between i and j in y coordinate
        real*8              :: ryij
        ! Distances between i and j in z coordinate
        real*8              :: rzij
        ! Distances between i and j
        real*8              :: rij
        ! Box length in x coordinate in Angstroms
        real*8              :: box_length_x
        ! Box length in y coordinate in Angstroms
        real*8              :: box_length_y
        ! Box length in z coordinate in Angstroms
        real*8              :: box_length_z
        ! Absolute temperature in K
        real*8              :: temperature
        ! Molar mass in g/mol
        real*8              :: molar_mass
        ! Sum of velocities
        real*8              :: sumv
        ! Normalization factor
        real*8              :: norm
        ! Average velocity in x coordinate
        real*8              :: avg_vx
        ! Average velocity in y coordinate
        real*8              :: avg_vy
        ! Average velocity in z coordinate
        real*8              :: avg_vz
        ! Gaussian random number
        real*8              :: gauss
        ! Random number
        real*8              :: random_n
        ! Sigma from LJ potential in Angstroms
        real*8              :: sigma
        ! Conversion factor from meters to Angstroms
        real*8, parameter   :: mtoang    = 1.d10
        ! Avogadro's number in 1/mol
        real*8, parameter   :: avogadro  = 6.0221409d23
        ! Boltzmann constant in J/K
        real*8, parameter   :: boltzmann = 1.38064852d-23
        ! Universal gas constant in J/mol/K
        real*8, parameter   :: r         = avogadro*boltzmann
        ! x coordinate position in Angstroms
        real*8, allocatable :: rx(:)
        ! y coordinate position in Angstroms
        real*8, allocatable :: ry(:)
        ! z coordinate position in Angstroms
        real*8, allocatable :: rz(:)
        ! x velocity in Angstroms/fs
        real*8, allocatable :: vx(:)
        ! y velocity in Angstroms/fs
        real*8, allocatable :: vy(:)
        ! z velocity in Angstroms/fs
        real*8, allocatable :: vz(:)
        ! Text
        character           :: comment*18
        ! Logical variable to state if there is no overlap
        logical             :: no_overlap

        open(1,file="file.inp")
        read(1,'(a18,i6)') comment,n_particles
        read(1,'(a18,f12.3)') comment,temperature
        read(1,'(a18,f12.3)') comment,box_length_x
        read(1,'(a18,f12.3)') comment,box_length_y
        read(1,'(a18,f12.3)') comment,box_length_z
        read(1,'(a18,f12.3)') comment,molar_mass
        read(1,'(a18,f12.3)') comment,sigma
        close(1)

        allocate(rx(n_particles))
        allocate(ry(n_particles))
        allocate(rz(n_particles))
        allocate(vx(n_particles))
        allocate(vy(n_particles))
        allocate(vz(n_particles))


        ! Initial position coordinates at the center of the simulation box
        rx(:) = 0.5*box_length_x
        ry(:) = 0.5*box_length_y
        rz(:) = 0.5*box_length_z

        cont  = 1
        trial = 1
        k     = 1
        do while (cont <= n_particles)
           ! Stating no overlaps
           no_overlap = .true.
           ! Calling a random number evenly distributed between 0 and 1
           call random_number(random_n)
           ! Computing the trial x coordinate in Angstroms
           rxtrial = random_n*box_length_x
           ! Calling a random number evenly distributed between 0 and 1
           call random_number(random_n)
           ! Computing the trial y coordinate in Angstroms
           rytrial = random_n*box_length_y
           ! Calling a random number evenly distributed between 0 and 1
           call random_number(random_n)
           ! Computing the trial z coordinate in Angstroms, considering the slit walls.
           rztrial = random_n*(box_length_z-2.d0*sigma)+sigma
           ! Testing overlaps with other particles
           do particle=1,k
              rxij = rx(particle)-rxtrial
              ryij = ry(particle)-rytrial
              rzij = rz(particle)-rztrial
              rij  = dsqrt(rxij*rxij+ryij*ryij+rzij*rzij)
              if (rij <= 2.d0*sigma) then
                 no_overlap = .false.
              end if
           end do
           if (no_overlap) then
              write(*,*) 'After',trial,' trials',cont, &
                         'particles successfully inserted! ...'
              rx(k) = rxtrial
              ry(k) = rytrial
              rz(k) = rztrial
              k     = k+1
              cont  = cont+1
              trial = 1
           else
              trial = trial+1
           end if
        end do 

        ! Centralizing the box
        rx(:)  = rx(:)-0.5d0*box_length_x
        ry(:)  = ry(:)-0.5d0*box_length_y

        ! Converting molar mass from g/mol to kg/mol
        molar_mass = 1.d-3*molar_mass

        ! Generating velocities according Maxwell-Boltzmann distribution
        do particle=1,n_particles
           vx(particle)  = gauss()
           vy(particle)  = gauss()
           vz(particle)  = gauss()
        end do
        sumv   = dsqrt(sum(vx(:)**2.d0+vy(:)**2.d0+vz(:)**2.d0))
        norm   = 1d-5*dsqrt(3.d0*dble(n_particles)*r*temperature       &
                 /molar_mass)
        vx(:)  = vx(:)*norm/sumv
        vy(:)  = vy(:)*norm/sumv
        vz(:)  = vz(:)*norm/sumv
        avg_vx = sum(vx)/dble(n_particles)
        avg_vy = sum(vy)/dble(n_particles)
        avg_vz = sum(vz)/dble(n_particles)
        vx(:)  = vx(:)-avg_vx
        vy(:)  = vy(:)-avg_vy
        vz(:)  = vz(:)-avg_vz
       
        ! Writing initial configuration in conf.xyz 
        open(1,file="conf.xyz")
        write(1,*) n_particles
        write(1,*) ''
        do i=1,n_particles
           write(1,*) 'C',rx(i),ry(i),rz(i),vx(i),vy(i),vz(i)
        end do
        close(1)

end program init_conf

function gauss() result (fun)
        implicit none
        integer*4         :: i
        real*8, parameter :: a1  = 3.949846138d0
        real*8, parameter :: a3  = 0.252408784d0
        real*8, parameter :: a5  = 0.076542912d0
        real*8, parameter :: a7  = 0.008355968d0
        real*8, parameter :: a9  = 0.029899776d0
        real*8            :: summ 
        real*8            :: r
        real*8            :: r2
        real*8            :: fun
        real*8            :: random_n

        summ = 0.d0
        do i=1,12
           call random_number(random_n)
           summ = summ+random_n
        end do

        r   = 0.25d0*(summ-6.d0)
        r2  = r*r
        fun = ((((a9*r2+a7)*r2+a5)*r2+a3)*r2+a1)*r
        
end function gauss
