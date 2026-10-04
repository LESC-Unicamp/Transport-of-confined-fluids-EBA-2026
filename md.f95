!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code to execute a NVE MD simulation with velocity-Verlet algorithm  !
! for LJ particles confined by slit pores in the z direction          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! University of Campinas                                              !
! School of Chemical Engineering                                      !
! Prof. Luis Fernando Mercier Franco                                  !
! Updated: Sep. 27th, 2026                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disclaimer:                                                         !
! The author does not accept any liability for the use of this code   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module globalvar
        ! Number of particles
        integer*8           :: n_particles
        ! Maximum number of steps in time integration
        integer*8           :: max_steps
        integer*8           :: n_save
        ! Number of equilibration steps
        integer*8           :: n_equil
        ! Avogadro's number in 1/mol
        real*8, parameter   :: avogadro  = 6.0221409d23
        ! Boltzmann constant in J/K
        real*8, parameter   :: boltzmann = 1.38064852d-23
        ! Universal gas constant in J/mol/K
        real*8, parameter   :: r         = avogadro*boltzmann
        ! Molar mass in g/mol
        real*8              :: molar_mass
        ! Density in kg/m³
        real*8              :: density
        ! Temperature in K
        real*8              :: temperature
        ! Potential energy in J/mol
        real*8              :: potential
        ! Auxiliary variable for computing the intermolecular potential energy
        real*8              :: v
        ! Auxiliary variable for computing the external potential energy
        real*8              :: vext
        ! Auxiliary variable for computing the pair potential energy
        real*8              :: vij
        ! Auxiliary variable for computing the external force
        real*8              :: fext
        ! Auxiliary variable for computing the intermolecular virial
        real*8              :: w
        ! Auxiliary variable for computing the pair virial
        real*8              :: wij
        ! Position in the z coordinate in Angstroms
        real*8              :: rzi
        ! Relative distance between i and j squared in Angstroms²
        real*8              :: rijsq
        ! Box length in the x coordinate in Angstroms
        real*8              :: box_length_x
        ! Box length in the y coordinate in Angstroms
        real*8              :: box_length_y
        ! Box length in the z coordinate in Angstroms
        real*8              :: box_length_z
        ! Time step in fs
        real*8              :: time_step
        real*8              :: eps
        real*8              :: epsfw
        real*8              :: eps24
        real*8              :: epsfw203
        real*8              :: sigma
        real*8              :: sigmasq
        real*8              :: rcut
        real*8              :: rcutsq
        ! x coordinates in Angstroms
        real*8, allocatable :: rx(:)
        ! y coordinates in Angstroms
        real*8, allocatable :: ry(:)
        ! z coordinates in Angstroms
        real*8, allocatable :: rz(:)
        ! x velocity in Angstroms/fs
        real*8, allocatable :: vx(:)
        ! y velocity in Angstroms/fs
        real*8, allocatable :: vy(:)
        ! z velocity in Angstroms/fs
        real*8, allocatable :: vz(:)
        ! x acceleration in Angstroms/fs²
        real*8, allocatable :: ax(:)
        ! y acceleration in Angstroms/fs²
        real*8, allocatable :: ay(:)
        ! z acceleration in Angstroms/fs²
        real*8, allocatable :: az(:)
        ! x force per mass in Angstroms/fs²
        real*8, allocatable :: fx(:)
        ! y force per mass in Angstroms/fs²
        real*8, allocatable :: fy(:)
        ! z force per mass in Angstroms/fs²
        real*8, allocatable :: fz(:)

        contains

        subroutine init_var()
               implicit none
               character :: comment*19


               open(1,file="file.inp")
               read(1,'(a19,i6)') comment,n_particles
               read(1,'(a19,f12.3)') comment,temperature
               read(1,'(a19,f12.3)') comment,box_length_x
               read(1,'(a19,f12.3)') comment,box_length_y
               read(1,'(a19,f12.3)') comment,box_length_z
               read(1,'(a19,f12.3)') comment,molar_mass
               read(1,'(a19,f12.3)') comment,sigma
               read(1,'(a19,f12.3)') comment,rcut
               read(1,'(a19,f12.3)') comment,eps
               read(1,'(a19,f12.3)') comment,epsfw
               read(1,'(a19,f12.3)') comment,time_step
               read(1,'(a19,i10)') comment,max_steps
               read(1,'(a19,i10)') comment,n_save
               read(1,'(a19,i10)') comment,n_equil
               close(1)

               allocate(rx(n_particles))
               allocate(ry(n_particles))
               allocate(rz(n_particles))
               allocate(vx(n_particles))
               allocate(vy(n_particles))
               allocate(vz(n_particles))
               allocate(ax(n_particles))
               allocate(ay(n_particles))
               allocate(az(n_particles))
               allocate(fx(n_particles))
               allocate(fy(n_particles))
               allocate(fz(n_particles))

               ! Converting from g/mol to kg/mol
               molar_mass     = molar_mass*1.d-3

               ! Epsilon in J/mol
               eps24          = 24.d-10*eps*r/molar_mass
               epsfw203       = 20.d0*epsfw*r/3.0

               ! Sigma in Angstroms
               sigmasq        = sigma*sigma

               ! Cut-off radius in Angstroms
               rcutsq         = rcut*rcut
        
        end subroutine init_var

end module globalvar   
      
      
program md
         use globalvar
         implicit none
         integer*8           :: i
         integer*8           :: steps
         ! Kinetic energy in J/mol
         real*8              :: kinetic
         ! Total energy in J/mol
         real*8              :: total_energy
         character           :: atom*1

         call init_var()

         ! Reading initial configuration
         open(1,file="conf.xyz")
         read(1,*) n_particles
         do i=1,n_particles
            read(1,*) atom,rx(i),ry(i),rz(i),vx(i),vy(i),vz(i)
         end do
         close(1)
        
         ! Computing kinetic energy 
         kinetic     = 0.5d10*molar_mass*sum(vx(:)**2.d0+vy(:)**2.d0   &
                       +vz(:)**2.d0)/dble(n_particles)
         ! Computing temperature
         temperature = 2.d0*kinetic/r/3.d0

         ! Computing accelerations
         call compute_acceleration()
    
         ! Computing total energy
         total_energy = kinetic+potential

         open(2,file="traj.xyz")
         open(3,file="thermo.dat")
         write(*,*) '    Steps       T(K)        K(J/mol)',&
                    '     U(J/mol)   E(J/mol)'
         ! Time integration
         do steps=1,max_steps
    
           write(*,'(i10,4f12.2)') steps,temperature,kinetic, &
                                   potential,total_energy
           write(3,'(i10,4e15.7)') steps,temperature,kinetic, &
                                   potential,total_energy

           ! Velocity-Verlet algorithm
           vx(:) = vx(:)+0.5d0*ax(:)*time_step
           vy(:) = vy(:)+0.5d0*ay(:)*time_step
           vz(:) = vz(:)+0.5d0*az(:)*time_step
           rx(:) = rx(:)+vx(:)*time_step
           ry(:) = ry(:)+vy(:)*time_step
           rz(:) = rz(:)+vz(:)*time_step
           call compute_acceleration()
           vx(:) = vx(:)+0.5d0*ax(:)*time_step
           vy(:) = vy(:)+0.5d0*ay(:)*time_step
           vz(:) = vz(:)+0.5d0*az(:)*time_step

           kinetic      = 0.5d10*molar_mass*sum(vx(:)**2.d0+vy(:)**2.d0&
                          +vz(:)**2.d0)/dble(n_particles)
           temperature  = 2.d0*kinetic/r/3.d0
           total_energy = kinetic+potential

           ! Storing trajectory every n_save steps
           if (mod(steps,n_save) == 0 .and. steps > n_equil) then
              write(2,*) n_particles
              write(2,*) ''
              do i=1,n_particles
                 write(2,*) 'C',rx(i),ry(i),rz(i)
              end do
           end if
          
    end do
    close(2)
    close(3)

end program md

subroutine compute_acceleration()
        use globalvar
        implicit none
        integer*8 :: i
        integer*8 :: j
        real*8    :: rxi
        real*8    :: ryi
        real*8    :: rxj
        real*8    :: ryj
        real*8    :: rzj
        real*8    :: rxij
        real*8    :: ryij
        real*8    :: rzij
        real*8    :: aij
        real*8    :: axi
        real*8    :: ayi
        real*8    :: azi
        real*8    :: fxi
        real*8    :: fyi
        real*8    :: fzi
        real*8    :: fij
        real*8    :: external_potential
    
        ax(:) = 0.d0
        ay(:) = 0.d0
        az(:) = 0.d0
        fx(:) = 0.d0
        fy(:) = 0.d0
        fz(:) = 0.d0
        v     = 0.d0
    
        do i=1,n_particles-1
           rxi = rx(i)
           ryi = ry(i)
           rzi = rz(i)
           fxi = fx(i)
           fyi = fy(i)
           fzi = fz(i)
           do j=i+1,n_particles
              rxj   = rx(j)
              ryj   = ry(j)
              rzj   = rz(j)
              rxij  = rxi-rxj
              ryij  = ryi-ryj
              rzij  = rzi-rzj
              ! Minimum image convention
              rxij  = rxij-box_length_x*dnint(rxij/box_length_x)
              ryij  = ryij-box_length_y*dnint(ryij/box_length_y)
              rijsq = rxij*rxij+ryij*ryij+rzij*rzij

              if (rijsq <= rcutsq) then
                 call compute_potential()
                 v     = v+vij
                 fij   = wij/rijsq
                 fxi   = fxi+fij*rxij
                 fyi   = fyi+fij*ryij
                 fzi   = fzi+fij*rzij
                 fx(j) = fx(j)-fij*rxij
                 fy(j) = fy(j)-fij*ryij
                 fz(j) = fz(j)-fij*rzij
              end if
           end do
           fx(i) = fxi
           fy(i) = fyi
           fz(i) = fzi
        end do

        ax(:)              = eps24*fx(:)
        ay(:)              = eps24*fy(:)
        external_potential = 0.0
        do i=1,n_particles
           rzi = rz(i)
           call compute_external_potential()
           az(i)              = eps24*fz(i)+1d-10*epsfw203*fext/molar_mass
           external_potential = external_potential+vext
        end do

        potential = 4.d0*eps*v/dble(n_particles)   &
                    +epsfw203*external_potential/dble(n_particles) 

end subroutine compute_acceleration

subroutine compute_potential()
        ! Subroutine to compute LJ potential and virial
        use globalvar
        implicit none
        real*8 :: sr2
        real*8 :: sr6
        real*8 :: sr12

        sr2  = sigmasq/rijsq
        sr6  = sr2*sr2*sr2
        sr12 = sr6*sr6
        vij  = sr12-sr6
        wij  = 2.d0*sr12-sr6

end subroutine compute_potential

subroutine compute_external_potential()
        ! Subroutine to compute solid-fluid potential and force
        use globalvar
        implicit none
        real*8 :: rz2
        real*8 :: sr2
        real*8 :: sr4
        real*8 :: sr10
        real*8 :: vrz
        real*8 :: vrh
        real*8 :: frz
        real*8 :: fh

        ! Potential and force from left wall (z = z0)
        rz2  = rzi*rzi
        sr2  = 0.25d0*sigmasq/rz2
        sr4  = sr2*sr2
        sr10 = sr4*sr4*sr2
        vrz  = 0.1d0*sr10-0.25d0*sr4
        frz  = (sr10-sr4)/rzi

        ! Potential and force from right wall (z = H)
        rz2  = (box_length_z-rzi)*(box_length_z-rzi)
        sr2  = 0.25d0*sigmasq/rz2
        sr4  = sr2*sr2
        sr10 = sr4*sr4*sr2
        vrh  = 0.1d0*sr10-0.25d0*sr4
        fh   = (sr10-sr4)/(box_length_z-rzi)

        vext = vrz+vrh

        fext = frz-fh

end subroutine compute_external_potential
