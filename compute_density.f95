!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code to compute the number density profile in z direction           !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! University of Campinas                                              !
! School of Chemical Engineering                                      !
! Prof. Luis Fernando Mercier Franco                                  !
! Updated: Sep. 27th, 2026                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disclaimer:                                                         !
! The author does not accept any liability for the use of this code   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program compute_density 
        implicit none
        integer*8                            :: i,j
        integer*8                            :: bin
        integer*8                            :: step
        integer*8                            :: n_particles
        integer*8                            :: max_bin 
        integer*8                            :: max_steps
        integer*8                            :: n_save
        integer*8                            :: n_equil
        integer*8, dimension(:), allocatable :: hist
        real*8                               :: time_step
        real*8                               :: density
        real*8                               :: temperature
        real*8                               :: molar_mass
        real*8                               :: eps
        real*8                               :: epsfw
        real*8                               :: sigma
        real*8                               :: rcut
        real*8                               :: box_length_x
        real*8                               :: box_length_y
        real*8                               :: box_length_z
        real*8                               :: delz
        real*8                               :: pi
        real*8, dimension(:), allocatable    :: rx
        real*8, dimension(:), allocatable    :: ry
        real*8, dimension(:), allocatable    :: rz
        real*8                               :: rzi
        real*8                               :: zlower
        real*8                               :: dummy
        character                            :: atom*1
        character                            :: comment*19

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
        read(1,'(a19,i10)') comment,max_bin
        close(1)

        allocate(hist(max_bin))

        max_steps  = (max_steps-n_equil)/n_save
        delz       = box_length_z/dble(max_bin)
        hist(:)    = 0

        allocate(rx(n_particles))
        allocate(ry(n_particles))
        allocate(rz(n_particles))

        open(1,file="traj.xyz")
        do step=1,max_steps
           read(1,*) n_particles
           do i=1,n_particles
             read(1,*) atom,rx(i),ry(i),rz(i)
           end do

           do i=1,n_particles
              rzi = rz(i)
              bin = dnint(rzi/delz)
              if (bin <= max_bin) then
                 hist(bin) = hist(bin)+1
              end if
           end do
        end do
        close(1)

        deallocate(rx,ry,rz)

        open(1,file="density.dat")
        do bin=1,max_bin
           zlower  = delz*dble(bin)
           density = dble(hist(bin))/dble(max_steps)/dble(n_particles)
           write(1,*) zlower+0.5d0*delz,density
        end do
        close(1)

end program compute_density
