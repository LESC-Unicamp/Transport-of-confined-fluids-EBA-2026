!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Code to compute the survival probability                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! University of Campinas                                              !
! School of Chemical Engineering                                      !
! Prof. Luis Fernando Mercier Franco                                  !
! Updated: Sep. 27th, 2026                                            !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Disclaimer:                                                         !
! The author does not accept any liability for the use of this code   !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program compute_survival_probability
        implicit none
        integer*8                            :: i,j
        integer*8                            :: time_origin
        integer*8                            :: step
        integer*8                            :: time
        integer*8                            :: n_particles
        integer*8                            :: max_steps
        integer*8                            :: n_save
        integer*8                            :: n_equil
        integer*8                            :: maxbin
        integer*8                            :: n_part_zero
        real*8                               :: time_step
        real*8                               :: probability
        real*8                               :: temperature
        real*8                               :: molar_mass
        real*8                               :: eps
        real*8                               :: epsfw
        real*8                               :: sigma
        real*8                               :: rcut
        real*8                               :: box_length_x
        real*8                               :: box_length_y
        real*8                               :: box_length_z
        real*8                               :: rx
        real*8                               :: ry
        real*8, dimension(:,:), allocatable  :: rz
        real*8                               :: rzi
        real*8                               :: lower
        real*8                               :: upper
        real*8                               :: dummy
        real*8, dimension(:), allocatable    :: prob
        character                            :: atom*1
        character                            :: comment*19
        logical, dimension(:), allocatable   :: in_old
        logical, dimension(:), allocatable   :: in_new

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
        read(1,'(a19,i10)') comment,maxbin
        read(1,'(a19,f12.3)') comment,lower
        read(1,'(a19,f12.3)') comment,upper
        close(1)

        allocate(in_old(n_particles))
        allocate(in_new(n_particles))

        max_steps      = (max_steps-n_equil)/n_save

        allocate(prob(max_steps))

        prob(:)        = 0.d0

        allocate(rz(n_particles,max_steps))

        open(1,file="traj.xyz")
        do step=1,max_steps
           write(*,*) 'Reading step = ',step
           read(1,*) n_particles
           do i=1,n_particles
              read(1,*) atom,rx,ry,rz(i,step)
           end do
        end do
        close(1)

        do time_origin=1,max_steps/2

           write(*,*) 'Computing at time origin = ',dble(time_origin)*dble(n_save)*time_step,' fs'

           in_old(:)      = .false.
           n_part_zero    = 0
           do i=1,n_particles
              rzi = rz(i,time_origin)
              if (rzi >= lower .and. rzi <= upper) then
                 n_part_zero = n_part_zero+1
                 in_old(i)   = .true.
              end if
           end do

           do step=1,max_steps/2

              do i=1,n_particles
                 rzi = rz(i,step+time_origin-1)
                 if (rzi >= lower .and. rzi <= upper) then
                    in_new(i) = .true.
                    if (in_new(i) .and. in_old(i)) then
                       prob(step) = prob(step)+1.d0/n_part_zero
                    end if
                 else 
                    in_old(i) = .false.
                 end if
              end do

           end do

        end do

        deallocate(rz)

        open(1,file="survival_probability.dat")
        do step=1,max_steps/2
           write(1,*) dble(step-1)*dble(n_save)*time_step,2.d0*prob(step)/dble(max_steps)
        end do
        close(1)

end program compute_survival_probability
