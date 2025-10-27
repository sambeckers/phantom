!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to generate training data for MACE using KROME
! Tracks particles paths and compute the chemical evolution along the path
!
! :References: None
!
! :Owner: Camille Landri
!
! :Runtime parameters: None
!
! :Dependencies: None
!
 use krome_user, only: krome_nmols
 use part,       only: maxp
 use raytracer,  only: get_all_tau
 implicit none
 character(len=20), parameter, public :: analysistype = 'krome'
 public :: do_analysis

 real, allocatable    :: abundance(:,:), abundance_prev(:,:), one(:)
 character(len=16)    :: abundance_label(krome_nmols)
 integer(8), allocatable :: iorig_old(:)
 integer, allocatable :: iprev(:)
 logical :: done_init = .false.
 real :: AuvAv = 4.65, albedo = 0.5

 integer :: ntrack = 0
 integer, allocatable :: track_id(:), mask(:)
 character(len=256) :: dir

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only: isdead_or_accreted, iorig, rhoh, nptmass, xyzmh_ptmass, iReff, iboundary, igas, iphase, iamtype, maxp
 use linklist,   only: set_linklist
 use units,      only: utime,unit_density,udist
 use physcon,    only: atomic_mass_unit
 use eos,        only: get_temperature, ieos, gamma,gmw, init_eos
 use io,         only: fatal, iverbose
 use krome_main, only: krome_init, krome
 use krome_user, only: krome_get_names,krome_set_user_Auv,krome_set_user_xi,&
                       krome_set_user_alb,krome_set_user_AuvAv
use krome_user, only: krome_idx_He,krome_idx_C,krome_idx_N,krome_idx_O,&
       krome_idx_H,krome_idx_S,krome_idx_Fe,krome_idx_Si,krome_idx_Mg,&
       krome_idx_Na,krome_idx_P,krome_idx_F,krome_idx_CO,krome_idx_C2H2,&
       krome_idx_C2H,krome_idx_H2,krome_idx_SiNC,krome_idx_e
 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real, save    :: tprev = 0.
 integer, save :: nprev = 0
 real          :: dt_cgs, rho_cgs, numberdensity, T_gas, gammai, mui, AUV, xi &
                , rholist(npart), Tlist(npart), mulist(npart), Auvlist(npart), xilist(npart)
 real          :: abundance_part(krome_nmols), Y(krome_nmols), column_density(npart), xyzh_copy(4,npart)
 real          :: max_radius, radius
 integer       :: i, j, k, i_radius, ierr, completed_iterations, npart_copy = 0
 integer       :: iu=10,ios
 character(len=9) :: filename
 integer :: isize

 if (.not.done_init) then
    ! find current directory
    dir = dumpfile
    do i=len_trim(dir),1,-1
       if (dir(i:i) == '/') then
          dir = dir(1:i)
          exit
       endif
    enddo

    ! read file with particle IDs to track and store IDs in array
    open(iu, file=trim(dir)//'particle_IDs.txt', status='old', action='read', iostat=ios)
    if (ios /= 0) call fatal(analysistype, "Could not open particle ID file")
    print*, "Reading particle IDs from particle_IDs.txt"

    ! make directory to store individual particle chemistry files in dumpfile directory
    ! get path from dumpfile
    dir = trim(dir)//'chem_output/'
    print *, "Creating directory for chemistry output in ", dir
    inquire(file=dir, exist=ios)
    if (ios == 0) then
       call system('mkdir '//dir)
    end if

    ! find number of tracked particles
    do
       read(iu,*,iostat=ios)
       if (ios /= 0) exit
       ntrack = ntrack + 1
    end do
    allocate(track_id(ntrack))
    allocate(mask(maxp))
    ! store particle IDs
    rewind(iu)
    do i=1,ntrack
       read(iu,*) track_id(i)
    end do
   close(iu) 
   print*, "Tracking ", ntrack, " particles"
   
   ! Initialise KROME and abundances
    done_init = .true.
    print*, "initialising KROME"
    call krome_init()
    print*, "Initialised KROME"
    abundance_label(:) = krome_get_names()
    allocate(abundance(krome_nmols,maxp))
    abundance = 0.
    allocate(abundance_prev(krome_nmols,maxp))
    abundance_prev = 0.
    allocate(one(maxp))
    one = 1.
    allocate(iorig_old(maxp))
    iorig_old = 0
    allocate(iprev(maxp))
    iprev = 0
    print*, "setting abundances"

   !$omp parallel do default(none) &
   !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev) &
   !$omp shared(abundance, abundance_prev,particlemass,unit_density) &
   !$omp shared(ieos,rho_cgs,T_gas,j) &
   !$omp private(i,abundance_part)
    do i=1, npart
      if (.not.isdead_or_accreted(xyzh(4,i))) then
          call chem_init(abundance_part)
          abundance(:,i) = abundance_part
       endif
    enddo
    call init_eos(ieos, ierr)
    if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")
 
 else
    dt_cgs = (time - tprev)*utime
    completed_iterations = 0
    print*, dumpfile, ": not first step data, timestep = ",dt_cgs, "npart = ",npart, "nprev = ",nprev
    xyzmh_ptmass(iReff,1) = 2.
    npart_copy = npart
    xyzh_copy = xyzh(:,:npart)
    call set_linklist(npart_copy,npart_copy,xyzh_copy,vxyzu)
    call get_all_tau(npart, nptmass, xyzmh_ptmass, xyzh, one, 5, .false., column_density)
    max_radius = 0.0
    do i = 1, npart
       if (.not.isdead_or_accreted(xyzh(4, i))) then
          radius = sqrt(xyzh(1, i)**2 + xyzh(2, i)**2 + xyzh(3, i)**2)
          if (radius > max_radius) then
             max_radius = radius
             i_radius = i
          endif
       endif
    enddo

    column_density = column_density + rhoh(xyzh(4,i_radius),particlemass)*unit_density * max_radius * udist
    ! update mask array to only compute chemistry for tracked particles
    mask = .false.
    do i=1,ntrack
       do j=1,npart
          if (iorig(j) == track_id(i)) then
             mask(j) = .true.
          endif
       enddo
    enddo
    !$omp parallel do default(none) &
    !$omp shared(npart,xyzh,vxyzu,dt_cgs,nprev,iorig,iorig_old,iprev,iverbose, dir) &
    !$omp shared(abundance,abundance_label,abundance_prev,particlemass,unit_density, mask, time, utime) &
    !$omp shared(ieos,gamma,gmw,completed_iterations,column_density,AuvAv,albedo) &
    !$omp shared(rholist,Tlist,mulist,Auvlist,xilist,iphase) &
    !$omp private(i,j,k,abundance_part,Y,rho_cgs,numberdensity,T_gas,gammai,mui,AUV,xi,radius,filename,iu,isize)
    outer: do i=1,npart
       if (mask(i) ==.true. .and. .not.isdead_or_accreted(xyzh(4,i))) then
          inner: do j=1,nprev
             if (iorig(i) == iorig_old(j)) then
                iprev(i) = j
                exit inner
             endif
          enddo inner
          if (j == iprev(i)) then
             abundance_part(:) = abundance_prev(:,iprev(i))
          else
             call chem_init(abundance_part)
          endif
          if (iamtype(iphase(i)) /= iboundary .and. i > 2460) then ! 2460 is the amount of boundary particles
             !Thermodynamic quantities
             rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
             gammai = gamma
             mui    = gmw
             numberdensity = rho_cgs / (mui * atomic_mass_unit)
             T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gammai,mui)
             T_gas = max(T_gas,20.0d0)
             radius = sqrt(xyzh(1, i)**2 + xyzh(2, i)**2 + xyzh(3, i)**2)

             !Radiation quantities
             AUV = AuvAv * column_density(i) / (mui * atomic_mass_unit) / 1.87e21
             xi = get_xi(AUV)
             call krome_set_user_Auv(AUV)
             call krome_set_user_xi(xi)
             call krome_set_user_alb(albedo)
             call krome_set_user_AuvAv(AuvAv)
             Y = abundance_part*numberdensity
             call krome(Y,T_gas,dt_cgs)
             abundance_part = Y/numberdensity
             abundance(:,i) = abundance_part
      
            !$omp critical
            write(filename, '(i9)') iorig(i)
            inquire(file=trim(dir)//trim(adjustl(filename))//'.chem', size=isize)
            if (isize == -1) then
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.chem', status='new', action='write')
               print *, 'Creating new file for particle ', iorig(i)
               write(iu, *) '# time(s)   radius(AU)   n(cm-3)   T(K)   mu   A_UV   xi   ', (abundance_label(k), k=1,krome_nmols) 
            else if (isize == 0) then
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.chem', status='old', action='write')
               print *, 'Filling empty file for particle ', iorig(i)
               write(iu, *) '# time(s)   radius(AU)   n(cm-3)   T(K)   mu   A_UV   xi   ', (abundance_label(k), k=1,krome_nmols)
            else
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.chem', status='old', action='write', position='append')
            endif 
             ! write physical parameters to file
             write(iu, '(ES16.8,1x,ES14.7,1x,ES14.7,1x,F8.2,1x,F6.3,1x,F8.3,1x,F8.3,1x)',advance="no") time*utime, radius, numberdensity, T_gas, mui, AUV, xi
             ! write abundances to file
            do k=1,krome_nmols
               write(iu, '(ES14.7,1x)', advance="no") abundance_part(k)
            enddo
            write(iu,*) ! new line
            close(iu)
            !$omp end critical
          endif

       endif
       if (iverbose > 1) then
          !$omp atomic
          completed_iterations = completed_iterations + 1
          print*, 'Completed ', completed_iterations, ' of ', npart
       endif


    enddo outer
 endif
   
 ! store current step data before moving on to next step
 nprev = npart
 tprev = time
 iorig_old(1:npart) = iorig(1:npart)
 abundance_prev(:,1:npart) = abundance(:,1:npart)
end subroutine do_analysis

real function get_xi(AUV)
 use physcon, only: pi
 real, intent(in) :: AUV
 real :: xi
 real :: W(6), GA(6), ceta
 integer :: i

 W(1) = 0.17132449
 W(2) = 0.36076157
 W(3) = 0.46791393
 W(4) = W(1)
 W(5) = W(2)
 W(6) = W(3)
 GA(1) = 0.93246951
 GA(2) = 0.66120939
 GA(3) = 0.23861919
 GA(4) = -GA(1)
 GA(5) = -GA(2)
 GA(6) = -GA(3)

 xi = 0.0
 do i=1,6
   ceta = (pi*GA(i)+pi)/2.0
   xi=xi+(W(i)*(sin(ceta)*exp((-AUV*ceta)/sin(ceta))))
 enddo
 xi = (pi/4.0)*xi
 
 get_xi = xi

end function get_xi


subroutine write_chem(npart, dumpfile, rholist, Tlist, mulist, Auvlist, xilist,mask)
 use krome_user, only: krome_idx_He,krome_idx_C,krome_idx_N,krome_idx_O,&
       krome_idx_H,krome_idx_S,krome_idx_Fe,krome_idx_Si,krome_idx_Mg,&
       krome_idx_Na,krome_idx_P,krome_idx_F,krome_idx_CO,krome_idx_C2H2,&
       krome_idx_C2H,krome_idx_H2,krome_idx_SiNC,krome_idx_e
 integer, intent(in)          :: npart
 character(len=*), intent(in) :: dumpfile
 real, intent(in)             :: rholist(*), Tlist(*), mulist(*), Auvlist(*), xilist(*)
 logical, intent(in)          :: mask(*)
 integer :: i, iu

 open(newunit=iu, file=dumpfile//'.comp', status='replace', action='write')
 write(iu, *) '# H, He, C, N, O, S, Fe, Si, Mg, Na, P, F, CO, C2H2, C2H, H2, SiNC, e-'
 do i=1, npart
   if (mask(i)) then
      write(iu, *) abundance(krome_idx_H, i),  abundance(krome_idx_He, i),   abundance(krome_idx_C, i),   &
                  abundance(krome_idx_N, i),  abundance(krome_idx_O, i),    abundance(krome_idx_S, i),   &
                  abundance(krome_idx_Fe, i), abundance(krome_idx_Si, i),   abundance(krome_idx_Mg, i),  &
                  abundance(krome_idx_Na, i), abundance(krome_idx_P, i),    abundance(krome_idx_F, i),   &
                  abundance(krome_idx_CO, i), abundance(krome_idx_C2H2, i), abundance(krome_idx_C2H, i), &
                  abundance(krome_idx_H2, i), abundance(krome_idx_SiNC, i), abundance(krome_idx_e, i), rholist(i), &
                  Tlist(i), mulist(i), Auvlist(i), xilist(i)
   endif
 enddo
 close(iu)
 
end subroutine write_chem

subroutine chem_init(abundance_part)
 use krome_user, only: krome_idx_H2,krome_idx_He,krome_idx_CO,krome_idx_C2H2,&
       krome_idx_HCN,krome_idx_N2,krome_idx_SiC2,krome_idx_CS,&
       krome_idx_SiS,krome_idx_SiO,krome_idx_CH4,krome_idx_H2O,&
       krome_idx_HCl,krome_idx_C2H4,krome_idx_NH3,krome_idx_HCP,&
       krome_idx_HF,krome_idx_H2S,krome_idx_e,krome_get_electrons
 real, intent(out) :: abundance_part(krome_nmols)

 ! Initial abundances for the krome model taken from Agúndez et al. (2020)
 ! H2, He, CO, C2H2, HCN, N2, SiC2, CS, SiS, SiO, CH4, H2O, HCl, C2H4, NH3, HCP, HF, H2S
 abundance_part(:)              = 0.
 abundance_part(krome_idx_H2)   = 0.5d0
 abundance_part(krome_idx_He)   = 8.5d-2
 abundance_part(krome_idx_CO)   = 4d-4
 abundance_part(krome_idx_C2H2) = 2.19d-5
 abundance_part(krome_idx_HCN)  = 2.045d-5
 abundance_part(krome_idx_N2)   = 2d-5
 abundance_part(krome_idx_SiC2) = 9.35d-6
 abundance_part(krome_idx_CS)   = 5.3d-6
 abundance_part(krome_idx_SiS)  = 2.99d-6
 abundance_part(krome_idx_SiO)  = 2.51d-6
 abundance_part(krome_idx_CH4)  = 1.75d-6
 abundance_part(krome_idx_H2O)  = 1.275d-6
 abundance_part(krome_idx_HCl)  = 1.625d-7
 abundance_part(krome_idx_C2H4) = 3.425d-8
 abundance_part(krome_idx_NH3)  = 3d-8
 abundance_part(krome_idx_HCP)  = 1.25d-8
 abundance_part(krome_idx_HF)   = 8.5d-9
 abundance_part(krome_idx_H2S)  = 2d-9
 abundance_part(krome_idx_e)    = krome_get_electrons(abundance_part(:))

end subroutine chem_init

end module analysis