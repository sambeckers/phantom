!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to generate data for chemistry
! Tracks particles paths and store the evolution in one file per particle
!
! :References: None
!
! :Owner: Camille Landri
!
! :Runtime parameters: None
!
! :Dependencies: None
!

 use part,       only: maxp
 use raytracer,  only: get_all_tau
 implicit none
 character(len=20), parameter, public :: analysistype = 'trace'
 public :: do_analysis

 real, allocatable    :: one(:)
 integer(8), allocatable :: iorig_old(:)
 integer, allocatable :: iprev(:)
 logical :: done_init = .false.
 real :: AuvAv = 4.65, albedo = 0.5

 integer :: ntrack = 0
 integer, allocatable :: track_id(:)
 logical, allocatable :: mask(:)
 character(len=256) :: dir

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only: isdead_or_accreted, iorig, rhoh, nptmass, xyzmh_ptmass, iReff, iboundary, igas, iphase, iamtype, maxp
 use linklist,   only: set_linklist
 use units,      only: utime,unit_density,udist
 use physcon,    only: atomic_mass_unit
 use eos,        only: get_temperature, ieos, gamma,gmw, init_eos
 use io,         only: fatal, iverbose

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 integer, save :: nprev = 0
 real          :: rho_cgs, numberdensity, T_gas, gammai, mui, AUV
 real          :: column_density(npart), xyzh_copy(4,npart)
 real          :: max_radius, radius
 integer       :: i, j, k, i_radius, completed_iterations, npart_copy = 0
 integer       :: iu=10,ios, ierr
 logical       :: iexist
 character(len=9) :: filename
 integer :: isize

 if (.not.done_init) then
      allocate(one(npart))
      one = 1.0
      allocate(iorig_old(maxp))
      allocate(iprev(maxp))
      done_init = .true.

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

      ! make directory to store individual particle trace files in dumpfile directory
      ! get path of dumpfile
      dir = trim(dir)//'trace_output/'
      print *, "Creating directory for tracked particles output in ", dir
      inquire(file=dir, exist=iexist)
      if (.not.iexist) then
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
      call init_eos(ieos, ierr)
      if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")
   
   else
      completed_iterations = 0
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
      !$omp shared(npart,xyzh,vxyzu,nprev,iorig,iorig_old,iprev,iverbose, dir) &
      !$omp shared(particlemass,unit_density, mask, time, utime, iphase) &
      !$omp shared(ieos,gamma,gmw,completed_iterations,column_density,AuvAv,albedo) &
      !$omp private(i,j,k,rho_cgs,numberdensity,T_gas,gammai,mui,AUV,filename,iu,isize)
      outer: do i=1,npart
         if (mask(i) .eqv. .true. .and. .not.isdead_or_accreted(xyzh(4,i))) then
            inner: do j=1,nprev
               if (iorig(i) == iorig_old(j)) then
                  iprev(i) = j
                  exit inner
               endif
            enddo inner
            if (iamtype(iphase(i)) /= iboundary .and. i > 2460) then ! 2460 is the amount of boundary particles
               !Thermodynamic quantities
               rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
               gammai = gamma
               mui    = gmw
               numberdensity = rho_cgs / (mui * atomic_mass_unit)
               T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gammai,mui)
               T_gas = max(T_gas,20.0d0) ! Floor temperature at 20K

               !Radiation quantities
               AUV = AuvAv * column_density(i) / (mui * atomic_mass_unit) / 1.87e21

            !$omp critical
            write(filename, '(i9)') iorig(i)
            inquire(file=trim(dir)//trim(adjustl(filename))//'.phys', size=isize)
            if (isize == -1) then
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.phys', status='new', action='write')
               print *, 'Creating new file for particle ', iorig(i)
               write(iu, *) '# time(s)   x(AU)   Y(AU)   Z(AU)   n(cm-3)   T(K)   A_UV'
            else if (isize == 0) then
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.phys', status='old', action='write')
               print *, 'Filling empty file for particle ', iorig(i)
               write(iu, *) '# time(s)   x(AU)   Y(AU)   Z(AU)   n(cm-3)   T(K)   A_UV'
            else
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.phys', status='old', action='write', position='append')
            endif 
               ! write physical parameters to file
               write(iu, *) time*utime, xyzh(1, i), xyzh(2, i), xyzh(3, i), numberdensity, T_gas, AUV
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
 iorig_old(1:npart) = iorig(1:npart)
end subroutine do_analysis

end module analysis