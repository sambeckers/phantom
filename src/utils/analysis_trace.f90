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
 implicit none
 character(len=20), parameter, public :: analysistype = 'trace'
 public :: do_analysis

 real, allocatable    :: one(:)
 integer(8), allocatable :: iorig_old(:)
 integer, allocatable :: iprev(:)
 logical :: done_init = .false.
 integer :: ntrack = 0
 integer, allocatable :: track_id(:)
 logical, allocatable :: mask(:)
 integer :: id_start = 1, id_end = 0, n_boundary = 0
 character(len=256) :: dir
 character(len=256) :: cwd
 real(8), allocatable :: av_lookup(:)

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only: isdead_or_accreted, iorig, rhoh, iboundary, iphase, iamtype, maxp
 use units,      only: utime,unit_density
 use physcon,    only: atomic_mass_unit
 use eos,        only: get_temperature, ieos, gamma,gmw, init_eos
 use io,         only: fatal, iverbose

 character(len=*), intent(in) :: dumpfile
 integer,          intent(in) :: num,npart,iunit
 real,             intent(in) :: xyzh(:,:),vxyzu(:,:)
 real,             intent(in) :: particlemass,time
 real, save    :: tprev = 0.
 integer, save :: nprev = 0
 real          :: dt_cgs, rho_cgs, numberdensity, T_gas, gammai, mui, AV
 integer       :: i, j, k, ierr, completed_iterations
 integer       :: iu=10,ios
 logical       :: iexist
 character(len=9) :: filename
 integer :: isize
 integer(4) :: npart_av, dump_num, last_under
 integer(4), allocatable :: av_id_tmp(:)
 real(8), allocatable :: av_val_tmp(:)
 character(len=40) :: avfile
 integer :: iu_av
 namelist /trace_config/ id_start, id_end, n_boundary

 if (.not.done_init) then
      allocate(one(maxp))
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

      ! read trace configuration from trace.cfg
      open(iu, file='trace.cfg', status='old', action='read', iostat=ios)
      if (ios /= 0) call fatal(analysistype, &
         "Could not open trace.cfg. Create it with: &trace_config id_start=1, id_end=10000, n_boundary=0 /")
      read(iu, nml=trace_config, iostat=ios)
      if (ios /= 0) call fatal(analysistype, "Failed to read trace.cfg namelist &trace_config")
      close(iu)
      if (id_end < id_start) call fatal(analysistype, "trace.cfg: id_end must be >= id_start")
      ntrack = id_end - id_start + 1
      print*, "Reading config from trace.cfg"
      print*, "  Tracking particle ID range: ", id_start, " to ", id_end, " (", ntrack, " particles)"
      print*, "  Skipping first ", n_boundary, " particles (boundary)"

      ! make directory to store individual particle trace files in current working directory
      dir = 'trace_output/'
      call getcwd(cwd)
      print *, "Creating trace_output directory in: ", trim(cwd)//'/'//trim(dir)
      inquire(file=dir, exist=iexist)
      if (.not.iexist) then
         call system('mkdir '//dir)
      end if

      allocate(track_id(ntrack))
      allocate(mask(maxp))
      allocate(av_lookup(maxp))
      av_lookup = 0.0d0
      ! fill track_id array from range
      do i=1,ntrack
         track_id(i) = id_start + i - 1
      end do 
      print*, "Tracking ", ntrack, " particles"
      call init_eos(ieos, ierr)
      if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")
      print*, "=== Initialisation complete: tracking ", ntrack, " particles (IDs ", id_start, "-", id_end, ") ==="
   
   else
      dt_cgs = (time - tprev)*utime
      completed_iterations = 0
      print*, dumpfile, ": not first step data, timestep = ",dt_cgs, "npart = ",npart, "nprev = ",nprev
      ! Read A_V binary file for this dump (written by column_densities.py)
      last_under = scan(dumpfile, '_', back=.true.)
      read(dumpfile(last_under+1:), *) dump_num
      write(avfile, '("AV/AV_",i5.5)') dump_num
      iu_av = 11
      open(iu_av, file=trim(avfile), access='stream', form='unformatted', status='old', iostat=ios)
      if (ios /= 0) then
         print*, "WARNING: Could not open A_V file "//trim(avfile)//", using A_V = 0"
         av_lookup = 0.0d0
      else
         read(iu_av) npart_av
         allocate(av_id_tmp(npart_av), av_val_tmp(npart_av))
         read(iu_av) av_id_tmp
         read(iu_av) av_val_tmp
         close(iu_av)
         av_lookup = 0.0d0
         do k = 1, npart_av
            if (av_id_tmp(k) >= 1 .and. av_id_tmp(k) <= maxp) then
               av_lookup(av_id_tmp(k)) = av_val_tmp(k)
            endif
         enddo
         deallocate(av_id_tmp, av_val_tmp)
         print*, "  Loaded A_V for ", npart_av, " particles from "//trim(avfile)
      endif
      ! update mask: track all particles whose original ID falls in [id_start, id_end]
      mask = .false.
      do j=1,npart
         if (iorig(j) >= id_start .and. iorig(j) <= id_end) then
            mask(j) = .true.
         endif
      enddo
      !$omp parallel do default(none) &
      !$omp shared(npart,xyzh,vxyzu,nprev,iorig,iorig_old,iprev,iverbose,dir) &
      !$omp shared(particlemass,unit_density,mask,time,utime,iphase,n_boundary,id_start,id_end) &
      !$omp shared(ieos,gamma,gmw,completed_iterations,av_lookup) &
      !$omp private(i,j,k,rho_cgs,numberdensity,T_gas,gammai,mui,AV,filename,iu,isize)
      outer: do i=1,npart
         iu = 10
         if (mask(i) .eqv. .true. .and. .not.isdead_or_accreted(xyzh(4,i))) then
            ! never write particles outside configured ID range.
            if (iorig(i) < id_start .or. iorig(i) > id_end) cycle outer
            inner: do j=1,nprev
               if (iorig(i) == iorig_old(j)) then
                  iprev(i) = j
                  exit inner
               endif
            enddo inner
            if (iamtype(iphase(i)) /= iboundary .and. i > n_boundary) then ! skip boundary particles
               !Thermodynamic quantities
               rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
               gammai = gamma
               mui    = gmw
               numberdensity = rho_cgs / (mui * atomic_mass_unit)
               T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gammai,mui)
               T_gas = max(T_gas,20.0d0) ! Floor temperature at 20K

               ! A_V from pre-computed binary file
               AV = real(av_lookup(iorig(i)))

            !$omp critical
            write(filename, '(i9)') iorig(i)
            inquire(file=trim(dir)//trim(adjustl(filename))//'.phys', size=isize)
            if (isize == -1) then
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.phys', status='new', action='write')
               print *, 'Creating new file for particle ', iorig(i)
               write(iu, *) '# time(s)   X(AU)   Y(AU)   Z(AU)   n(cm-3)   T(K)   A_V'
            else if (isize == 0) then
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.phys', status='old', action='write')
               print *, 'Filling empty file for particle ', iorig(i)
               write(iu, *) '# time(s)   X(AU)   Y(AU)   Z(AU)   n(cm-3)   T(K)   A_V'
            else
               open(iu, file=trim(dir)//trim(adjustl(filename))//'.phys', status='old', action='write', position='append')
            endif 
               ! write physical parameters to file
               write(iu, '(ES16.8,1x,ES14.7,1x,ES14.7,1x,ES14.7,1x,ES14.7,1x,F8.2,1x,F10.8,1x)')&
                time*utime, xyzh(1, i), xyzh(2, i), xyzh(3, i), numberdensity, T_gas, AV
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
      print*, "=== Analysis of ", trim(dumpfile), " complete ==="
   endif
   
 ! store current step data before moving on to next step
 nprev = npart
 iorig_old(1:npart) = iorig(1:npart)
end subroutine do_analysis

end module analysis