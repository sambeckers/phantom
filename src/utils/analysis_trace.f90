!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to generate data for chemistry
! Tracks particle paths and stores evolution in shared per-batch binary files
!
! :References: None
!
! :Owner: Camille Landri
!
! :Runtime parameters: None
!
! :Dependencies: None
!

 use iso_fortran_env, only: int64, real64
 implicit none
 character(len=20), parameter, public :: analysistype = 'trace'
 public :: do_analysis

 logical :: done_init = .false.
 integer :: ntrack = 0
 integer :: id_start = 1, id_end = 0
 integer :: n_batches = 0
 integer, allocatable :: batch_for_id(:)
 integer, allocatable :: batch_units(:)
 character(len=256) :: batch_map_file = 'trace_batch_map.txt'
 character(len=256) :: dir
 character(len=256) :: cwd
 real(8), allocatable :: av_lookup(:)

contains

subroutine do_analysis(dumpfile,num,xyzh,vxyzu,particlemass,npart,time,iunit)
 use part,       only: isdead_or_accreted, iorig, rhoh, iboundary, iphase, iamtype
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
 integer       :: i, k, ierr, completed_iterations
 integer       :: iu=10,ios
 logical       :: iexist
 logical, allocatable :: mask(:)
 integer :: batch_idx
 integer(int64) :: pid_out
 real(real64) :: time_out, x_out, y_out, z_out, density_out, temp_out, av_out
 integer(4) :: npart_av, dump_num, last_under
 integer(4), allocatable :: av_id_tmp(:)
 real(8), allocatable :: av_val_tmp(:)
 character(len=40) :: avfile
 character(len=256) :: batch_file
 integer :: iu_av
 namelist /trace_config/ id_start, id_end, n_batches, batch_map_file

 if (.not.done_init) then
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
         "Could not open trace.cfg. Create it with: &trace_config id_start=1, id_end=10000, n_batches=1 /")
      read(iu, nml=trace_config, iostat=ios)
      if (ios /= 0) call fatal(analysistype, "Failed to read trace.cfg namelist &trace_config")
      close(iu)
      if (id_end < id_start) call fatal(analysistype, "trace.cfg: id_end must be >= id_start")
      if (n_batches <= 0) call fatal(analysistype, "trace.cfg: n_batches must be > 0")
      print*, "Reading config from trace.cfg"
      print*, "  Tracking particle ID range: ", id_start, " to ", id_end
      print*, "  Batch map file: ", trim(batch_map_file)
      print*, "  Number of batches: ", n_batches

      ! make directory to store batch trace binary files in current working directory
      dir = 'trace_output/'
      call getcwd(cwd)
      print *, "Creating trace_output directory in: ", trim(cwd)//'/'//trim(dir)
      inquire(file=dir, exist=iexist)
      if (.not.iexist) then
         call system('mkdir '//dir)
      end if

      ! Allocate arrays sized to id_end: these are indexed by original particle ID
      ! (iorig), which is stable for the whole run and bounded by id_end from trace.cfg.
      ! Using id_end instead of maxp avoids silent out-of-bounds access if Phantom
      ! internally reallocates its particle arrays (growing maxp) mid-run.
      allocate(av_lookup(id_end))
      allocate(batch_for_id(id_end))
      allocate(batch_units(n_batches))
      av_lookup = 0.0d0
      batch_for_id = -1

      open(iu, file=trim(batch_map_file), status='old', action='read', iostat=ios)
      if (ios /= 0) call fatal(analysistype, "Could not open batch map file listed in trace.cfg")
      ntrack = 0
      do
         read(iu, *, iostat=ios) pid_out, batch_idx
         if (ios /= 0) exit
         if (pid_out >= int(id_start, int64) .and. pid_out <= int(id_end, int64)) then
            if (batch_idx >= 0 .and. batch_idx < n_batches) then
               batch_for_id(int(pid_out)) = batch_idx
               ntrack = ntrack + 1
            endif
         endif
      enddo
      close(iu)
      if (ntrack <= 0) call fatal(analysistype, "Batch map contained zero valid tracked particles")

      do i=1,n_batches
         write(batch_file, '(a,"trace_batch_",i5.5,".bin")') trim(dir), i-1
         open(50+i, file=trim(batch_file), access='stream', form='unformatted', status='replace', action='write', iostat=ios)
         if (ios /= 0) call fatal(analysistype, "Could not create batch trace binary file")
         batch_units(i) = 50 + i
      enddo

      print*, "Tracking ", ntrack, " particles from batch map"
      call init_eos(ieos, ierr)
      if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")
      print*, "=== Initialisation complete: tracking ", ntrack, " particles across ", n_batches, " batch files ==="
      done_init = .true.

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
            ! Guard against IDs outside the tracked range before indexing av_lookup
            if (av_id_tmp(k) >= id_start .and. av_id_tmp(k) <= id_end) then
               av_lookup(av_id_tmp(k)) = av_val_tmp(k)
            endif
         enddo
         deallocate(av_id_tmp, av_val_tmp)
         print*, "  Loaded A_V for ", npart_av, " particles from "//trim(avfile)
      endif

      ! mask is indexed by current slot j (1..npart), not by original ID,
      ! so it must be sized to npart for this dump and allocated fresh each call.
      allocate(mask(npart))
      mask = .false.
      do i=1,npart
         if (iorig(i) >= id_start .and. iorig(i) <= id_end) then
            if (batch_for_id(iorig(i)) >= 0) mask(i) = .true.
         endif
      enddo

      !$omp parallel do default(none) &
      !$omp shared(npart,xyzh,vxyzu,iorig,iverbose,dir) &
      !$omp shared(particlemass,unit_density,mask,time,utime,iphase,id_start,id_end) &
      !$omp shared(ieos,gamma,gmw,completed_iterations,av_lookup,batch_for_id,batch_units,n_batches) &
      !$omp private(i,k,rho_cgs,numberdensity,T_gas,gammai,mui,AV,iu,batch_idx) &
      !$omp private(pid_out,time_out,x_out,y_out,z_out,density_out,temp_out,av_out)
      do i=1,npart
         iu = 10
         if (mask(i) .and. .not.isdead_or_accreted(xyzh(4,i))) then
            if (iamtype(iphase(i)) /= iboundary) then ! skip boundary particles by particle type
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
            batch_idx = batch_for_id(iorig(i))
            if (batch_idx >= 0 .and. batch_idx < n_batches) then
               pid_out      = int(iorig(i), int64)
               time_out     = real(time*utime, real64)
               x_out        = real(xyzh(1, i), real64)
               y_out        = real(xyzh(2, i), real64)
               z_out        = real(xyzh(3, i), real64)
               density_out  = real(numberdensity, real64)
               temp_out     = real(T_gas, real64)
               av_out       = real(AV, real64)
               write(batch_units(batch_idx+1)) pid_out, time_out, x_out, y_out, z_out, density_out, temp_out, av_out
            endif
            !$omp end critical
            endif
         endif

         if (iverbose > 1) then
            !$omp atomic
            completed_iterations = completed_iterations + 1
            print*, 'Completed ', completed_iterations, ' of ', npart
         endif
      enddo
      !$omp end parallel do

      deallocate(mask)
      print*, "=== Analysis of ", trim(dumpfile), " complete ==="
   endif

 ! store current step data before moving on to next step
 nprev = npart
 tprev = time
end subroutine do_analysis

end module analysis