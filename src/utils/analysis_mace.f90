!--------------------------------------------------------------------------!
! The Phantom Smoothed Particle Hydrodynamics code, by Daniel Price et al. !
! Copyright (c) 2007-2024 The Authors (see AUTHORS)                        !
! See LICENCE file for usage and distribution conditions                   !
! http://phantomsph.github.io/                                             !
!--------------------------------------------------------------------------!
module analysis
!
! Analysis routine to run MACE on phantom dumps
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
 use ftorch, only : torch_model, torch_tensor, torch_kCPU, torch_delete, &
                    torch_tensor_from_array, torch_model_load, torch_model_forward
 use dvode_module
 use part,       only: maxp
 use raytracer,  only: get_all_tau
 implicit none
 character(len=20), parameter, public :: analysistype = 'mace'
 public :: do_analysis

 integer, parameter :: abs_size = 468
 integer, parameter :: real_size = 472 ! 4 physical parameters + 468 species, input tensor size
 integer :: latent_size = 16
 character(len=8) :: epoch = '14'
 character(len=256) :: encoder_file, decoder_file


 real, dimension(real_size), target :: in_data
 real, dimension(abs_size), target :: in_data2
 real, allocatable, target :: latent_abs(:)
 real, allocatable, target :: latent_abs_evolved(:) 
 real, dimension(abs_size), target :: out_data
 real, dimension(abs_size), target :: expected
 character(len=32), dimension(abs_size) :: species
 character(len=256000000) :: line
 integer :: pos, start, nvals
 real, allocatable :: values(:)
 real, allocatable    :: abundance(:,:), abundance_prev(:,:), one(:)
 character(len=16)    :: abundance_label(abs_size)
 integer(8), allocatable :: iorig_old(:)
 integer, allocatable :: iprev(:)
 logical :: done_init = .false.
 real :: AuvAv = 4.65, albedo = 0.5

 integer :: ntrack = 0
 integer, allocatable :: track_id(:)
 logical, allocatable :: mask(:)
 character(len=256) :: dir

 ! Minmaxes for scaling
 real, dimension(12) :: minmax
 real :: rho_min, rho_max
 real :: T_min, T_max
 real :: delta_min, delta_max
 real :: Av_min, Av_max
 real :: dt_max, dt_fract
 real :: n_min, n_max

 ! Set up Torch data structures
 ! The net, a vector of input tensors (in this case we only have one), and the output tensor
 type(torch_model) :: model, encoder, decoder
 type(torch_tensor), dimension(1) :: in_tensors
 type(torch_tensor), dimension(1) :: latent_tensors
 type(torch_tensor), dimension(1) :: latent_evolved_tensors
 type(torch_tensor), dimension(1) :: out_tensors

 ! ODE parameters
 real           :: atol, rtol
 real, allocatable :: atol_arr(:), rtol_arr(:)
 real, allocatable :: C(:)
 real, allocatable :: A(:, :)
 real, allocatable :: B(:, :, :)
 character(len=16) :: label

 ! ODE solver variables
 type(dvode_t) :: solver
 integer  :: mf, istate,iopt, itol, itask, neq
 integer, parameter  :: lrw = 4096, liw = 46
 real :: t, tout, dt
 integer, dimension(liw)  :: iwork
 real, dimension(lrw) :: rwork
 real, allocatable :: y(:)

 ! Flag for testing
 logical :: test_pass
 logical :: verbose

 ! Timing variables
 real :: startTime, stopTime

 ! Model name
 character(len=256) :: model_file = '20251009_092808'


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
 real, save    :: tprev = 0.
 integer, save :: nprev = 0
 real          :: dt_cgs, rho_cgs, numberdensity, T_gas, gammai, mui, AUV, xi &
                , rholist(npart), Tlist(npart), mulist(npart), Auvlist(npart), xilist(npart)
 real          :: abundance_part(abs_size), column_density(npart), xyzh_copy(4,npart)
 real          :: max_radius, radius
 integer       :: i, j, k, i_radius, ierr, completed_iterations, npart_copy = 0
 integer       :: iu=10,ios
 character(len=9) :: filename
 integer :: isize

 if (.not.done_init) then
    done_init = .true.
    print*, "Initialising MACE"
    ! read in model metadata
    call read_metadata(trim(model_file)//'/'//trim('meta.txt'), latent_size, dt_fract, epoch)
    call read_species_labels(trim(model_file)//'/'//trim('species.txt'), species, abs_size)
    ! Path to torchscript models
    encoder_file = trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_encoder.pt')
    decoder_file = trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_decoder.pt')
    ! load autoencoder
    print*, " - Loading Torch models"
    call torch_model_load(encoder, encoder_file, torch_kCPU)
    call torch_model_load(decoder, decoder_file, torch_kCPU)

    ! load minmax values
    call read_minmax(trim(model_file)//'/'//trim('minmax.txt'), 12, minmax)
    rho_min = log10(minmax(1))
    rho_max = log10(minmax(2))
    T_min = log10(minmax(3))
    T_max = log10(minmax(4))
    delta_min = log10(minmax(5))
    delta_max = log10(minmax(6))
    Av_min = log10(minmax(7))
    Av_max = log10(minmax(8))
    n_min = log10(minmax(9))
    n_max = log10(minmax(10))
    dt_max = minmax(11)

    ! Initialise abundance arrays
    print*, " - Allocating arrays"
    maxp = maxp / 10
    allocate(abundance(abs_size,maxp))
    abundance = 0.
    allocate(abundance_prev(abs_size,maxp))
    abundance_prev = 0.
    allocate(latent_abs(latent_size))
    latent_abs = 0.
    allocate(latent_abs_evolved(latent_size))
    latent_abs_evolved = 0.

    ! Initialise other arrays
    allocate(one(maxp))
    one = 1.
    allocate(iorig_old(maxp))
    iorig_old = 0
    allocate(iprev(maxp))
    iprev = 0

    ! Allocate ODE solver arrays
    allocate(atol_arr(latent_size))
    allocate(rtol_arr(latent_size))
    allocate(C(latent_size))
    allocate(A(latent_size, latent_size))
    allocate(B(latent_size, latent_size, latent_size))
    allocate(y(latent_size))

    ! Initialise ODE parameters
    print*, " - Initialising ODE parameters"
    call read_ode_params(model_file, epoch, atol, rtol)
    neq    = latent_size ! number of first order odes.
    itol   = 2 !or 2 according as atol (below) is a scalar or array.
    itask  = 1 !for normal computation of output values of y at t = tout.
    istate = 1 !integer flag (input and output).  set istate = 1.
    iopt   = 0 !to indicate no optional input used.
    mf = 22
    atol_arr = atol
    rtol_arr = rtol

    print*, " - Setting abundances"
    do i=1, npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          call chem_init(abundance_part,species)
          abundance(:,i) = abundance_part
       endif
    enddo
    call init_eos(ieos, ierr)
    if (ierr /= 0) call fatal(analysistype, "Failed to initialise EOS")
    print*, "MACE initialisation complete"
 else
    dt_cgs = (time - tprev)*utime
    completed_iterations = 0
    print*, " - Not first step data, timestep = ",dt_cgs, "npart = ",npart, "nprev = ",nprev
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

    outer: do i=1,npart
       if (.not.isdead_or_accreted(xyzh(4,i))) then
          inner: do j=1,nprev
             ! get previous index of particle
             if (iorig(i) == iorig_old(j)) then
                iprev(i) = j
                exit inner
             endif
          enddo inner
          ! Thermodynamic quantities
          rho_cgs = rhoh(xyzh(4,i),particlemass)*unit_density
          gammai = gamma
          mui    = gmw
          numberdensity = rho_cgs / (mui * atomic_mass_unit)
          T_gas = get_temperature(ieos,xyzh(1:3, i),rhoh(xyzh(4,i),particlemass),vxyzu(:,i),gammai,mui)
          T_gas = max(T_gas,20.0d0)
          ! Radiation quantities
          AUV = 4.65 * column_density(i) / 1.87e21
          xi = get_xi(AUV)
          if (j == iprev(i)) then
            ! if particle existed in previous dump, evolve abundances
             abundance_part(:) = abundance_prev(:,iprev(i))
             in_data(1:4) = [numberdensity, T_gas, xi, AUV]
             in_data(5:real_size) = abundance_part(:)
             print*, "Particle ", iorig(i), " found in previous step, evolving abundances"
             in_data2 = in_data(5:real_size)

             ! Scale input data
             call cpu_time(startTime)
             call log_and_scale (in_data(1), rho_min, rho_max) ! density
             call log_and_scale (in_data(2), T_min, T_max)     ! temperature
             call log_and_scale (in_data(3), delta_min, delta_max) ! delta
             call log_and_scale (in_data(4), Av_min, Av_max)   ! Av
             !dt = 5.52500167E+07 ! from 5001.chem
             dt = dt_cgs / dt_max * dt_fract           ! dt ! dt will need to be time of dump - time of previous dump (first dump is just initialisation)
             do j = 6, real_size
                call log_and_scale (in_data(j), n_min, n_max) ! abundances
             enddo

             ! Create Torch input/output tensors from the above arrays
             call torch_tensor_from_array(in_tensors(1), in_data, torch_kCPU)
             call torch_tensor_from_array(latent_tensors(1), latent_abs, torch_kCPU)
             call torch_tensor_from_array(latent_evolved_tensors(1), latent_abs_evolved, torch_kCPU)
             call torch_tensor_from_array(out_tensors(1), out_data, torch_kCPU)
             ! Run MACE model
             ! Start timing
             ! Encode
             call torch_model_forward(encoder, in_tensors, latent_tensors)
               ! Evolve latent space
             ! Initialise ODE solver variables
             y      = latent_abs ! initial value of the dependent variable.
             t      = 0.0d0 ! initial value of the independent variable.
             tout   = dt ! first point where output is desired (/= t).
             call solver%initialize(f=ode)

             ! Call the ODE solver to evolve latent_abs
             call solver%solve(neq,y,t,tout,itol,rtol_arr,atol_arr,itask,istate,&
                                 iopt,rwork,lrw,iwork,liw,mf)
             write(*,*) "After ODE solve, istate =", istate
             latent_abs_evolved = y

             ! Decode
             call torch_model_forward(decoder, latent_evolved_tensors, out_tensors) 
             call cpu_time(stopTime)

             ! Unscale abundances
             do j = 1, abs_size
                call unscale_and_unlog(out_data(j), n_min, n_max)
             end do
             !do j =1, abs_size
             !   write(*,*) j, trim(species(j)), in_data2(j), out_data(j)
             !end do
             write(*, '(A, F8.6)') 'Elapsed time, s : ',  (stopTime - startTime)
             ! Update abundances
             abundance(:,i) = out_data(:)
          else
          ! initialise abundances if no previous match found (i.e. new particle)
             print*, "Particle ", iorig(i), " not found in previous step, initializing abundances"
             call chem_init(abundance_part,species)
          endif 
       endif
    enddo outer
   endif

 ! Cleanup
 !call torch_delete(model)
 !call torch_delete(in_tensors)
 !call torch_delete(latent_tensors)
 !call torch_delete(latent_evolved_tensors)
 !call torch_delete(out_tensors)

 ! store current step data before moving on to next step
 nprev = npart
 tprev = time
 iorig_old = iorig
 abundance_prev = abundance
 print*, species(72),abundance(72,1)
 print*, species(72),abundance_prev(72,1)
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


subroutine write_chem(i,abundance_part,time,radius,numberdensity,T_gas,mui,AUV,xi,abundance_label,dir)
   integer(8), intent(in) :: i
   real, intent(in) :: abundance_part(:), time, radius, numberdensity, T_gas, mui, AUV, xi
   character(len=16), intent(in) :: abundance_label(:)
   character(len=*), intent(in) :: dir
   integer :: iu, isize, k
   character(len=9) :: filename

   write(filename, '(i9)') i
   inquire(file=trim(dir)//trim(adjustl(filename))//'.chem', size=isize)
   if (isize == -1) then
      open(iu, file=trim(dir)//trim(adjustl(filename))//'.chem', status='new', action='write')
      print *, 'Creating new file for particle ', i
      write(iu, *) '# time(s)   radius(cm)   n(cm-3)   T(K)   mu   A_UV   xi   ', (abundance_label(k), k=1,abs_size) 
   else if (isize == 0) then
      open(iu, file=trim(dir)//trim(adjustl(filename))//'.chem', status='old', action='write')
      print *, 'Filling empty file for particle ', i
      write(iu, *) '# time(s)   radius(cm)   n(cm-3)   T(K)   mu   A_UV   xi   ', (abundance_label(k), k=1,abs_size)
   else
      open(iu, file=trim(dir)//trim(adjustl(filename))//'.chem', status='old', action='write', position='append')
   endif 
      ! write physical parameters to file
      write(iu, '(ES16.8,1x,ES14.7,1x,ES14.7,1x,F8.2,1x,F6.3,1x,F8.3,1x,F8.3,1x)',advance="no") time, radius, numberdensity, T_gas, mui, AUV, xi
      ! write abundances to file
   do k=1,abs_size
      write(iu, '(ES14.7,1x)', advance="no") abundance_part(k)
   enddo
   write(iu,*) ! new line
   close(iu)

 close(iu)
 
end subroutine write_chem

subroutine ode(me, neq, t, y, ydot)
      ! Defines the ODE system dyi/dt = Ci + Aij*yj + Bijk*yi*yk
      class(dvode_t),intent(inout) :: me
      integer :: neq
      real :: t
      real :: y(neq)
      real :: ydot(neq)
      integer :: i, j, k

      do j = 1, neq
         ydot(j) = C(j)
         do i = 1, neq
               ydot(j) = ydot(j) + A(j, i) * y(i)
         end do
         do i = 1, neq
               do k = 1, neq
                  ydot(j) = ydot(j) + B(j, i, k) * y(i) * y(k)
               end do
         end do
      end do
   end subroutine ode

   subroutine read_ode_params(model_file, epoch, atol, rtol)
      ! Subroutine to read ODE parameters from files
      implicit none
      character(len=*), intent(in), optional :: model_file, epoch
      real, intent(out) :: atol, rtol
      integer :: i, j, k

      ! Load ODE parameters from file
      print*, " - Reading ODE parameters from " // trim(trim(adjustl(epoch)) //'_ODE_params.txt')
      open(unit=10, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_params.txt'), status='old')
      read(10, *) label, atol
      read(10, *) label, rtol
      close(10)
      write(*,*) "      > atol =", atol
      write(*,*) "      > rtol =", rtol

      ! Read C (vector)
      print*, " - Reading Coefficient C from " // trim(trim(adjustl(epoch)) //'_ODE_C.txt')
      open(unit=11, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_C.txt'), status='old')
      do i = 1, latent_size
         read(11, *) C(i)
      end do
      close(11)

      ! Read A (matrix)
      print*, " - Reading Coefficient A from " // trim(trim(adjustl(epoch)) //'_ODE_A.txt')
      open(unit=12, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_A.txt'), status='old')
      do i = 1, latent_size
         do j = 1, latent_size
               read(12, *) A(i, j)
         end do
      end do
      close(12)

      ! Read B (tensor)
      print*, " - Reading Coefficient B from " // trim(trim(adjustl(epoch)) //'_ODE_B.txt')
      open(unit=13, file=trim(model_file)//'/'//trim(trim(adjustl(epoch))//'_ODE_B.txt'), status='old')
      do i = 1, latent_size
         do j = 1, latent_size
               do k = 1, latent_size
                  read(13, *) B(i, j, k)
               end do
         end do
      end do
      close(13)
   end subroutine read_ode_params

   subroutine log_and_scale (x, xmin, xmax)
      ! scale parameters
      ! First clip to min
      ! log transform then normalise
      ! xmin and xmax are log10 of min and max values for parameter

      real, intent(inout) :: x
      real, intent(in) :: xmin, xmax
      ! clip to min
      if (x < 10.0d0**xmin) then
         x = 10.0d0**xmin
      end if
      x = log10(x)
      x = (x - xmin) / ABS(xmax - xmin)
   end subroutine log_and_scale


   subroutine unscale_and_unlog (x, xmin, xmax)
      ! unscale parameters 
      ! denormalise then exp10 transform
      ! xmin and xmax are log10 of min and max values for parameter
      real, intent(inout) :: x
      real, intent(in) :: xmin, xmax
      x = x * ABS(xmax - xmin) + xmin
      x = 10.0d0**x
   end subroutine unscale_and_unlog

   subroutine read_minmax(filename, filesize, minmax)
      ! Subroutine to read minmax values from a file
      implicit none
      character(len=*), intent(in) :: filename
      character(len=128) :: dummy
      integer, intent(in) :: filesize
      real, dimension(filesize), intent(out) :: minmax
      integer :: i, unit_num,ios

      write(*,*) " - Reading minmax values from ", filename
      unit_num = 20
      open(unit=unit_num, file=filename, status='old')
         do i = 1, filesize
            read(unit_num, *, iostat=ios) dummy, dummy, minmax(i)
            if (ios /= 0) then
               write(*,*) "Error reading minmax file."
               stop 1
            end if
         end do
      close(unit_num)
   end subroutine read_minmax

   subroutine read_species_labels(filename,labels, num_labels)
      ! Subroutine to read species species labels from a file
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(in) :: num_labels
      character(len=32), dimension(num_labels), intent(out) :: labels
      integer :: i, unit_num, ios
      write(*,*) " - Reading ", num_labels, "species labels from ", filename
      unit_num = 30
      open(unit=unit_num, file=filename, status='old')
      do i = 1, num_labels
         read(unit_num, *, iostat=ios) labels(i)
         if (ios /= 0) then
            write(*,*) "Error reading species labels file."
            stop 1
         end if
      end do
      close(unit_num)
   end subroutine read_species_labels

   subroutine read_metadata(filename, latent_size, dt_fract, epoch)
      ! Subroutine to read metadata from a file
      implicit none
      character(len=*), intent(in) :: filename
      integer, intent(inout) :: latent_size
      character(len=8), intent(inout) :: epoch
      real, intent(inout) :: dt_fract
      character(len=128) :: label, metadata
      integer :: i, unit_num, ios, metadata_size
      write(*,*) " - Reading metadata from ", filename
      unit_num = 40
      open(unit=unit_num, file=filename, status='old')
      ! get length of metadata
      metadata_size = 0
      do
         read(unit_num, *, iostat=ios)
         if (ios /= 0) exit
         metadata_size = metadata_size + 1
      end do
      rewind(unit_num)
      ! read metadata
      do i = 1, metadata_size
         read(unit_num, *, iostat=ios) label, metadata
         if (trim(label) == 'z_dim:') then
            read(metadata, *) latent_size
            write(*,*) "      > latent_size =", latent_size
         else if (trim(label) == 'dt_fract:') then
            read(metadata, *) dt_fract
            write(*,*) "      > dt_fract =", dt_fract
         else if (trim(label) == 'epoch:') then
            read (metadata, *) epoch
            write(*,*) "      > epoch =", epoch
         end if

         if (ios /= 0) then
            write(*,*) "Error reading metadata file."
            stop 1
         end if
      end do
      close(unit_num)
   end subroutine read_metadata

   subroutine chem_init(abundance_part, species)
      ! Subroutine to initialise chemical abundances
      implicit none
      character(len=32), dimension(:), intent(in) :: species
      real, dimension(size(species)+4), intent(out) :: abundance_part
      integer :: i
      ! Initial abundances for the krome model taken from Agúndez et al. (2020)
      ! H2, He, CO, C2H2, HCN, N2, SiC2, CS, SiS, SiO, CH4, H2O, HCl, C2H4, NH3, HCP, HF, H2S
      do i = 1, abs_size
         if (species(i) == 'H2') then
            abundance_part(i) = 0.5d0
         else if (species(i) == 'He') then
            abundance_part(i) = 8.5d-2
         else if (species(i) == 'CO') then
            abundance_part(i) = 4.0d-4
         else if (species(i) == 'C2H2') then
            abundance_part(i) = 2.19d-5
         else if (species(i) == 'HCN') then
            abundance_part(i) = 2.045d-5
         else if (species(i) == 'N2') then
            abundance_part(i) = 2.0d-5
         else if (species(i) == 'SiC2') then
            abundance_part(i) = 9.35d-6
         else if (species(i) == 'CS') then
            abundance_part(i) = 5.3d-6
         else if (species(i) == 'SiS') then
            abundance_part(i) = 2.99d-6
         else if (species(i) == 'SiO') then
            abundance_part(i) = 2.51d-6
         else if (species(i) == 'CH4') then
            abundance_part(i) = 1.75d-6
         else if (species(i) == 'H2O') then
            abundance_part(i) = 1.275d-6
         else if (species(i) == 'HCl') then
            abundance_part(i) = 1.625d-7
         else if (species(i) == 'C2H4') then
            abundance_part(i) = 3.425d-8
         else if (species(i) == 'NH3') then
            abundance_part(i) = 3.0d-8
         else if (species(i) == 'HCP') then
            abundance_part(i) = 1.25d-8
         else if (species(i) == 'HF') then
            abundance_part(i) = 8.5d-9
         else if (species(i) == 'H2S') then
            abundance_part(i) = 2.0d-9
         else if (species(i) == 'E') then
            abundance_part(i) = 2.0d-11
         else
            abundance_part(i) = 1.0d-20
         end if
      end do
   end subroutine chem_init

end module analysis