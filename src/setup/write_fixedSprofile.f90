program write_fixed_S_profile
 use fixedSprofile, only:calc_rho_and_pres
 use readwrite_mesa,only:write_mesa
 use eos,           only:calc_temp_and_ene,ieos
 use table_utils,   only:flip_array
 implicit none
 real, allocatable, dimension(:) :: m,r,rho,pres,ene,temp
 character(len=120)              :: outputpath
 integer                         :: ierr,i

 call calc_rho_and_pres(m,r,rho,pres)
 allocate(ene(0:size(m)-1), temp(0:size(m)-1))
 print*,'Calculating thermal profile'
 do i = 0,size(m)-1
     call calc_temp_and_ene(ieos,rho(i),pres(i),ene(i),temp(i),ierr)
 enddo

 print*,'Flipping arrays to sort from surface to centre'
 call flip_array(m)
 call flip_array(pres)
 call flip_array(temp)
 call flip_array(r)
 call flip_array(rho)
 call flip_array(ene)

 outputpath = 'fixedSprofile.dat'
 call write_mesa(outputpath,m,pres,temp,r,rho,ene)
 print*,'Profile written to ',outputpath

 print*,'May this envelope stay put.'

end program write_fixed_S_profile