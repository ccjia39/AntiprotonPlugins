program cippres_runci
  use general
  implicit none
  BEGIN_DOC
! CIPPRES stands for Configuration Interaction Plugin for Photoionized and Resonant Electronic States.
! ORMAS-CI calculations are performed to compute bound and approximate continuum states. Hamiltonian and Dipole coupling matrix elements can then be computed between these states. Finally, Stieltjes imaging technique is applied to recover the correct continuum coupling matrix elements.
! cippres_runci performs ORMAS-CI calculations using lists of CSFs generated previously (with cippres_gencsf)
  END_DOC

  character(len=lenmax) :: finput
  integer :: ilen, jlen
  logical :: file_e

  integer :: i, j, ista
  double precision, dimension(:,:), allocatable :: densmat
  double precision, dimension(:,:), allocatable :: densmatvec
  double precision, dimension(:), allocatable :: densmatval

  double precision :: vNentrop,linEntrop,Entrop

! TODO Read lists of CSFs from EZFIO (see cippres_gencsf first)


! GENERAL
! TODO Compute dipole matrix elements between two different CI runs
! TODO Include Stieltjes in qp

  print*,'orthonorm'
  call orthonormalize_mos
  print*,'done'

  if(ifcsf==1) then

    do i = 1, n_ciruns_cippres
!      print*,'ncsfs',i, n_csf_cippres(i) 
      eigvalues_cippres(1:n_csf_cippres(i),i) += nuclear_repulsion
      print*,'CI eigval =', eigvalues_cippres(1:n_csf_cippres(i),i)
    enddo

   call ezfio_set_cippres_eigvalues_cippres(eigvalues_cippres)
   call ezfio_set_cippres_eigvectors_cippres(eigvectors_cippres)
   call ezfio_set_cippres_ifcsf(2)
 
   integer :: irun
   irun = 1
!  print*,'irun=',irun
   call from_csf_to_det(irun)

   allocate(densmat(mo_num,mo_num))
   allocate(densmatvec(mo_num,mo_num))
   allocate(densmatval(mo_num))
   do ista = 1, N_states
   print*,""
   print*,"DENSITY MATRIX", ista
   print*,""
   write(30,*)ista
   do i = 1, mo_num
    do j = 1, mo_num
!       print*,j,i,(one_e_dm_mo_beta(j,i,ista)+one_e_dm_mo_alpha(j,i,ista))*0.5
       densmat(j,i) = 0.5d0*(one_e_dm_mo_beta(j,i,ista)+one_e_dm_mo_alpha(j,i,ista))
     enddo
    enddo
    call lapack_diagd(densmatval,densmatvec,densmat,mo_num,mo_num) 
       vNentrop = 0d0
       linEntrop = 1d0
       Entrop = 0d0
       densmatval(:) = densmatval(:)+1d-15
     do i = 1, mo_num
       write(30,*)densmatval(i)
       vNentrop = vNentrop - abs(densmatval(i)) * log(abs(densmatval(i)))/log(2d0)
       linEntrop = linEntrop - abs(densmatval(i))**2
       Entrop = Entrop - ( abs(densmatval(i)) * log(abs(densmatval(i))) + (1d0-abs(densmatval(i))) * log( (1d0-abs(densmatval(i)) ) ) )
     enddo
    write(20,*)eigvalues_cippres(ista,1),vNentrop,linEntrop, Entrop
   enddo
   deallocate(densmat,densmatvec,densmatval)

  else

   print*, "Please run cippres_gencsf first or type qp set cippres ifcsf 1"
   print*, "Note that if you rerun cippres_runci, you should first delete the eigval/vec from EZFIO"

  endif

end program cippres_runci
