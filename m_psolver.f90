module psolver_wrapper  
!
!######################################################################!
! Module taha gahters all utilities to copy and write data.            !
!######################################################################!
! PLT@CFM(2018)                                                        !
!######################################################################!
  implicit none 
  private
  public :: eenergy_from_sum, psolver
  contains 

!!#####################################################################!
!----------------------------------------------------------------------!
  subroutine eenergy_from_sum(rho,pot,hgrid,energy)
!----------------------------------------------------------------------!
! This subroutine calculates the electrostatic energy between a given  ! 
! electronic density, rho, and a potential, pot, using a direct sum.   !
!----------------------------------------------------------------------!
  real(8), dimension(:,:,:,:), intent(in):: rho, pot
  real(8),intent(in) :: hgrid(1,3)
  real(8),intent(out):: energy 
  real(8)            :: dV
  integer :: i, j, k 
  integer :: ispin, nspin
  integer :: grid(1,3)
!----------------------------------------------------------------------!
! 
! Grid definintion: 
!
  nspin=size(rho,4)
  do i=1,3  
   grid(1,i)=size(rho,i)
  enddo
!
  energy=0.0D0
  dV=hgrid(1,1)*hgrid(1,2)*hgrid(1,3)
! 
  do ispin=1, nspin
    do i=1,grid(1,1) 
      do j=1,grid(1,2) 
        do k=1,grid(1,3) 
          energy=energy+rho(i,j,k,ispin)*pot(i,j,k,ispin)*dV
        enddo 
      enddo 
    enddo 
  enddo 
  energy=0.5D0*energy
  write(*,*) "Eelectrostatic E reconstructed by rho*V:", energy 
  write(*,*) "Using a grid and dV of:", grid, dV 
!----------------------------------------------------------------------!
  return 
  end subroutine eenergy_from_sum 
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
  subroutine psolver(rho, pot, hgrid, Eelectros, iproc, nproc,   &
                     gdistrib, potgthrd)
!----------------------------------------------------------------------!
! Calculation of electrostatic potential (V_Ne) Poisson solver of      !
! BigDFT package.                                                      !
! Arguments:                                                           !
! rho(in) = density on which calculation is going to be performed.     !
!            It can be either a total density or just a 'piece'.       ! 
! Pot(OUT) = calculated potential on rho1.                             !
! grid1(IN) = grid mesh division of rho1.                              !
! hgrid(IN) = grid mesh separation of rho1.                            !
! Eelectros(OUT) = electrostatic energy.                               !
! iprock(IN) = process rank.                                           !
! nprock(IN) = total number of MPI processes.                          !
! gdistrib(IN) = logical argument to indicate whether information of   !
!               rho1 is known by all processors or not.                !
! rho(IN) = bigger array used to collect the potential in case of      !
!           rho1 is only known by each process (ie. gdistrib=.false.). !
! grid(IN) = grid mesh division of rho.                                !
!----------------------------------------------------------------------!
! We use the psolver program of BigDFT package:                        !
  use poisson_Solver,  only: electrostatic_solver,            &
                              coulomb_operator, pkernel_set
  use pstypes,         only: pkernel_init, pkernel_free
  use psbox,           only: ps_gather, ps_gather_asimmetric
  use yaml_output
  use utils,           only: wrtarray 
  use futile
  use dictionaries, dict_set => set
!----------------------------------------------------------------------!
  real*8, intent(in)   :: rho(:,:,:,:)   ! Electronic density.
  real*8, allocatable  :: rho_tmp(:,:,:) ! Working array carrying only one spin component at a time.  
  real*8, allocatable  :: rhogthrd(:,:,:,:) ! Working array to carry the full electronic density.  
  real*8, intent(out)  :: pot(:,:,:,:)   ! Output potential.
  integer, intent(in)  :: iproc, nproc   ! ID process, number of processes.
  real*8, intent(out), optional :: potgthrd(:,:,:,:)  ! if gdistrib=.false., it is used to gather total rho 
  real*8, intent(in)  :: hgrid(1,3)      ! Grid separation.
  integer  :: grid(1,3)                  ! Grid mesh for rho and pot 
  integer  :: gridgthrd(1,3)                   ! grid mesh for rho and pot 
  real*8, intent(out) :: Eelectros       ! output electrostatic energy
  logical, intent(in) :: gdistrib        ! type of distribution.
  type(coulomb_operator) :: pkernel      ! BigDFT main object.
  type(dictionary), pointer :: dict      ! dictionary containing input parameters
  integer :: isf_order, fd_order         ! BigDFT input parameters 
  real*8  :: alpha,beta,gamma            ! possible angles of the box axis
  integer :: i, j, k 
  integer :: ispin, nspin
! -------------------------------------------------------------------- C
!
!!! Pring header:
!
  if (iproc==0) then
    write(*,*) ""
    write(*,*) "#---------------- PSOLVER CALCULATION -------------------#"
  endif
!
!!! Define working array grids used in the kernel definition: 
!
  nspin = size(rho, 4)
!
! global grid for rho and pot:
!
  do i = 1, 3  
   grid(1,i) = size(rho, i)
  enddo
  write(*,*) "Grid of working array", grid 
!
! if gathered array is define:
!
  if (present(potgthrd)) then
    do i = 1, 3
      gridgthrd(1, i) = size(potgthrd, i)
    enddo
    write(*,*) "Grid of array for gathering:", gridgthrd
  endif
! 
  write(*,*) "Just before Psolver sum(rho)", iproc, sum(rho)
  if (present(potgthrd)) write(*,*) "Just before Psolver sum(potgthrd)", iproc, sum(potgthrd)
!
!!! Now it comes the proper calculation. Two cases are defined depen-  !
!   ding on whether all processors know all the information or not.    !
!
  if (gdistrib) then    ! global distribution, everyone knows everything
!
    allocate(rho_tmp(grid(1,1),grid(1,2),grid(1,3)))
!
! Dictionary set-up 
!
    isf_order=16 !  2, 4, 6, 8, 14, 16, 20, 24, 30, 40, 50, 60, 100
    call dict_init(dict)
    call set(dict//'setup'//'global_data',gdistrib)
!!  call set(dict//'setup'//'taskgroup_size',2)
!!  dict=> dict_new('kernel' .is. dict_new('isf_order' .is. isf_order))
!
    do ispin=1, nspin      ! loop on spin component 
!
! Kernel initialisation:
!
       pkernel=pkernel_init(iproc, nproc, dict, 'F', grid, hgrid)
       call dict_free(dict)
       call pkernel_set(pkernel)
       rho_tmp=rho(:,:,:,ispin)           ! Potential enters containing rho at first
!
! Main calling to solver: 
!
       write(*,*) "Entering Electrostatic_Solver:"
!
       CALL Electrostatic_Solver(pkernel,rho_tmp,ehartree=Eelectros)
       pot(:,:,:,ispin)=rho_tmp
!
       write(*,*) "After Psolver sum(Pot)", iproc, sum(pot), "for spin component:", ispin
       write(*,*) "After Psolver Eelectros", iproc, Eelectros, "for spin component:", ispin 
!!
      CALL pkernel_free(pkernel)
      deallocate(rho_tmp)
!
    enddo                   ! end loop on spin component
!
  else                  ! the entire rho is not known by all processors. 
!
    allocate(rho_tmp(grid(1,1),grid(1,2),grid(1,3)))
    allocate(rhogthrd(gridgthrd(1,1),gridgthrd(1,2),gridgthrd(1,3),nspin))
!
! Dictionary set-up 
!
    isf_order=16 !  2, 4, 6, 8, 14, 16, 20, 24, 30, 40, 50, 60, 100
    call dict_init(dict)
    call set(dict//'setup'//'global_data', .true.)
!
    do ispin=1, nspin      ! loop on spin component 
!
      write(*,*) "Distributed array information:"
      write(*,*) "grid for process", iproc,  grid
      write(*,*) "size(Pot1)", size(rho) 
!
      write(*,*) "Entering Electrostatic_Solver:"
      rho_tmp=rho(:,:,:,ispin)
      pkernel=pkernel_init(iproc, nproc, dict, 'F', gridgthrd, hgrid)
      call dict_free(dict)
      call pkernel_set(pkernel)
!      call PS_gather(rho_tmp, pkernel, rhogthrd, ispin)
      call PS_gather_asimmetric(rho_tmp, pkernel, rhogthrd, grid) 
      if (iproc == 0) call wrtarray(rhogthrd, 70,'short','rhogathered.dat')! 
      if (iproc == 0) write(*,*) "node", iproc, "says sum(rhogthrd) before Psolver", sum(rhogthrd(:,:,:,ispin))
!
       write(*,*) "Electrostatic potential in processor", iproc, "out of:", nproc, "will be collected."  
       rho_tmp=rhogthrd(:,:,:,ispin)   ! rho_tmp enters the density per spin component
       CALL Electrostatic_Solver(pkernel,rho_tmp,ehartree=Eelectros)
       potgthrd(:,:,:,ispin)=rho_tmp(:,:,:)
!
      if (iproc == 0) then
        write(*,*) "After Psolver sum(rhogthrd)", iproc, sum(rhogthrd(:,:,:,ispin)), "spin component:", ispin 
        write(*,*) "After Psolver sum(potgthrd)", iproc, sum(potgthrd(:,:,:,ispin)), "spin component:", ispin  
        write(*,*) "After Psolver Eelectros", iproc, Eelectros, "spin component:", ispin  
      endif
      CALL pkernel_free(pkernel)
    enddo                   ! end loop on spin component
      deallocate(rho_tmp)
!
  endif   ! end of  main condition
! 
  if (iproc==nproc)  then
   write(*,*) "#--------------- END PSOLVER CALCULATION -----------------#"
   write(*,*) ""
  endif
!
!----------------------------------------------------------------------!
 return 
 end subroutine psolver 
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
end module  psolver_wrapper
