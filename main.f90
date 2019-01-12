Program Test_Psolver_BigDFT 
!#####################################################################!
!                                                                     !
! This program reads an external rho.dat density and cuts it          !
! into two pieces according to a XYZcut given in one axis.            ! 
! It can also nullify values up/down to a XYZzero value in a given    ! 
! axis.                                                               ! 
!                                                                     !
! Remember: if cut is done, two smaller boxes are generated, if nulli-!
! fying is active, two boxes of the same size are generated (with     !
! zeros at certain values.                                            !
!                                                                     !
!##### HOW IT WORKS:                                                  !
!                                                                     !
! You must define the grid of the entry data, ie: gridX, gridY, gridZ !
! You must deifine Xcut, Ycut and Zcut values (zeros mean no cut)     !
! OR you must define Xzero, Yzero and Zzero to nullify values up/down !
! that limit.                                                         !
! eg:                                                                 ! 
!   parameter(Xcut=0, Ycut=0, Zcut=0)                                 !
!   parameter(Xzero=34, Yzero=0, Zzero=0)                             !
!                                                                     !
! it will set grids to make no cut but nullify values obove/below 34. !
!                                                                     !
! You must change the tag argument of the copy calling, i.e.:         ! 
! call cpyarray(rho,grid,rho_cut,grid1,Zzero,'nzu')                      ! 
! to your needs. 'nzu' means 'n'ullify in 'z' axis 'u'p of a Zzero    ! 
! value. Other options for instance: 'cxd' cut X-axis down to a value.!
!                                                                     !
!##### INPUT/OUTPUT FILES:                                            ! 
!                                                                     !
! Input file name: rho.dat                                            !
! Output file names: cutup.dat and cutdn.dat                          !
!                                                                     !
!#####################################################################!
!#####################################################################!
! PLT@CFM(2018)                                                       !
!#####################################################################!
!
  use mpi, only: mpi_bcast
  use utils, only: readarray, wrtarray, cpyarray, set_grid
  use psolver_wrapper, only: eenergy_from_sum, psolver 
!
! Variables: 
!
  IMPLICIT NONE
!
  include "mpif.h"
!
  real*8, allocatable :: rho_cut(:,:,:,:), rho2(:,:,:,:)  ! partial densities 
  real*8, allocatable :: rho_cut_copy(:,:,:,:), rho2_copy(:,:,:,:)  ! auxiliary arrays
  real*8, allocatable :: rho(:,:,:,:)    ! total density
  real*8, allocatable :: rho_copy(:,:,:,:)      ! auxiliary arrray 
  real*8, allocatable :: pot(:,:,:,:)             ! total potential
  real*8, allocatable :: pot_cut(:,:,:,:), pot2(:,:,:,:) ! potentials of rho pieces
  real*8, allocatable :: potgthrd(:,:,:,:)     ! potential reconstructed by PSolver gathering both pieces
  real*8  ::  length, L                ! box length in angstrom, box lenght in a.u.
  integer, dimension(1,3) :: grid1     ! grid of rho_cut
  integer, dimension(1,3) :: grid2     ! grid of rho2
  integer, dimension(1,3) :: grid      ! total grid dimensions
  real*8, dimension(1,3)  :: hgrid     ! total grid separatations, supposed to be equal to all arrays 
  real*8  :: Uharrs, Uharrs1, Uharrs2  ! Hartree energy 
  real*8  :: Eelectros, Eelectros1, Eelectros2 ! Hartree energy from manual rho*V*dr solution 
  real*8  :: Egathered                 ! ... the same but for a gathered potential via PS_gather 
  integer :: grid1X, grid1Y, grid1Z    ! grid divisions of rho_cut
  integer :: grid2X, grid2Y, grid2Z    ! grid divisions of rho2
  integer :: gridX, gridY, gridZ       ! rho total grid divisions
  real*8  :: Volume, dV, rho_sum, pot_sum, convert
  integer :: nt1, nt2, nt              ! total number of grid points for each array
  integer :: Xcut, Ycut, Zcut          ! to make a cut in a given direction
  integer :: Xzero, Yzero, Zzero          ! to make a cut in a given direction
  integer :: i, j, k, w
  logical :: gdistrib=.true.           ! global distribution? how arrays are distributed for BigDFT solver
                                       ! if .true. all mpi_threads know all the array. 
  integer :: nspin, ispin              ! number of spin components, counter on spin components
!
! MPI version 
!
  integer ierr
  integer Nodes 
  integer node 
!
!
  parameter(length=10.0D0, convert=0.529177249D0)
  parameter(gridX=64, gridY=64, gridZ=64)
  parameter(Xcut=0, Ycut=0, Zcut=32)     ! For cutting densities. 
  parameter(Xzero=0, Yzero=0, Zzero=0)  ! For nullifying values.  
  parameter(nspin=1)                    ! Restricted or LSD calculation? 
!#####################################################################!
!
! Library initialisation: 
! 
  call MPI_Init ( ierr )
  call f_lib_initialize()
!
!  Determine this process's rank.
!
  call MPI_Comm_rank ( MPI_COMM_WORLD, node, ierr )
!
!  Find out the number of processes available.
!
  call MPI_Comm_size ( MPI_COMM_WORLD, Nodes, ierr )
!
!  Have Process 0 say hello.
!
  if ( node == 0 ) then
    call timestamp ( )
    write ( *, '(a)' ) ' Cut_and_gather_rho program: '
    write ( *, '(a)' ) '  A MPI test of Psolver usage'
    write ( *, '(a,i8,a)' ) '  Running on:', Nodes , ' processors'
  endif
! 
! We work on a.u.:
! 
  L=length/convert
  Volume= L**3
!  
! Grid points of total rho:
! 
   grid(1,1) = gridX
   grid(1,2) = gridY
   grid(1,3) = gridZ
!
   nt = grid(1,1)*grid(1,2)*grid(1,3)
!
! Array allocation
!
!
 allocate(rho(grid(1,1),grid(1,2),grid(1,3),nspin))
 allocate(pot(grid(1,1),grid(1,2),grid(1,3),nspin))
 allocate(potgthrd(grid(1,1),grid(1,2),grid(1,3),nspin))
!
! Set partial grids according to the cut defined by Xcut, Ycut and Zcut:
!
 call set_grid(grid1,grid2,grid,Xcut,Ycut,Zcut,nt1,nt2)
! 
! Allocation of partial arrays:
!
 allocate(rho_cut(grid1(1, 1), grid1(1, 2), grid1(1, 3), nspin), pot_cut(grid1(1, 1), grid1(1, 2), grid1(1, 3), nspin))
 nt1 = grid1(1, 1) * grid1(1, 2) * grid1(1, 3)
!
! total grid separation, the same for all arrays: 
!
  hgrid = L/grid
! 
! We read rho: 
! 
 if (node==0) call readarray(rho, 10, 'rho.dat')
 if (Nodes.gt.1) call MPI_BCAST(rho, nt, MPI_DOUBLE, 0, MPI_COMM_WORLD, ierr)
!
! At this point both processors carry rho total 
! 
!!! TOTAL RHO
! Electrostastatic potential calculation for rho (total density)
!
  gdistrib=.true.   ! because all processors know rho
  call psolver(rho,pot, hgrid, Uharrs, node, Nodes, gdistrib)
  call eenergy_from_sum(rho, pot, hgrid, Eelectros)  
  if (node==0) call wrtarray(pot, 40, 'inraw', 'V.dat')
!
!!! CREATING PARTIAL RHO
!
! Now we split in such a way that tasks know only a part of the density: 
!
 if (node==0) then
   call cpyarray(rho,rho_cut,Zcut,'czd') ! eg. nzu= 'n'ullify 'z'-values 'u'p from a Zcut
   call wrtarray(rho_cut,20,'inraw','rho0.dat')
 endif
!
 if (node==1) then
   call cpyarray(rho,rho_cut,Zcut,'czu') ! eg. nzu= 'n'ullify 'z'-values 'u'p from a Zcut
   call wrtarray(rho_cut,30,'inraw','rho2.dat')
 endif
! 
! 
! If you want to broadcast info between nodes: 
!
!  call MPI_BCAST(rho_cut, nt1, MPI_DOUBLE, 1, MPI_COMM_WORLD, ierr)
!  write(6,*) "sum(rho_cut)", node , sum(rho1)
! 
! Electrostastatic potential calculation of partial densities
!
   gdistrib=.false.   ! because only one processor knows rho_cut or rho2
   call psolver(rho_cut, pot_cut, hgrid, Uharrs1, node, nodes, gdistrib, potgthrd)
   call eenergy_from_sum(rho_cut, pot_cut, hgrid, Eelectros)  
   write(*,*) "Uharrs1", node, Uharrs1
   CALL MPI_ALLGATHER(pot_cut, grid1(1,1)*grid1(1,2)*grid1(1,3), MPI_DOUBLE, pot, & 
                      grid(1, 1) * grid(1, 2) * grid(1, 3), &
                      MPI_DOUBLE, MPI_COMM_WORLD, ierr)
!
! Final test, potential from the 'entire' density must be equal to the one reconstructed: 
!
   if (node == 0) write(*,*) "sum(potgthrd) must be equal to sum(pot)", node, sum(potgthrd), sum(pot)
   if (node == 0) call eenergy_from_sum(rho, potgthrd, hgrid, Eelectros)  
!
! Deallocation 
! 
  deallocate(rho, rho_cut, pot, pot_cut)
!
! File closing 
!
  call MPI_Finalize ( ierr )
  if (node==0) call timestamp ( )
  write(*,*) ""
!
!!#####################################################################!
End Program Test_Psolver_BigDFT
!----------------------------------------------------------------------!
