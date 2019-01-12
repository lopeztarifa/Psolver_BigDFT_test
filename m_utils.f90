module utils 
!
!######################################################################!
! Module taha gahters all utilities to copy and write data.            !
!######################################################################!
! PLT@CFM(2018)                                                        !
!######################################################################!
  implicit none 
  private
  public :: readarray, wrtarray, cpyarray, set_grid
  CONTAINS

! -------------------------------------------------------------------- !
!     subroutineS #####################################################!
! -------------------------------------------------------------------- !
      subroutine set_grid(grid1,grid2,grid,Xcut,Ycut,Zcut,nt1,nt2)
! -------------------------------------------------------------------- !
! Given a total grid (grid), the subroutine sets a two partial grids   ! 
! (grid1,grid2) defined by the Xcut, Ycut and Zcut values.             !
! -------------------------------------------------------------------- !
      integer:: I,J,K
      integer,intent(OUT):: grid1(1,3), grid2(1,3)
      integer,intent(IN) :: grid(1,3), Xcut, Ycut, Zcut
      integer,intent(OUT):: nt1, nt2  ! total number of grid points  
! -------------------------------------------------------------------- !
   if (Xcut.ne.0) then
    write(*,*) "A cut on X-axis will be made. Xcut:", Xcut
    grid1(1,1)=Xcut
    grid1(1,2)=grid(1,2)
    grid1(1,3)=grid(1,3)
!
    grid2(1,1)=grid(1,1)-Xcut
    grid2(1,2)=grid(1,2)
    grid2(1,3)=grid(1,3)
   elseif (Ycut.ne.0) then
    write(*,*) "A cut on Y-axis will be made. Ycut:", Ycut
    grid1(1,1)=grid(1,1)
    grid1(1,2)=Ycut
    grid1(1,3)=grid(1,3)
!
    grid2(1,1)=grid(1,1)
    grid2(1,2)=grid(1,2)-Ycut
    grid2(1,3)=grid(1,3)
   elseif (Zcut.ne.0) then
    write(*,*) "A cut on Z-axis will be made. Zcut:", Zcut
    grid1(1,1)=grid(1,1)
    grid1(1,2)=grid(1,2)
    grid1(1,3)=Zcut
!
    grid2(1,1)=grid(1,1)
    grid2(1,2)=grid(1,2)
    grid2(1,3)=grid(1,3)-Zcut
   else
    grid1=grid
    grid2=grid
   endif
   nt1=grid1(1,1)*grid1(1,2)*grid1(1,3)
   nt2=grid2(1,1)*grid2(1,2)*grid2(1,3)
      
! -------------------------------------------------------------------- !
  end subroutine set_grid
! -------------------------------------------------------------------- !
! -------------------------------------------------------------------- !
  subroutine readarray(m,unitnumber,fname)
! -------------------------------------------------------------------- !
! This subroutines reads a 1-column ASCII array, into an array.        !
! -------------------------------------------------------------------- !
  real(8), dimension(:,:,:,:),intent(out):: M
  character(len=*), intent(in) :: fname
  integer,intent(in):: unitnumber 
  integer :: i, j, k 
  integer :: nspin, ispin
  integer :: grid(1,3)
! -------------------------------------------------------------------- !
! 
!!! Assigning dimension:
!
  open(unitnumber, file=trim(fname))
  nspin=size(m,4)
  do i=1,3
   grid(1,i)=size(m,i)
  enddo
!
  do ispin=1,nspin
    do i = 1, grid(1,1)
      do j = 1, grid(1,2)
        do k = 1, grid(1,3)
          read(unitnumber,*) m(i,j,k,ispin)
        enddo
      enddo
    enddo
  enddo
!
  close(unitnumber)
! -------------------------------------------------------------------- !
  end subroutine readarray
! -------------------------------------------------------------------- !
! -------------------------------------------------------------------- !
  subroutine wrtarray(m,unitnumber,tag,fname)
! -------------------------------------------------------------------- !
! This subroutine write a M matrix.                                    !
! -------------------------------------------------------------------- C
  real(8), dimension(:,:,:,:),intent(in):: M
  character(len=*), intent(in) :: fname
  integer,intent(in):: unitnumber 
  character(5),intent(in):: tag
  integer :: i, j, k 
  integer :: nspin, ispin
  integer :: grid(1,3)
! -------------------------------------------------------------------- C
!
  open(unitnumber, file=trim(fname))
  nspin=size(m,4)
  do i=1,3
   grid(1,i)=size(m,i)
  enddo
!  
  write(*,*) "Writing array"
!
  if (tag == 'short') then  
    do ispin = 1, nspin 
      do i = 1, grid(1,1)
        do j = 1, grid(1,2)
          do k = 1, grid(1,3)
            write(unitnumber,10) M(I,J,K,ispin)
          enddo
        enddo
      enddo
    enddo
  endif
!
  if (tag == 'large') then  
    do ispin = 1, nspin 
      do i = 1, grid(1,1)
        do j = 1, grid(1,2)
          do k = 1, grid(1,3)
            write(unitnumber,20) M(I,J,K,ispin)
          enddo
        enddo
      enddo
    enddo
  endif
!
  if (tag == 'inraw') then  
    do ispin = 1, nspin 
      do i = 1, grid(1,1)
        do j = 1, grid(1,2)
          do k = 1, grid(1,3)
            write(unitnumber,*) M(I,J,K,ispin)
          enddo
        enddo
      enddo
    enddo
  endif
!
!
  close(unitnumber)
  10 format (100F10.4)
  20 format (100F10.6)
! -------------------------------------------------------------------- C
  end subroutine wrtarray 
!----------------------------------------------------------------------!
!----------------------------------------------------------------------!
!!!#####################################################################!
!! -------------------------------------------------------------------- C
!      subroutine REORDERARRAY(M, dimM, tag)
!      IMPLICIT NONE
!! -------------------------------------------------------------------- C
!! This subroutine reorders an array according to tag.                  ! 
!! -------------------------------------------------------------------- C
!      integer:: I, J, K, W
!      integer,intent(IN):: dimM(1,3)
!      real(8), dimension(dimM(1,1),dimM(1,2),dimM(1,3)),intent(INOUT) :: M
!      real(8), dimension(dimM(1,1),dimM(1,2),dimM(1,3)) :: N
!      character(3),intent(IN):: tag
!! -------------------------------------------------------------------- C
!      N=M
!      if (tag=='zxy') then
!      do I=1, DIM(1,1)
!       do J=1, DIM(1,2)
!        do K=1, DIM(1,3)
!        M(I,J,K)=N(I,J,K)
!        enddo
!       enddo
!      enddo
!      endif
!
!! -------------------------------------------------------------------- C
!      end subroutine REORDERARRAY 
!!#####################################################################!
! -------------------------------------------------------------------- C
      subroutine cpyarray(M, N, cut, tag)
! -------------------------------------------------------------------- C
! Subroutine that copies partial values of M into N according to cut   ! 
! and tag.                                                             ! 
! -------------------------------------------------------------------- C
      real(8), dimension(:, :, :, :),intent(IN) :: M
      real(8), dimension(:, :, :, :),intent(OUT):: N
      character(3), intent(in):: tag
      integer, intent(IN):: cut 
      integer :: gridM(1, 3), gridN(1, 3)
      integer :: I, J, K, W
      integer :: nspin, ispin 
      logical :: consistency  = .false.
! -------------------------------------------------------------------- C
!
! Grid set-up:
!
  nspin=size(m,4)
  do i=1,3
   gridM(1,i)=size(M,i)
   gridN(1,i)=size(N,i)
  enddo
!
! Consistency check:
!
  if ((gridM(1,1).eq.gridN(1,1)).and.(gridM(1,2).eq.gridN(1,2)).and.(gridM(1,3).eq.gridN(1,3))) &
  consistency=.true.
!
  write(*,*) "Consistency of arrays:", consistency
!
! Cuting ----------------------------!
!
! cut on Xaxis 'cx'
! down to a Xcut 'cxd'
  if (tag == 'cxd') then
    write(*,*) "Cuting values of rho.dat on X-axis down to a cut:", cut
    do ispin = 1, nspin
      W=1
      do i = 1, cut
        do j = 1, gridM(1, 2)
          do k= 1, gridM(1, 3)
            N(w, j, k, ispin)=M(i, j, k, ispin)
          enddo
        enddo
        W=W+1
      enddo
    enddo
  endif
! up to a Xcut 'cxu'
  if (tag == 'cxu') then
    write(*,*) "Cuting values of rho.dat on X-axis up to a cut:", cut
    do ispin = 1, nspin
      do i = cut + 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1, gridM(1, 3)
            N(I, J, K, ispin) = M(I, J, K, ispin)
          enddo
        enddo
      enddo
    enddo
  endif
! cut on Yaxis 'cy'
! down to a Ycut 'cyd'
  if (tag == 'cyd') then
    write(*,*) "Cuting values of rho.dat on Y-axis down to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, cut
          do k = 1, gridM(1, 3)
            N(i, j, k, ispin) = M(i, j, k, ispin)
          enddo
        enddo
       enddo
     enddo
   endif
! up to a Ycut 'cyu'
  if (tag == 'cyu') then
    write(*,*) "Cuting values of rho.dat on Y-axis up to a cut:", cut
    do ispin = 1, nspin
      do I = 1, gridM(1, 1)
        w = 1
        do J = cut + 1, gridM(1, 2)
          do k = 1, gridM(1, 3)
            N(i, w, k, ispin) = M(i, j, k, ispin)
          enddo
          w = w + 1
        enddo
      enddo
    enddo
  endif
! cut on Zaxis 'cz'
! down to a Zcut 'czd'
  if (tag == 'czd') then
    write(*,*) "Cuting values of rho.dat on Z-axis down to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1, cut
            N(i, j, k, ispin) = M (i, j, k, ispin)
          enddo
        enddo
      enddo
    enddo
  endif
! up to a Zcut 'czu'
  if (tag == 'czu') then
    write(*,*) "Cuting values of rho.dat on Z-axis up to a cut:", cut
    write(*,*) "sum inside:", sum(M)
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          w = 1
          do k = cut + 1, gridM(1, 3)
            N(i, j, w, ispin) = M(i, j, k, ispin)
            w = w + 1
         enddo
        enddo
       enddo
       enddo
      endif
!
! nullifying ------------------------!
!
! Nullify on Xaxis 'nx'
! down to a Xcut 'nxd'
  if (tag == 'nxd') then
    if (.not. consistency) then
      write(*,*) "For nullifying values arrays must be consistent!"
      stop
    endif
    write(*,*) "Nullifying values of rho.dat on X-axis down to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1, gridM(1,3)
            if (i .le. cut) then
              N(i, j, k, ispin) = M(i, j, k, ispin) 
            else
              N(i, j, k, nspin) = 0.0D0
            endif
          enddo
        enddo
      enddo
    enddo
  endif
! up to a Xcut 'nxu'
  if (tag == 'nxu') then
    if (.not. consistency) then
      write(*,*) "For nullifying values arrays must be consistent!"
      stop
    endif
    write(*,*) "Nullifying values of rho.dat on X-axis up to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1, gridM(1, 3)
            if (i .gt. cut) then
              N(i, j, k, ispin) = M(i, j, k, ispin)
            else
              N(i, j, k, ispin) = 0.0D0
           endif
          enddo
        enddo
      enddo
    enddo
  endif
! Nullify on Yaxis 'ny'
! down to a Ycut 'nyd'
  if (tag == 'nyd') then
    if (.not. consistency) then
      write(*,*) "For nullifying values arrays must be consistent!"
      stop
    endif
    write(*,*) "Nullifying values of rho.dat on Y-axis down to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1, gridM(1, 3)
            if (j .le. cut) then
              N(i, j, k, ispin) = M(i, j, k, ispin)
            else
              N(i, j, k, ispin) = 0.0D0
            endif
          enddo
        enddo
      enddo
    enddo
  endif
! up to a Ycut 'nyu'
  if (tag == 'nyu') then
    if (.not. consistency) then
      write(*,*) "For nullifying values arrays must be consistent!"
      stop
    endif
    write(*,*) "Nullifying values of rho.dat on Y-axis up to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1, gridM(1, 3)
            if (j .gt. cut) then
              N(i, j, k, ispin) = M(i, j, k, ispin)
            else
             N(i, j, k, ispin) = 0.0D0
            endif
          enddo
        enddo
      enddo
    enddo
  endif
! Nullify on Zaxis 'nz'
! down to a Zcut 'nzd'
  if (tag == 'nzd') then
    if (.not. consistency) then
      write(*,*) "For nullifying values arrays must be consistent!"
      stop
    endif
    write(*,*) "Nullifying values of rho.dat on Z-axis down to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1, gridM(1, 3)
            if (k .le. cut) then
              N(i, j, k, ispin) = M(i, j, k, ispin)
            else
              N(i, j, k, ispin) = 0.0D0
            endif
          enddo
        enddo
      enddo
    enddo
  endif
! up to a Zcut 'nzu'
  if (tag == 'nzu') then
    if (.not. consistency) then
      write(*,*) "For nullifying values arrays must be consistent!"
      stop
    endif
    write(*,*) "Nullifying values of rho.dat on Z-axis up to a cut:", cut
    do ispin = 1, nspin
      do i = 1, gridM(1, 1)
        do j = 1, gridM(1, 2)
          do k = 1,gridM(1, 3)
            if (k .gt. cut) then
               N(i, j, k, ispin) = M(i, j, k, ispin)
             else
               N(i, j, k, ispin) = 0.0D0
            endif
          enddo
        enddo
       enddo
     enddo
   endif
! -------------------------------------------------------------------- C
      end subroutine CPYARRAY
!----------------------------------------------------------------------!
end module utils 
