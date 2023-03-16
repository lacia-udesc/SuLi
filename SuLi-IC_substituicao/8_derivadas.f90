SUBROUTINE derivax(a1,dimx,dimy,dimz,campo_saida)

USE disc

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida
!real(8), dimension(dimx,dimy,dimz) :: a1

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i+1,j,k)-a1(i-1,j,k))/(2.*dx)
  enddo
 enddo
enddo

END SUBROUTINE derivax

!##############################################################

SUBROUTINE derivay(a1,dimx,dimy,dimz,campo_saida)

USE disc

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j+1,k)-a1(i,j-1,k))/(2.*dy)
  enddo
 enddo
enddo

END SUBROUTINE derivay

!##############################################################

SUBROUTINE derivaz(a1,dimx,dimy,dimz,campo_saida)

USE disc


IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))/(2.*dz)
  enddo
 enddo
enddo

END SUBROUTINE derivaz












!##############################################################

SUBROUTINE derivaxu2p(a1,dimx,dimy,dimz,campo_saida)

USE disc

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida
!real(8), dimension(dimx,dimy,dimz) :: a1



!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx-1
    campo_saida(i,j,k)=0.5*(-a1(i+2,j,k)+4*a1(i+1,j,k)-3*a1(i,j,k))/dx
  enddo
 enddo
enddo

i = dimx
do k=1,dimz
 do j=1,dimy
    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
  enddo
 enddo




END SUBROUTINE derivaxu2p

!##############################################################

SUBROUTINE derivayu2p(a1,dimx,dimy,dimz,campo_saida)

USE disc

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy-1
  do i=1,dimx
    campo_saida(i,j,k)=0.5*(-a1(i,j+2,k)+4*a1(i,j+1,k)-3*a1(i,j,k))/dy
  enddo
 enddo
enddo

j = dimy
do k=1,dimz
 do i=1,dimx
    campo_saida(i,j,k)=0.5*(a1(i,j+1,k)-a1(i,j-1,k))/dy
  enddo
 enddo


END SUBROUTINE derivayu2p

!##############################################################

SUBROUTINE derivazu2p(a1,dimx,dimy,dimz,campo_saida)

USE disc


IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz-1
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=0.5*(-a1(i,j,k+2)+4*a1(i,j,k+1)-3*a1(i,j,k))/dz
  enddo
 enddo
enddo

k = dimz
do j=1,dimy
 do i=1,dimx
    campo_saida(i,j,k)=0.5*(a1(i,j,k+1)-a1(i,j,k-1))/dz
  enddo
 enddo



END SUBROUTINE derivazu2p












!##############################################################

SUBROUTINE derivaxu2n(a1,dimx,dimy,dimz,campo_saida)

USE disc

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida
!real(8), dimension(dimx,dimy,dimz) :: a1


!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=2,dimx
    campo_saida(i,j,k)=0.5*(a1(i-2,j,k)-4*a1(i-1,j,k)+3*a1(i,j,k))/dx
  enddo
 enddo
enddo

i = 1
do k=1,dimz
 do j=1,dimy
    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
  enddo
 enddo




END SUBROUTINE derivaxu2n

!##############################################################

SUBROUTINE derivayu2n(a1,dimx,dimy,dimz,campo_saida)

USE disc

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida


!a1=campo_entrada
do k=1,dimz
 do j=2,dimy
  do i=1,dimx
    campo_saida(i,j,k)=0.5*(a1(i,j-2,k)-4*a1(i,j-1,k)+3*a1(i,j,k))/dy
  enddo
 enddo
enddo

j = 1
do k=1,dimz
 do i=1,dimx
    campo_saida(i,j,k)=0.5*(a1(i,j+1,k)-a1(i,j-1,k))/dy
  enddo
 enddo


END SUBROUTINE derivayu2n

!##############################################################

SUBROUTINE derivazu2n(a1,dimx,dimy,dimz,campo_saida)

USE disc


IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=2,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=0.5*(a1(i,j,k-2)-4*a1(i,j,k-1)+3*a1(i,j,k))/dz
  enddo
 enddo
enddo

k = 1
do j=1,dimy
 do i=1,dimx
    campo_saida(i,j,k)=0.5*(a1(i,j,k+1)-a1(i,j,k-1))/dz
  enddo
 enddo



END SUBROUTINE derivazu2n











!##############################################################

SUBROUTINE interpx_cf(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(dimx,dimy,dimz)  :: a1
real(8), intent(OUT), dimension(dimx+1,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=2,dimx
    campo_saida(i,j,k)=(a1(i,j,k)+a1(i-1,j,k))*0.5
  enddo
  i=1
   !campo_saida(i,j,k)=2.*campo_saida(i+1,j,k)-campo_saida(i+2,j,k)
   campo_saida(i,j,k)=campo_saida(i+1,j,k)
  i=dimx+1
   !campo_saida(i,j,k)=2.*campo_saida(i-1,j,k)-campo_saida(i-2,j,k)
   campo_saida(i,j,k)=campo_saida(i-1,j,k)
 enddo
enddo

END SUBROUTINE interpx_cf

!##############################################################

SUBROUTINE interpx_fc(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(dimx,dimy,dimz)  :: a1
real(8), intent(OUT), dimension(dimx-1,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx-1
    campo_saida(i,j,k)=(a1(i,j,k)+a1(i+1,j,k))*0.5
  enddo
 enddo
enddo


END SUBROUTINE interpx_fc

!##############################################################

SUBROUTINE interpy_cf(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(dimx,dimy,dimz)  :: a1
real(8), intent(OUT), dimension(dimx,dimy+1,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=2,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j-1,k))*0.5
  enddo
 enddo
enddo
do k=1,dimz
 do i=1,dimx
  j=1
   !campo_saida(i,j,k)=2.*campo_saida(i,j+1,k)-campo_saida(i,j+2,k)
   campo_saida(i,j,k)=campo_saida(i,j+1,k)
  j=dimy+1
   !campo_saida(i,j,k)=2.*campo_saida(i,j-1,k)-campo_saida(i,j-2,k)
   campo_saida(i,j,k)=campo_saida(i,j-1,k)
 enddo
enddo

END SUBROUTINE interpy_cf

!##############################################################

SUBROUTINE interpy_fc(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(dimx,dimy,dimz)  :: a1
real(8), intent(OUT), dimension(dimx,dimy-1,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy-1
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j+1,k))*0.5
  enddo
 enddo
enddo


END SUBROUTINE interpy_fc

!##############################################################

SUBROUTINE interpz_cf(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(dimx,dimy,dimz)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz+1) :: campo_saida

!a1=campo_entrada
do k=2,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j,k-1))*0.5
  enddo
 enddo
enddo
do j=1,dimy
 do i=1,dimx
  k=1
   !campo_saida(i,j,k)=2.*campo_saida(i,j,k+1)-campo_saida(i,j,k+2)
   campo_saida(i,j,k)=campo_saida(i,j,k+1)                
  k=dimz+1
   !campo_saida(i,j,k)=2.*campo_saida(i,j,k-1)-campo_saida(i,j,k-2)
   campo_saida(i,j,k)=campo_saida(i,j,k-1)

 enddo
enddo

END SUBROUTINE interpz_cf

!##############################################################

SUBROUTINE interpz_fc(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(dimx,dimy,dimz)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz-1) :: campo_saida

!a1=campo_entrada
do k=1,dimz-1
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j,k+1))*0.5
  enddo
 enddo
enddo

END SUBROUTINE interpz_fc
