! derivadas centradas
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

!####################

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

!####################

SUBROUTINE derivaz_no(a1,dimx,dimy,dimz,campo_saida)

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
    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))/(2.*dz(i,j,k))
  enddo
 enddo
enddo


END SUBROUTINE derivaz_no

!####################

SUBROUTINE derivaz_x(a1,dimx,dimy,dimz,campo_saida)

USE disc
USE dzs

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))/(2.*dzx(i,j,k))
  enddo
 enddo
enddo


END SUBROUTINE derivaz_x

!####################

SUBROUTINE derivaz_y(a1,dimx,dimy,dimz,campo_saida)

USE disc
USE dzs

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))/(2.*dzy(i,j,k))
  enddo
 enddo
enddo


END SUBROUTINE derivaz_y

!####################

SUBROUTINE derivaz_z(a1,dimx,dimy,dimz,campo_saida)

USE disc
USE dzs

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))/(2.*dzz(i,j,k))
  enddo
 enddo
enddo


END SUBROUTINE derivaz_z

!####################
!####################
!####################
!####################
!####################

SUBROUTINE interpx_cf(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro
USE cond

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
 enddo
enddo

if (ccx0 == 0) then
 do k=1,dimz
  do j=1,dimy
   i=1
    campo_saida(i,j,k)=(a1(dimx,j,k)+a1(1,j,k))*0.5
   i=dimx+1
    campo_saida(i,j,k)=(a1(dimx,j,k)+a1(1,j,k))*0.5
  enddo
 enddo
else
 do k=1,dimz
  do j=1,dimy
   i=1
    campo_saida(i,j,k)=2.*campo_saida(i+1,j,k)-campo_saida(i+2,j,k)
   i=dimx+1
    campo_saida(i,j,k)=2.*campo_saida(i-1,j,k)-campo_saida(i-2,j,k)
  enddo
 enddo
endif

END SUBROUTINE interpx_cf

!####################

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

!####################

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
   campo_saida(i,j,k)=2.*campo_saida(i,j+1,k)-campo_saida(i,j+2,k)
  j=dimy+1
   campo_saida(i,j,k)=2.*campo_saida(i,j-1,k)-campo_saida(i,j-2,k)
 enddo
enddo


END SUBROUTINE interpy_cf

!####################

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

!####################

SUBROUTINE interpz_cf(a1,dimx,dimy,dimz,campo_saida) !cf --> centro pra face fc--> face pro centro

USE disc

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
   campo_saida(i,j,k)=2.*campo_saida(i,j,k+1)-campo_saida(i,j,k+2)
			!a1(i,j,k+1)+(a1(i,j,k+1)-a1(i,j,k+2))&
                      !*dz(i,j,k+1)/dz(i,j,k)
  k=dimz+1
   campo_saida(i,j,k)=2.*campo_saida(i,j,k-1)-campo_saida(i,j,k-2)
			!a1(i,j,k-1)+(a1(i,j,k-1)-a1(i,j,k-2))&
                      !*dz(i,j,k-1)/dz(i,j,k)
 enddo
enddo


END SUBROUTINE interpz_cf

!####################

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


!####################

SUBROUTINE interpx_fp(a1,dimx,dimy,dimz,campo_saida) !fp --> face para ponto (para plotagem)
USE cond

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx,dimy+1,dimz+1) :: campo_saida

!a1=campo_entrada
do k=1,dimz+1
 do j=1,dimy+1
  do i=1,dimx
    campo_saida(i,j,k)=(a1(i,j,k) + a1(i,j-1,k) + a1(i,j,k-1) + a1(i,j-1,k-1)) * 0.25
  enddo
 enddo
enddo


END SUBROUTINE interpx_fp

!####################

SUBROUTINE interpy_fp(a1,dimx,dimy,dimz,campo_saida) !fp --> face para ponto (para plotagem)
USE cond

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx+1,dimy,dimz+1) :: campo_saida

!a1=campo_entrada
do k=1,dimz+1
 do j=1,dimy
  do i=1,dimx+1
    campo_saida(i,j,k)=(a1(i,j,k) + a1(i,j,k-1) + a1(i-1,j,k) + a1(i-1,j,k-1)) * 0.25
  enddo
 enddo
enddo


END SUBROUTINE interpy_fp


!####################

SUBROUTINE interpz_fp(a1,dimx,dimy,dimz,campo_saida) !fp --> face para ponto (para plotagem)
USE cond

IMPLICIT NONE
integer :: i, j, k
integer, intent(IN) :: dimx,dimy,dimz
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
real(8), intent(OUT), dimension(dimx+1,dimy+1,dimz) :: campo_saida

!a1=campo_entrada
do k=1,dimz
 do j=1,dimy+1
  do i=1,dimx+1
    campo_saida(i,j,k)=(a1(i,j,k) + a1(i,j-1,k) + a1(i-1,j,k) + a1(i-1,j-1,k)) * 0.25
  enddo
 enddo
enddo


END SUBROUTINE interpz_fp

!####################

SUBROUTINE interp2d_fp(a1,dimx,dimy,campo_saida) !fp --> face para ponto (para plotagem)
USE cond

IMPLICIT NONE
integer :: i, j
integer, intent(IN) :: dimx,dimy
real(8), intent(IN), dimension(0:dimx+1,0:dimy+1)  :: a1
real(8), intent(OUT), dimension(dimx+1,dimy+1) :: campo_saida

!a1=campo_entrada
do j=1,dimy+1
 do i=1,dimx+1
   campo_saida(i,j)=(a1(i,j) + a1(i,j-1) + a1(i-1,j) + a1(i-1,j-1)) * 0.25
 enddo
enddo


END SUBROUTINE interp2d_fp


