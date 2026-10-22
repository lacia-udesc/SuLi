SUBROUTINE derivax(a1,dimx,dimy,dimz,campo_saida)
! derivada centrada de segunda ordem na direção x

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida
	!real(8), dimension(dimx,dimy,dimz) :: a1

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i+1,j,k)-a1(i-1,j,k))*inv_2dx
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE derivax

!##############################################################

SUBROUTINE derivay(a1,dimx,dimy,dimz,campo_saida)
! derivada centrada de segunda ordem na direção y

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

    !a1=campo_entrada

	!$omp parallel do
    do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j+1,k)-a1(i,j-1,k))*inv_2dy
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE derivay

!##############################################################

SUBROUTINE derivaz(a1,dimx,dimy,dimz,inv2d,campo_saida)
! derivada centrada de segunda ordem na direção z

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: inv2d
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))*inv2d(k)
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE derivaz


!##############################################################

SUBROUTINE derivaxu2p(a1,dimx,dimy,dimz,campo_saida)
! derivada por frontward de segunda ordem na direção x (usado para upwind)

	USE disc
	USE cond, only : ccx0
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida
	!real(8), dimension(dimx,dimy,dimz) :: a1

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx-1
	    campo_saida(i,j,k)=(-a1(i+2,j,k)+4*a1(i+1,j,k)-3*a1(i,j,k))*inv_2dx
	  enddo
	 enddo
	enddo
    !$omp end parallel do


	!i = dimx
	!do k=1,dimz
	! do j=1,dimy
	!    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
	! enddo
	!enddo

	i = dimx
	if (ccx0 .ne. 0) then
	 !$omp parallel do
	  do k=1,dimz
	   do j=1,dimy
	    campo_saida(i,j,k)=(a1(i+1,j,k)-a1(i-1,j,k))*inv_2dx
	   enddo
	  enddo
	 !$omp end parallel do
	else
	!$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	    campo_saida(i,j,k)=(-a1(3,j,k)+4*a1(2,j,k)-3*a1(i,j,k))*inv_2dx
	 enddo
	enddo
	!$omp end parallel do
	endif

END SUBROUTINE derivaxu2p

!##############################################################

SUBROUTINE derivayu2p(a1,dimx,dimy,dimz,campo_saida)
! derivada por frontward de segunda ordem na direção y (usado para upwind)

	USE disc
	USE cond, only : ccy0
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=1,dimy-1
	  do i=1,dimx
	    campo_saida(i,j,k)=(-a1(i,j+2,k)+4*a1(i,j+1,k)-3*a1(i,j,k))*inv_2dy
	  enddo
	 enddo
	enddo
    !$omp end parallel do


	if (ccy0 .ne. 0) then
		j = dimy
		!$omp parallel do
		do k=1,dimz
		do i=1,dimx
			campo_saida(i,j,k)=(a1(i,j+1,k)-a1(i,j-1,k))*inv_2dy
		enddo
		enddo
		!$omp end parallel do
	else
		j = dimy
		!$omp parallel do
		do k=1,dimz
		do i=1,dimx
			campo_saida(i,j,k)=(-a1(i,3,k)+4*a1(i,2,k)-3*a1(i,j,k))*inv_2dy
		enddo
		enddo
		!$omp end parallel do
	endif

END SUBROUTINE derivayu2p

!##############################################################

SUBROUTINE derivazu2p(a1,dimx,dimy,dimz,inv2d,campo_saida)
! derivada por frontward de segunda ordem na direção z (usado para upwind)

!! NOTE QUE MESMO PRA MAIOR ORDEM SE USOU APENAS O DZ DO PONTO CALCULADO

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: inv2d
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz-1
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(-a1(i,j,k+2)+4*a1(i,j,k+1)-3*a1(i,j,k))*inv2d(k)
	  enddo
	 enddo
	enddo
    !$omp end parallel do

    k = dimz
    !$omp parallel do
	do j=1,dimy
	 do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))*inv2d(k)
	 enddo
	enddo
    !$omp end parallel do



END SUBROUTINE derivazu2p

!##############################################################

SUBROUTINE derivaxu2n(a1,dimx,dimy,dimz,campo_saida)
! derivada por backward de segunda ordem na direção x (usado para upwind)

	USE disc
	USE cond, only : ccx0
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN), dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida
	!real(8), dimension(dimx,dimy,dimz) :: a1

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	  do i=2,dimx
	    campo_saida(i,j,k)=(a1(i-2,j,k)-4*a1(i-1,j,k)+3*a1(i,j,k))*inv_2dx
	  enddo
	 enddo
	enddo
	!$omp end parallel do

	!i = 1
	!do k=1,dimz
	! do j=1,dimy
	!    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
	! enddo
	!enddo

	if (ccx0 .ne. 0) then
        i = 1
	!$omp parallel do
	  do k=1,dimz
	   do j=1,dimy
	    campo_saida(i,j,k)=(a1(i+1,j,k)-a1(i-1,j,k))*inv_2dx
	   enddo
	  enddo
	!$omp end parallel do
	else
	
	i = 1
	!$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	    campo_saida(i,j,k)=(a1(dimx-2,j,k)-4*a1(dimx-1,j,k)+3*a1(i,j,k))*inv_2dx
	 enddo
	enddo
	!$omp end parallel do
	endif

END SUBROUTINE derivaxu2n

!##############################################################

SUBROUTINE derivayu2n(a1,dimx,dimy,dimz,campo_saida)
! derivada por backward de segunda ordem na direção y (usado para upwind)

	USE disc
	USE cond, only : ccy0
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=2,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j-2,k)-4*a1(i,j-1,k)+3*a1(i,j,k))*inv_2dy
	  enddo
	 enddo
	enddo
    !$omp end parallel do

	j = 1
	if (ccy0 .ne. 0) then
		!$omp parallel do
		do k=1,dimz
		do i=1,dimx
			campo_saida(i,j,k)=(a1(i,j+1,k)-a1(i,j-1,k))*inv_2dy
		enddo
		enddo
		!$omp end parallel do
	else
		!$omp parallel do
		do k=1,dimz
		do i=1,dimx
			campo_saida(i,j,k)=(a1(i,dimy-2,k)-4*a1(i,dimy-1,k)+3*a1(i,j,k))*inv_2dy
		enddo
		enddo
		!$omp end parallel do
	endif

END SUBROUTINE derivayu2n

!##############################################################

SUBROUTINE derivazu2n(a1,dimx,dimy,dimz,inv2d,campo_saida)
! derivada por backward de segunda ordem na direção z (usado para upwind)

!! NOTE QUE MESMO PRA MAIOR ORDEM SE USOU APENAS O DZ DO PONTO CALCULADO

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: inv2d
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

	!a1=campo_entrada

    !$omp parallel do
	do k=2,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k-2)-4*a1(i,j,k-1)+3*a1(i,j,k))*inv2d(k)
	  enddo
	 enddo
	enddo
    !$omp end parallel do

	k = 1
    !$omp parallel do
	do j=1,dimy
	 do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))*inv2d(k)
	  enddo
	 enddo
    !$omp end parallel do

END SUBROUTINE derivazu2n


!##############################################################

SUBROUTINE interpx_cf(a1,dimx,dimy,dimz,campo_saida) 
! interpolação linear: cf --> centro pra face fc--> face pro centro

	USE cond, only : ccx0
	!$ USE omp_lib
	
	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx+1,dimy,dimz) :: campo_saida

    !a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	  do i=2,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i-1,j,k))*0.5d0
	  enddo
	 enddo
	enddo
    !$omp end parallel do

	if (ccx0 .ne. 0) then
		!$omp parallel do
		do k=1,dimz
		do j=1,dimy
			!i=1
			!campo_saida(1,j,k)=2.*campo_saida(2,j,k)-campo_saida(3,j,k)
			campo_saida(1,j,k)=campo_saida(2,j,k)
			!i=dimx+1
			!campo_saida(dimx+1,j,k)=2.*campo_saida(dimx,j,k)-campo_saida(dimx-1,j,k)
			campo_saida(dimx+1,j,k)=campo_saida(dimx,j,k)
		enddo
		enddo
		!$omp end parallel do
	else
		!$omp parallel do
		do k=1,dimz
		do j=1,dimy
			campo_saida(1,j,k)	=(a1(1,j,k)+a1(dimx,j,k))*0.5d0
			campo_saida(dimx+1,j,k) =(a1(1,j,k)+a1(dimx,j,k))*0.5d0
		enddo
		enddo
		!$omp end parallel do
	endif

END SUBROUTINE interpx_cf

!##############################################################

SUBROUTINE interpx_fc(a1,dimx,dimy,dimz,campo_saida) 
! interpolação linear: cf --> centro pra face fc--> face pro centro

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx-1,dimy,dimz) :: campo_saida

  	!a1=campo_entrada
	
	!$omp parallel do
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx-1
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i+1,j,k))*0.5d0
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpx_fc

!##############################################################

SUBROUTINE interpy_cf(a1,dimx,dimy,dimz,campo_saida)
! interpolação linear: cf --> centro pra face fc--> face pro centro

	USE cond, only : ccy0
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy+1,dimz) :: campo_saida

	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=2,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j-1,k))*0.5d0
	  enddo
	 enddo
	enddo
    !$omp end parallel do

	
	if (ccy0 .ne. 0) then

		!$omp parallel do
		do k=1,dimz
		do i=1,dimx
			!j=1
			!campo_saida(i,1,k)=2.*campo_saida(i,2,k)-campo_saida(i,3,k)
			campo_saida(i,1,k)=campo_saida(i,2,k)
			!j=dimy+1
			!campo_saida(i,dimy+1,k)=2.*campo_saida(i,dimy,k)-campo_saida(i,dimy-1,k)
			campo_saida(i,dimy+1,k)=campo_saida(i,dimy,k)
		enddo
		enddo
		!$omp end parallel do

	else
	
		!$omp parallel do
		do k=1,dimz
		do i=1,dimx
			campo_saida(i,1,k)	=(a1(i,1,k)+a1(i,dimy,k))*0.5d0
			campo_saida(i,dimy+1,k) =(a1(i,1,k)+a1(i,dimy,k))*0.5d0
		enddo
		enddo
		!$omp end parallel do


	endif
END SUBROUTINE interpy_cf

!##############################################################

SUBROUTINE interpy_fc(a1,dimx,dimy,dimz,campo_saida)
! interpolação linear: cf --> centro pra face fc--> face pro centro

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy-1,dimz) :: campo_saida


	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz
	 do j=1,dimy-1
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j+1,k))*0.5d0
	  enddo
	 enddo
	enddo
    !$omp end parallel do



END SUBROUTINE interpy_fc

!##############################################################

SUBROUTINE interpz_cf(a1,dimx,dimy,dimz,d,campo_saida) 
! interpolação linear: cf --> centro pra face fc--> face pro centro

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz+1) :: campo_saida

	!a1=campo_entrada

    !$omp parallel do
	do k=2,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)*d(k)+a1(i,j,k-1)*d(k-1))*inv_2dzz(k)
	  enddo
	 enddo
	enddo
    !$omp end parallel do

    !$omp parallel do
    do j=1,dimy
    do i=1,dimx
        !k=1
        !campo_saida(i,j,1)=2.*campo_saida(i,j,2)-campo_saida(i,j,3)
        campo_saida(i,j,1)=campo_saida(i,j,2)
        !k=dimz+1
        !campo_saida(i,j,dimz+1)=2.*campo_saida(i,j,dimz)-campo_saida(i,j,dimz-1)
        campo_saida(i,j,dimz+1)=campo_saida(i,j,dimz)
    enddo
	enddo
    !$omp end parallel do

END SUBROUTINE interpz_cf

!##############################################################

SUBROUTINE interpz_fc(a1,dimx,dimy,dimz,d,campo_saida) 
! interpolação linear: cf --> centro pra face fc--> face pro centro

    USE disc
    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz-1) :: campo_saida


	!a1=campo_entrada

    !$omp parallel do
	do k=1,dimz-1
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)*d(k)+a1(i,j,k+1)*d(k+1))*inv_dzzsum(k)
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpz_fc

!##############################################################

SUBROUTINE interpxy_ff(u,nx,ny,nz,nx1,ny1,uy) !ff --> face pra face  !velocidade x na face y

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,nx1,ny,ny1,nz
	real(8),intent(IN),dimension(0:nx1+1,0:ny+1,0:nz+1)  :: u
	real(8),intent(OUT),dimension(nx,ny1,nz) :: uy

    !$omp parallel do
	do k=1,nz
	 do j=1,ny1
	  do i=1,nx
		uy(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * 0.25d0
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpxy_ff

!##############################################################

SUBROUTINE interpxz_ff(u,nx,ny,nz,nx1,nz1,d,uz) !ff --> face pra face  !velocidade x na face z

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,nx1,ny,nz,nz1
	real(8),intent(IN),dimension(0:nz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(0:nx1+1,0:ny+1,0:nz+1)  :: u
	real(8),intent(OUT),dimension(nx,ny,nz1) :: uz


    !$omp parallel do
	do k=1,nz1
	 do j=1,ny
	  do i=1,nx
		uz(i,j,k) = ((u(i,j,k) + u(i+1,j,k))*d(k) + (u(i,j,k-1) + u(i+1,j,k-1))*d(k-1))*0.5/(d(k)+d(k-1))
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpxz_ff

!##############################################################

SUBROUTINE interpyx_ff(v,nx,ny,nz,nx1,ny1,vx) !ff --> face pra face  !velocidade y na face x

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,nx1,ny,ny1,nz
	real(8),intent(IN),dimension(0:nx+1,0:ny1+1,0:nz+1)  :: v
	real(8),intent(OUT),dimension(nx1,ny,nz) :: vx

    !$omp parallel do
	do k=1,nz
	 do j=1,ny
	  do i=1,nx1
		vx(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * 0.25d0
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpyx_ff

!##############################################################

SUBROUTINE interpyz_ff(v,nx,ny,nz,ny1,nz1,d,vz) !ff --> face pra face  !velocidade y na face z

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,ny,ny1,nz,nz1
	real(8),intent(IN),dimension(0:nz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(0:nx+1,0:ny1+1,0:nz+1)  :: v
	real(8),intent(OUT),dimension(nx,ny,nz1) :: vz

    !$omp parallel do
	do k=1,nz1
	 do j=1,ny
	  do i=1,nx
		vz(i,j,k) = ((v(i,j,k) + v(i,j+1,k))*d(k) + (v(i,j,k-1) + v(i,j+1,k-1))*d(k-1))*0.5/(d(k)+d(k-1))
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpyz_ff

!##############################################################

SUBROUTINE interpzx_ff(w,nx,ny,nz,nx1,nz1,d,wx) !ff --> face pra face  !velocidade z na face x

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,nx1,ny,nz,nz1
	real(8),intent(IN),dimension(0:nz1+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(0:nx+1,0:ny+1,0:nz1+1) :: w
	real(8),intent(OUT),dimension(nx1,ny,nz) :: wx


    !$omp parallel do
	do k=1,nz
	 do j=1,ny
	  do i=1,nx1
		wx(i,j,k) = ((w(i,j,k) + w(i-1,j,k))*d(k) + (w(i,j,k+1)  + w(i-1,j,k+1))*d(k+1))*0.5/(d(k)+d(k+1))
	  enddo
	 enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpzx_ff

!##############################################################

SUBROUTINE interpzy_ff(w,nx,ny,nz,ny1,nz1,d,wy) !ff --> face pra face  !velocidade z na face y

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,ny,ny1,nz,nz1
	real(8),intent(IN),dimension(0:nz1+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(0:nx+1,0:ny+1,0:nz1+1) :: w
	real(8),intent(OUT),dimension(nx,ny1,nz) :: wy


    !$omp parallel do
	do k=1,nz
        do j=1,ny1
        do i=1,nx
                wy(i,j,k) = ((w(i,j,k) + w(i,j-1,k))*d(k) + (w(i,j,k+1) + w(i,j-1,k+1))*d(k+1))*0.5/(d(k)+d(k+1))
        enddo
        enddo
	enddo
    !$omp end parallel do


END SUBROUTINE interpzy_ff

!##############################################################

SUBROUTINE fantasmas_cf(a1,dimx,dimy,dimz,campo_saida)
! interpolação linear: cf --> centro pra face fc--> face pro centro

	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: campo_saida

	!a1=campo_entrada
    !$omp parallel do
    do k=1,dimz
    do j=1,dimy
    do i=1,dimx
    	campo_saida(i,j,k) = a1(i,j,k)
    enddo
    enddo
    enddo
    !$omp end parallel do

	CALL cc_geral(campo_saida,dimx,dimy,dimz,1)	

END SUBROUTINE fantasmas_cf

!##############################################################

SUBROUTINE wenox(a,dimx,dimy,dimz,dx1,dphidxp,dphidxn,ihs)

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i,j,k,kk,ii,ihs,dimx,dimy,dimz
	real(8),intent(in) :: dx1
	real(8),intent(in),dimension(dimx,dimy,dimz) :: a

	real(8),intent(out),dimension(dimx,dimy,dimz) :: dphidxp,dphidxn
	real(8), parameter :: aux1 = 13.d0/12.d0
	real(8), parameter :: aux2 = 0.25d0
	real(8), parameter :: aux3 = 1.d0/6.d0
	real(8), parameter :: aux6 = 1.d-8
	real(8), parameter, dimension(3) :: auxx = (/0.1d0, 0.6d0, 0.3d0/)
	real(8) :: is1up, is2up, is3up, is1un, is2un, is3un
	real(8) :: al1up, al2up, al3up, al1un, al2un, al3un
	real(8) :: alsum_up, alsum_un
	real(8) :: om1up, om2up, om3up, om1un, om2un, om3un

	real(8),dimension(-2:dimx+3,dimy,dimz) :: phi1
	real(8),dimension(dimx+4,dimy,dimz)    :: un
	real(8),dimension(-3:dimx,dimy,dimz)   :: up


    phi1(1:dimx, 1:dimy, 1:dimz) = a(1:dimx, 1:dimy, 1:dimz)


	if (ihs == 0) then !Contorno periódico
    !$omp parallel do
		do k = 1, dimz
		do j = 1, dimy
			phi1(0,j,k)  = phi1(dimx,j,k)
			phi1(-1,j,k) = phi1(dimx-1,j,k)
			phi1(-2,j,k) = phi1(dimx-2,j,k)

			phi1(dimx+1,j,k) = phi1(1,j,k)
			phi1(dimx+2,j,k) = phi1(2,j,k)
			phi1(dimx+3,j,k) = phi1(3,j,k)

		enddo
		enddo
    !$omp end parallel do

	elseif (ihs == 1) then !Distance extrapolation
    !$omp parallel do
		do k = 1, dimz
		do j = 1, dimy
			phi1(0,j,k)  = 1./5. * (12.*phi1(1,j,k)  - 9.*phi1(2,j,k) + 2.*phi1(3,j,k) )
			phi1(-1,j,k) = 1./5. * (12.*phi1(0,j,k)  - 9.*phi1(1,j,k) + 2.*phi1(2,j,k) )
			phi1(-2,j,k) = 1./5. * (12.*phi1(-1,j,k) - 9.*phi1(0,j,k) + 2.*phi1(1,j,k) )


			phi1(dimx+1,j,k) = 1./11. * (18.*phi1(dimx,j,k)   - 9.*phi1(dimx-1,j,k) + 2.*phi1(dimx-2,j,k))
			phi1(dimx+2,j,k) = 1./11. * (18.*phi1(dimx+1,j,k) - 9.*phi1(dimx,j,k)   + 2.*phi1(dimx-1,j,k))
			phi1(dimx+3,j,k) = 1./11. * (18.*phi1(dimx+2,j,k) - 9.*phi1(dimx+1,j,k) + 2.*phi1(dimx,j,k)  )

			!Dirichlet
			!phi1(1,j,k)  = 0.
			!phi1(0,j,k)  = 1./11. * (12.*phi1(1,j,k)  - 3.*phi1(2,j,k) + 2.*phi1(3,j,k) )
			!phi1(-1,j,k) = 1./11. * (12.*phi1(0,j,k)  - 3.*phi1(1,j,k) + 2.*phi1(2,j,k) )
			!phi1(-2,j,k) = 1./11. * (12.*phi1(-1,j,k) - 3.*phi1(0,j,k) + 2.*phi1(1,j,k) )

			!phi1(i,dimy,k) = 0.
			!phi1(dimx+1,j,k) = 1./11. * (24.*phi1(dimx,j,k)   - 15.*phi1(dimx-1,j,k) + 2.*phi1(dimx-2,j,k))
			!phi1(dimx+2,j,k) = 1./11. * (24.*phi1(dimx+1,j,k) - 15.*phi1(dimx,j,k)   + 2.*phi1(dimx-1,j,k))
			!phi1(dimx+3,j,k) = 1./11. * (24.*phi1(dimx+2,j,k) - 15.*phi1(dimx+1,j,k) + 2.*phi1(dimx,j,k)  )
  			!Dirichlet
		enddo
		enddo
    !$omp end parallel do

	elseif (ihs == 2) then !Derivative zero
    !$omp parallel do
		do k = 1, dimz
		do j = 1, dimy

			phi1(0,j,k)  = 1./11. * (18.*phi1(1,j,k)  - 9.*phi1(2,j,k) + 2.*phi1(3,j,k) )
			phi1(-1,j,k) = 1./11. * (18.*phi1(0,j,k)  - 9.*phi1(1,j,k) + 2.*phi1(2,j,k) )
			phi1(-2,j,k) = 1./11. * (18.*phi1(-1,j,k) - 9.*phi1(0,j,k) + 2.*phi1(1,j,k) )

			phi1(dimx+1,j,k) = 1./11. * (18.*phi1(dimx,j,k)   - 9.*phi1(dimx-1,j,k) + 2.*phi1(dimx-2,j,k))
			phi1(dimx+2,j,k) = 1./11. * (18.*phi1(dimx+1,j,k) - 9.*phi1(dimx,j,k)   + 2.*phi1(dimx-1,j,k))
			phi1(dimx+3,j,k) = 1./11. * (18.*phi1(dimx+2,j,k) - 9.*phi1(dimx+1,j,k) + 2.*phi1(dimx,j,k)  )

		enddo
		enddo
    !$omp end parallel do
	endif


    !$omp parallel do
	do k = 1, dimz
	do j = 1, dimy
	do i =-3, dimx
	        up(i,j,k)=(phi1(i+3,j,k)-phi1(i+2,j,k))/dx1
	enddo
	enddo
	enddo
    !$omp end parallel do


    !$omp parallel do
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx+4
	        un(i,j,k)=(phi1(i-2,j,k)-phi1(i-3,j,k))/dx1
	enddo
	enddo
	enddo
    !$omp end parallel do




    !$omp parallel do private(is1up,is2up,is3up,is1un,is2un,is3un,al1up,al2up,al3up,al1un,al2un,&
	!$omp al3un,alsum_up,alsum_un,om1up,om2up,om3up,om1un,om2un,om3un)
	
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
	is1up = aux1*(up(i,j,k)-2*up(i-1,j,k)+up(i-2,j,k))*(up(i,j,k)-2*up(i-1,j,k)+up(i-2,j,k)) &
	      + aux2*(up(i,j,k)-4*up(i-1,j,k)+3*up(i-2,j,k))*(up(i,j,k)-4*up(i-1,j,k)+3*up(i-2,j,k))
	is1un = aux1*(un(i,j,k)-2*un(i+1,j,k)+un(i+2,j,k))*(un(i,j,k)-2*un(i+1,j,k)+un(i+2,j,k)) &
	      + aux2*(un(i,j,k)-4*un(i+1,j,k)+3*un(i+2,j,k))*(un(i,j,k)-4*un(i+1,j,k)+3*un(i+2,j,k))
	is2up = aux1*(up(i-1,j,k)-2*up(i-2,j,k)+up(i-3,j,k))*(up(i-1,j,k)-2*up(i-2,j,k)+up(i-3,j,k)) &
	      + aux2*(up(i-1,j,k)-up(i-3,j,k))*(up(i-1,j,k)-up(i-3,j,k))
	is2un = aux1*(un(i+1,j,k)-2*un(i+2,j,k)+un(i+3,j,k))*(un(i+1,j,k)-2*un(i+2,j,k)+un(i+3,j,k)) &
	      + aux2*(un(i+1,j,k)-un(i+3,j,k))*(un(i+1,j,k)-un(i+3,j,k))
	is3up = aux1*(up(i-2,j,k)-2*up(i-3,j,k)+up(i-4,j,k))*(up(i-2,j,k)-2*up(i-3,j,k)+up(i-4,j,k)) &
	      + aux2*(3*up(i-2,j,k)-4*up(i-3,j,k)+up(i-4,j,k))*(3*up(i-2,j,k)-4*up(i-3,j,k)+up(i-4,j,k))
	is3un = aux1*(un(i+2,j,k)-2*un(i+3,j,k)+un(i+4,j,k))*(un(i+2,j,k)-2*un(i+3,j,k)+un(i+4,j,k)) &
	      + aux2*(3*un(i+2,j,k)-4*un(i+3,j,k)+un(i+4,j,k))*(3*un(i+2,j,k)-4*un(i+3,j,k)+un(i+4,j,k))

	al1up = auxx(1)/((aux6+is1up)*(aux6+is1up))
	al2up = auxx(2)/((aux6+is2up)*(aux6+is2up))
	al3up = auxx(3)/((aux6+is3up)*(aux6+is3up))
	al1un = auxx(1)/((aux6+is1un)*(aux6+is1un))
	al2un = auxx(2)/((aux6+is2un)*(aux6+is2un))
	al3un = auxx(3)/((aux6+is3un)*(aux6+is3un))

	alsum_up = al1up + al2up + al3up
	alsum_un = al1un + al2un + al3un

	om1up = al1up/alsum_up
	om2up = al2up/alsum_up
	om3up = al3up/alsum_up
	om1un = al1un/alsum_un
	om2un = al2un/alsum_un
	om3un = al3un/alsum_un

	dphidxp(i,j,k) = aux3*(om1up*(2*up(i,j,k)-7*up(i-1,j,k)+11*up(i-2,j,k)) + &
	om2up*(-up(i-1,j,k)+5*up(i-2,j,k)+2*up(i-3,j,k)) + om3up*(2*up(i-2,j,k)+5*up(i-3,j,k)-up(i-4,j,k)))
	dphidxn(i,j,k) = aux3*(om1un*(2*un(i,j,k)-7*un(i+1,j,k)+11*un(i+2,j,k)) + &
	om2un*(-un(i+1,j,k)+5*un(i+2,j,k)+2*un(i+3,j,k)) + om3un*(2*un(i+2,j,k)+5*un(i+3,j,k)-un(i+4,j,k)))
	enddo
	enddo
	enddo
    !$omp end parallel do


END SUBROUTINE wenox

!##############################################################

SUBROUTINE wenoy(b,dimx,dimy,dimz,dy1,dphidyp,dphidyn,ihs)

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i,j,k,kk,ii,ihs,dimx,dimy,dimz
	real(8),intent(in) :: dy1
	real(8),intent(in),dimension(dimx,dimy,dimz) :: b
				
	real(8),intent(out),dimension(dimx,dimy,dimz) :: dphidyp,dphidyn
	real(8), parameter :: aux1 = 13.d0/12.d0
	real(8), parameter :: aux2 = 0.25d0
	real(8), parameter :: aux3 = 1.d0/6.d0
	real(8), parameter :: aux6 = 1.d-8
	real(8), parameter, dimension(3) :: auxx = (/0.1d0, 0.6d0, 0.3d0/)
	real(8) :: is1vp, is2vp, is3vp, is1vn, is2vn, is3vn
	real(8) :: al1vp, al2vp, al3vp, al1vn, al2vn, al3vn
	real(8) :: alsum_vp, alsum_vn
	real(8) :: om1vp, om2vp, om3vp, om1vn, om2vn, om3vn
	
	real(8),dimension(dimx,-2:dimy+3,dimz) :: phi1
	real(8),dimension(dimx,dimy+4,dimz)    :: vn
	real(8),dimension(dimx,-3:dimy,dimz)   :: vp


    phi1(1:dimx, 1:dimy, 1:dimz) = b(1:dimx, 1:dimy, 1:dimz)


	if (ihs == 0) then !Contorno periódico
    !$omp parallel do
		do k = 1, dimz
		do i = 1, dimx
			phi1(i,0,k)  = phi1(i,dimy,k)
			phi1(i,-1,k) = phi1(i,dimy-1,k)
			phi1(i,-2,k) = phi1(i,dimy-2,k)

			phi1(i,dimy+1,k) = phi1(i,1,k)
			phi1(i,dimy+2,k) = phi1(i,2,k)
			phi1(i,dimy+3,k) = phi1(i,3,k)
			
		enddo
		enddo
    !$omp end parallel do
		
	elseif (ihs == 1) then !Distance extrapolation
    !$omp parallel do
		do k = 1, dimz
		do i = 1, dimx
			phi1(i,0,k)  = 1./5. * (12.*phi1(i,1,k)  - 9.*phi1(i,2,k) + 2.*phi1(i,3,k) )
			phi1(i,-1,k) = 1./5. * (12.*phi1(i,0,k)  - 9.*phi1(i,1,k) + 2.*phi1(i,2,k) )
			phi1(i,-2,k) = 1./5. * (12.*phi1(i,-1,k) - 9.*phi1(i,0,k) + 2.*phi1(i,1,k) )


			phi1(i,dimy+1,k) = 1./11. * (18.*phi1(i,dimy,k)   - 9.*phi1(i,dimy-1,k) + 2.*phi1(i,dimy-2,k))
			phi1(i,dimy+2,k) = 1./11. * (18.*phi1(i,dimy+1,k) - 9.*phi1(i,dimy,k)   + 2.*phi1(i,dimy-1,k))
			phi1(i,dimy+3,k) = 1./11. * (18.*phi1(i,dimy+2,k) - 9.*phi1(i,dimy+1,k) + 2.*phi1(i,dimy,k)  )

			!Dirichlet
			!phi1(i,1,k)  = 0.
			!phi1(i,0,k)  = 1./11. * (12.*phi1(i,1,k)  - 3.*phi1(i,2,k) + 2.*phi1(i,3,k) )
			!phi1(i,-1,k) = 1./11. * (12.*phi1(i,0,k)  - 3.*phi1(i,1,k) + 2.*phi1(i,2,k) )
			!phi1(i,-2,k) = 1./11. * (12.*phi1(i,-1,k) - 3.*phi1(i,0,k) + 2.*phi1(i,1,k) )

			!phi1(i,dimy,k) = 0.
			!phi1(i,dimy+1,k) = 1./11. * (24.*phi1(i,dimy,k)   - 15.*phi1(i,dimy-1,k) + 2.*phi1(i,dimy-2,k))
			!phi1(i,dimy+2,k) = 1./11. * (24.*phi1(i,dimy+1,k) - 15.*phi1(i,dimy,k)   + 2.*phi1(i,dimy-1,k))
			!phi1(i,dimy+3,k) = 1./11. * (24.*phi1(i,dimy+2,k) - 15.*phi1(i,dimy+1,k) + 2.*phi1(i,dimy,k)  )
			!Dirichlet
		enddo
		enddo
    !$omp end parallel do
		
	elseif (ihs == 2) then !Derivative zero
    !$omp parallel do
		do k = 1, dimz
		do i = 1, dimx
					
			phi1(i,0,k)  = 1./11. * (18.*phi1(i,1,k)  - 9.*phi1(i,2,k) + 2.*phi1(i,3,k) )
			phi1(i,-1,k) = 1./11. * (18.*phi1(i,0,k)  - 9.*phi1(i,1,k) + 2.*phi1(i,2,k) )
			phi1(i,-2,k) = 1./11. * (18.*phi1(i,-1,k) - 9.*phi1(i,0,k) + 2.*phi1(i,1,k) )

			phi1(i,dimy+1,k) = 1./11. * (18.*phi1(i,dimy,k)   - 9.*phi1(i,dimy-1,k) + 2.*phi1(i,dimy-2,k))
			phi1(i,dimy+2,k) = 1./11. * (18.*phi1(i,dimy+1,k) - 9.*phi1(i,dimy,k)   + 2.*phi1(i,dimy-1,k))
			phi1(i,dimy+3,k) = 1./11. * (18.*phi1(i,dimy+2,k) - 9.*phi1(i,dimy+1,k) + 2.*phi1(i,dimy,k)  )
			
		enddo
		enddo
    !$omp end parallel do
	endif
    
    !$omp parallel do
	do k = 1, dimz
	do j = -3, dimy
	do i = 1, dimx
                vp(i,j,k)=(phi1(i,j+3,k)-phi1(i,j+2,k))/dy1
	enddo
	enddo
	enddo
    !$omp end parallel do


    !$omp parallel do
	do k = 1, dimz
	do j = 1, dimy+4
	do i = 1, dimx
                vn(i,j,k)=(phi1(i,j-2,k)-phi1(i,j-3,k))/dy1
	enddo
	enddo
	enddo
    !$omp end parallel do

    !$omp parallel do private(is1vp,is2vp,is3vp,is1vn,is2vn,is3vn,al1vp,al2vp,al3vp,al1vn,al2vn,al3vn,&
	!$omp alsum_vp,alsum_vn,om1vp,om2vp,om3vp,om1vn,om2vn,om3vn)
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
	is1vp = aux1*(vp(i,j,k)-2*vp(i,j-1,k)+vp(i,j-2,k))*(vp(i,j,k)-2*vp(i,j-1,k)+vp(i,j-2,k)) &
	      + aux2*(vp(i,j,k)-4*vp(i,j-1,k)+3*vp(i,j-2,k))*(vp(i,j,k)-4*vp(i,j-1,k)+3*vp(i,j-2,k))
	is1vn = aux1*(vn(i,j,k)-2*vn(i,j+1,k)+vn(i,j+2,k))*(vn(i,j,k)-2*vn(i,j+1,k)+vn(i,j+2,k)) &
	      + aux2*(vn(i,j,k)-4*vn(i,j+1,k)+3*vn(i,j+2,k))*(vn(i,j,k)-4*vn(i,j+1,k)+3*vn(i,j+2,k))
	is2vp = aux1*(vp(i,j-1,k)-2*vp(i,j-2,k)+vp(i,j-3,k))*(vp(i,j-1,k)-2*vp(i,j-2,k)+vp(i,j-3,k)) &
	      + aux2*(vp(i,j-1,k)-vp(i,j-3,k))*(vp(i,j-1,k)-vp(i,j-3,k))
	is2vn = aux1*(vn(i,j+1,k)-2*vn(i,j+2,k)+vn(i,j+3,k))*(vn(i,j+1,k)-2*vn(i,j+2,k)+vn(i,j+3,k)) &
	      + aux2*(vn(i,j+1,k)-vn(i,j+3,k))*(vn(i,j+1,k)-vn(i,j+3,k))
	is3vp = aux1*(vp(i,j-2,k)-2*vp(i,j-3,k)+vp(i,j-4,k))*(vp(i,j-2,k)-2*vp(i,j-3,k)+vp(i,j-4,k)) &
	      + aux2*(3*vp(i,j-2,k)-4*vp(i,j-3,k)+vp(i,j-4,k))*(3*vp(i,j-2,k)-4*vp(i,j-3,k)+vp(i,j-4,k))
	is3vn = aux1*(vn(i,j+2,k)-2*vn(i,j+3,k)+vn(i,j+4,k))*(vn(i,j+2,k)-2*vn(i,j+3,k)+vn(i,j+4,k)) &
	      + aux2*(3*vn(i,j+2,k)-4*vn(i,j+3,k)+vn(i,j+4,k))*(3*vn(i,j+2,k)-4*vn(i,j+3,k)+vn(i,j+4,k))

	al1vp = auxx(1)/((aux6+is1vp)*(aux6+is1vp))
	al2vp = auxx(2)/((aux6+is2vp)*(aux6+is2vp))
	al3vp = auxx(3)/((aux6+is3vp)*(aux6+is3vp))
	al1vn = auxx(1)/((aux6+is1vn)*(aux6+is1vn))
	al2vn = auxx(2)/((aux6+is2vn)*(aux6+is2vn))
	al3vn = auxx(3)/((aux6+is3vn)*(aux6+is3vn))

	alsum_vp = al1vp + al2vp + al3vp
	alsum_vn = al1vn + al2vn + al3vn

	om1vp = al1vp/alsum_vp
	om2vp = al2vp/alsum_vp
	om3vp = al3vp/alsum_vp
	om1vn = al1vn/alsum_vn
	om2vn = al2vn/alsum_vn
	om3vn = al3vn/alsum_vn

	dphidyp(i,j,k) = aux3*(om1vp*(2*vp(i,j,k)-7*vp(i,j-1,k)+11*vp(i,j-2,k)) + &
	om2vp*(-vp(i,j-1,k)+5*vp(i,j-2,k)+2*vp(i,j-3,k)) + om3vp*(2*vp(i,j-2,k)+5*vp(i,j-3,k)-vp(i,j-4,k)))
	dphidyn(i,j,k) = aux3*(om1vn*(2*vn(i,j,k)-7*vn(i,j+1,k)+11*vn(i,j+2,k)) + &
	om2vn*(-vn(i,j+1,k)+5*vn(i,j+2,k)+2*vn(i,j+3,k)) + om3vn*(2*vn(i,j+2,k)+5*vn(i,j+3,k)-vn(i,j+4,k)))
	enddo
	enddo
	enddo
    !$omp end parallel do


END SUBROUTINE wenoy

!##############################################################

SUBROUTINE wenoz(d,dimx,dimy,dimz,dz1,dphidzp,dphidzn,ihs)

    !$ USE omp_lib

	IMPLICIT NONE
	integer :: i,j,k,kk,ii,ihs,dimx,dimy,dimz
	real(8),intent(in),dimension(0:dimz+1) :: dz1 !ter certeza desses limites!
	real(8),intent(in),dimension(dimx,dimy,dimz) :: d
				
	real(8),intent(out),dimension(dimx,dimy,dimz) :: dphidzp,dphidzn
	real(8), parameter :: aux1 = 13.d0/12.d0
	real(8), parameter :: aux2 = 0.25d0
	real(8), parameter :: aux3 = 1.d0/6.d0
	real(8), parameter :: aux6 = 1.d-8
	real(8), parameter, dimension(3) :: auxx = (/0.1d0, 0.6d0, 0.3d0/)
	real(8) :: is1wp, is2wp, is3wp, is1wn, is2wn, is3wn
	real(8) :: al1wp, al2wp, al3wp, al1wn, al2wn, al3wn
	real(8) :: alsum_wp, alsum_wn
	real(8) :: om1wp, om2wp, om3wp, om1wn, om2wn, om3wn

	real(8),dimension(-2:dimz+3) :: dz11	
	real(8),dimension(dimx,dimy,-2:dimz+3) :: phi1
	real(8),dimension(dimx,dimy,dimz+4)    :: wn
	real(8),dimension(dimx,dimy,-3:dimz)   :: wp


    phi1(1:dimx, 1:dimy, 1:dimz) = d(1:dimx, 1:dimy, 1:dimz)

	dz11(0:dimz+1) = dz1(0:dimz+1)
	
	if (ihs == 0) then !Contorno periódico
		dz11(0)  = dz11(dimz)
		dz11(-1) = dz11(dimz-1)
		dz11(-2) = dz11(dimz-2)
		dz11(dimz+1) = dz11(1)
		dz11(dimz+2) = dz11(2)
		dz11(dimz+3) = dz11(3)
    !$omp parallel do
		do j = 1, dimy
		do i = 1, dimx
			phi1(i,j,0)  = phi1(i,j,dimz)
			phi1(i,j,-1) = phi1(i,j,dimz-1)
			phi1(i,j,-2) = phi1(i,j,dimz-2)

			phi1(i,j,dimz+1) = phi1(i,j,1)
			phi1(i,j,dimz+2) = phi1(i,j,2)
			phi1(i,j,dimz+3) = phi1(i,j,3)
		enddo
		enddo
    !$omp end parallel do
		
	elseif (ihs == 1) then !Velocidade zero com derivada de 1a ordem 
	! no-slip de velocidade = 0 ou free-slip para velocidade perpendicular
		! dz11 ghost cells: independent of (i,j), compute serially first
		dz11(0)  = 1./5. * (12.*dz11(1)  - 9.*dz11(2) + 2.*dz11(3) )
		dz11(-1) = 1./5. * (12.*dz11(0)  - 9.*dz11(1) + 2.*dz11(2) )
		dz11(-2) = 1./5. * (12.*dz11(-1) - 9.*dz11(0) + 2.*dz11(1) )

		dz11(dimz+1) = 1./11. * (18.*dz11(dimz)   - 9.*dz11(dimz-1) + 2.*dz11(dimz-2))
		dz11(dimz+2) = 1./11. * (18.*dz11(dimz+1) - 9.*dz11(dimz)   + 2.*dz11(dimz-1))
		dz11(dimz+3) = 1./11. * (18.*dz11(dimz+2) - 9.*dz11(dimz+1) + 2.*dz11(dimz)  )
    !$omp parallel do
		do j = 1, dimy
		do i = 1, dimx
			!Dirichlet
			!phi1(i,j,1)  = 0.
			!phi1(i,j,0)  = 1./11. * (12.*phi1(i,j,1)  - 3.*phi1(i,j,2) + 2.*phi1(i,j,3) )
			!phi1(i,j,-1) = 1./11. * (12.*phi1(i,j,0)  - 3.*phi1(i,j,1) + 2.*phi1(i,j,2) )
			!phi1(i,j,-2) = 1./11. * (12.*phi1(i,j,-1) - 3.*phi1(i,j,0) + 2.*phi1(i,j,1) )

			!phi1(i,j,dimz) = 0.
			!phi1(i,j,dimz+1) = 1./11. * (24.*phi1(i,j,dimz)   - 15.*phi1(i,j,dimz-1) + 2.*phi1(i,j,dimz-2))
			!phi1(i,j,dimz+2) = 1./11. * (24.*phi1(i,j,dimz+1) - 15.*phi1(i,j,dimz)   + 2.*phi1(i,j,dimz-1))
			!phi1(i,j,dimz+3) = 1./11. * (24.*phi1(i,j,dimz+2) - 15.*phi1(i,j,dimz+1) + 2.*phi1(i,j,dimz)  )

			!Dirichlet
			

			phi1(i,j,0)  = 1./5. * (12.*phi1(i,j,1)  - 9.*phi1(i,j,2) + 2.*phi1(i,j,3) )
			phi1(i,j,-1) = 1./5. * (12.*phi1(i,j,0)  - 9.*phi1(i,j,1) + 2.*phi1(i,j,2) )
			phi1(i,j,-2) = 1./5. * (12.*phi1(i,j,-1) - 9.*phi1(i,j,0) + 2.*phi1(i,j,1) )

			phi1(i,j,dimz+1) = 1./11. * (18.*phi1(i,j,dimz)   - 9.*phi1(i,j,dimz-1) + 2.*phi1(i,j,dimz-2))
			phi1(i,j,dimz+2) = 1./11. * (18.*phi1(i,j,dimz+1) - 9.*phi1(i,j,dimz)   + 2.*phi1(i,j,dimz-1))
			phi1(i,j,dimz+3) = 1./11. * (18.*phi1(i,j,dimz+2) - 9.*phi1(i,j,dimz+1) + 2.*phi1(i,j,dimz)  )			
		enddo
		enddo
    !$omp end parallel do
		
	elseif (ihs == 2) then !Derivative zero
		! dz11 ghost cells: independent of (i,j), compute serially first
		dz11(0)  = 1./11. * (18.*dz11(1)  - 9.*dz11(2) + 2.*dz11(3) )
		dz11(-1) = 1./11. * (18.*dz11(0)  - 9.*dz11(1) + 2.*dz11(2) )
		dz11(-2) = 1./11. * (18.*dz11(-1) - 9.*dz11(0) + 2.*dz11(1) )

		dz11(dimz+1) = 1./11. * (18.*dz11(dimz)   - 9.*dz11(dimz-1) + 2.*dz11(dimz-2))
		dz11(dimz+2) = 1./11. * (18.*dz11(dimz+1) - 9.*dz11(dimz)   + 2.*dz11(dimz-1))
		dz11(dimz+3) = 1./11. * (18.*dz11(dimz+2) - 9.*dz11(dimz+1) + 2.*dz11(dimz)  )
    !$omp parallel do
		do j = 1, dimy
		do i = 1, dimx
			
			phi1(i,j,0)  = 1./11. * (18.*phi1(i,j,1)  - 9.*phi1(i,j,2) + 2.*phi1(i,j,3) )
			phi1(i,j,-1) = 1./11. * (18.*phi1(i,j,0)  - 9.*phi1(i,j,1) + 2.*phi1(i,j,2) )
			phi1(i,j,-2) = 1./11. * (18.*phi1(i,j,-1) - 9.*phi1(i,j,0) + 2.*phi1(i,j,1) )

			phi1(i,j,dimz+1) = 1./11. * (18.*phi1(i,j,dimz)   - 9.*phi1(i,j,dimz-1) + 2.*phi1(i,j,dimz-2))
			phi1(i,j,dimz+2) = 1./11. * (18.*phi1(i,j,dimz+1) - 9.*phi1(i,j,dimz)   + 2.*phi1(i,j,dimz-1))
			phi1(i,j,dimz+3) = 1./11. * (18.*phi1(i,j,dimz+2) - 9.*phi1(i,j,dimz+1) + 2.*phi1(i,j,dimz)  )
		enddo
		enddo
    !$omp end parallel do
	endif
    

    !$omp parallel do
	do k = -3, dimz
	do j = 1, dimy
	do i = 1, dimx
	wp(i,j,k)=(phi1(i,j,k+3)-phi1(i,j,k+2))/dz11(k+3)
	enddo
	enddo
	enddo
    !$omp end parallel do


    !$omp parallel do
	do k = 1, dimz+4
	do j = 1, dimy
	do i = 1, dimx
	wn(i,j,k)=(phi1(i,j,k-2)-phi1(i,j,k-3))/dz11(k-2)
	enddo
	enddo
	enddo
    !$omp end parallel do

    !$omp parallel do private(is1wp,is2wp,is3wp,is1wn,is2wn,is3wn,al1wp,al2wp,al3wp,al1wn,al2wn,al3wn,&
	!$omp alsum_wp,alsum_wn,om1wp,om2wp,om3wp,om1wn,om2wn,om3wn)
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
	is1wp = aux1*(wp(i,j,k)-2*wp(i,j,k-1)+wp(i,j,k-2))*(wp(i,j,k)-2*wp(i,j,k-1)+wp(i,j,k-2)) &
	      + aux2*(wp(i,j,k)-4*wp(i,j,k-1)+3*wp(i,j,k-2))*(wp(i,j,k)-4*wp(i,j,k-1)+3*wp(i,j,k-2))
	is1wn = aux1*(wn(i,j,k)-2*wn(i,j,k+1)+wn(i,j,k+2))*(wn(i,j,k)-2*wn(i,j,k+1)+wn(i,j,k+2)) &
	      + aux2*(wn(i,j,k)-4*wn(i,j,k+1)+3*wn(i,j,k+2))*(wn(i,j,k)-4*wn(i,j,k+1)+3*wn(i,j,k+2))
	is2wp = aux1*(wp(i,j,k-1)-2*wp(i,j,k-2)+wp(i,j,k-3))*(wp(i,j,k-1)-2*wp(i,j,k-2)+wp(i,j,k-3)) &
	      + aux2*(wp(i,j,k-1)-wp(i,j,k-3))*(wp(i,j,k-1)-wp(i,j,k-3))
	is2wn = aux1*(wn(i,j,k+1)-2*wn(i,j,k+2)+wn(i,j,k+3))*(wn(i,j,k+1)-2*wn(i,j,k+2)+wn(i,j,k+3)) &
	      + aux2*(wn(i,j,k+1)-wn(i,j,k+3))*(wn(i,j,k+1)-wn(i,j,k+3))
	is3wp = aux1*(wp(i,j,k-2)-2*wp(i,j,k-3)+wp(i,j,k-4))*(wp(i,j,k-2)-2*wp(i,j,k-3)+wp(i,j,k-4)) &
	      + aux2*(3*wp(i,j,k-2)-4*wp(i,j,k-3)+wp(i,j,k-4))*(3*wp(i,j,k-2)-4*wp(i,j,k-3)+wp(i,j,k-4))
	is3wn = aux1*(wn(i,j,k+2)-2*wn(i,j,k+3)+wn(i,j,k+4))*(wn(i,j,k+2)-2*wn(i,j,k+3)+wn(i,j,k+4)) &
	      + aux2*(3*wn(i,j,k+2)-4*wn(i,j,k+3)+wn(i,j,k+4))*(3*wn(i,j,k+2)-4*wn(i,j,k+3)+wn(i,j,k+4))

	al1wp = auxx(1)/((aux6+is1wp)*(aux6+is1wp))
	al2wp = auxx(2)/((aux6+is2wp)*(aux6+is2wp))
	al3wp = auxx(3)/((aux6+is3wp)*(aux6+is3wp))
	al1wn = auxx(1)/((aux6+is1wn)*(aux6+is1wn))
	al2wn = auxx(2)/((aux6+is2wn)*(aux6+is2wn))
	al3wn = auxx(3)/((aux6+is3wn)*(aux6+is3wn))

	alsum_wp = al1wp + al2wp + al3wp
	alsum_wn = al1wn + al2wn + al3wn

	om1wp = al1wp/alsum_wp
	om2wp = al2wp/alsum_wp
	om3wp = al3wp/alsum_wp
	om1wn = al1wn/alsum_wn
	om2wn = al2wn/alsum_wn
	om3wn = al3wn/alsum_wn

	dphidzp(i,j,k) = aux3*(om1wp*(2*wp(i,j,k)-7*wp(i,j,k-1)+11*wp(i,j,k-2)) + &
	om2wp*(-wp(i,j,k-1)+5*wp(i,j,k-2)+2*wp(i,j,k-3)) + om3wp*(2*wp(i,j,k-2)+5*wp(i,j,k-3)-wp(i,j,k-4)))
	dphidzn(i,j,k) = aux3*(om1wn*(2*wn(i,j,k)-7*wn(i,j,k+1)+11*wn(i,j,k+2)) + &
	om2wn*(-wn(i,j,k+1)+5*wn(i,j,k+2)+2*wn(i,j,k+3)) + om3wn*(2*wn(i,j,k+2)+5*wn(i,j,k+3)-wn(i,j,k+4)))
	enddo
	enddo
	enddo
    !$omp end parallel do
		
END SUBROUTINE wenoz
