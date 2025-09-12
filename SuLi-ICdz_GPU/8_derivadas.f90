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

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx !i=2,dimx-1
	    campo_saida(i,j,k)=(a1(i+1,j,k)-a1(i-1,j,k))/(2.*dx)
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

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

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
        !a1=campo_entrada
        do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j+1,k)-a1(i,j-1,k))/(2.*dy)
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE derivay

!##############################################################

SUBROUTINE derivaz(a1,dimx,dimy,dimz,d,campo_saida)
! derivada centrada de segunda ordem na direção z

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1, d) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k+1)-a1(i,j,k-1))/(2.*d(k))
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

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

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx-1
	    campo_saida(i,j,k)=0.5*(-a1(i+2,j,k)+4*a1(i+1,j,k)-3*a1(i,j,k))/dx
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	!i = dimx
	!do k=1,dimz
	! do j=1,dimy
	!    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
	! enddo
	!enddo

	if (ccx0 .ne. 0) then
	 i = dimx
	  do k=1,dimz
	   do j=1,dimy
	    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
	   enddo
	  enddo
	else
            !campo_saida(1,:,:)=campo_saida(dimx,:,:)
	i = dimx
	do k=1,dimz
	 do j=1,dimy
	    campo_saida(i,j,k)=0.5*(-a1(2,j,k)+4*a1(1,j,k)-3*a1(i,j,k))/dx
	 enddo
	enddo
	endif

END SUBROUTINE derivaxu2p

!##############################################################

SUBROUTINE derivayu2p(a1,dimx,dimy,dimz,campo_saida)
! derivada por frontward de segunda ordem na direção y (usado para upwind)

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy-1
	  do i=1,dimx
	    campo_saida(i,j,k)=0.5*(-a1(i,j+2,k)+4*a1(i,j+1,k)-3*a1(i,j,k))/dy
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	j = dimy
!       !$omp parallel do collapse(2) private(i, k)
	do k=1,dimz
	 do i=1,dimx
	    campo_saida(i,j,k)=0.5*(a1(i,j+1,k)-a1(i,j-1,k))/dy
	 enddo
	enddo
!       !$omp end parallel do


END SUBROUTINE derivayu2p

!##############################################################

SUBROUTINE derivazu2p(a1,dimx,dimy,dimz,d,campo_saida)
! derivada por frontward de segunda ordem na direção z (usado para upwind)

!! NOTE QUE MESMO PRA MAIOR ORDEM SE USOU APENAS O DZ DO PONTO CALCULADO

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1, d) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz-1
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=0.5*(-a1(i,j,k+2)+4*a1(i,j,k+1)-3*a1(i,j,k))/d(k)
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

        k = dimz
!       !$omp parallel do collapse(2) private(i, j)
	do j=1,dimy
	 do i=1,dimx
	    campo_saida(i,j,k)=0.5*(a1(i,j,k+1)-a1(i,j,k-1))/d(k)
	 enddo
	enddo
!       !$omp end parallel do



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

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy
	  do i=2,dimx
	    campo_saida(i,j,k)=0.5*(a1(i-2,j,k)-4*a1(i-1,j,k)+3*a1(i,j,k))/dx
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	!i = 1
	!do k=1,dimz
	! do j=1,dimy
	!    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
	! enddo
	!enddo

	if (ccx0 .ne. 0) then
        i = 1
!       !$omp parallel do collapse(2) private(j, k)
	  do k=1,dimz
	   do j=1,dimy
	    campo_saida(i,j,k)=0.5*(a1(i+1,j,k)-a1(i-1,j,k))/dx
	   enddo
	  enddo
!       !$omp end parallel do
	else
	
	i = 1
!       !$omp parallel do collapse(2) private(j, k)
	do k=1,dimz
	 do j=1,dimy
	    campo_saida(i,j,k)=0.5*(a1(dimx-1,j,k)-4*a1(dimx,j,k)+3*a1(i,j,k))/dx
	 enddo
	enddo
!       !$omp end parallel do
	endif

END SUBROUTINE derivaxu2n

!##############################################################

SUBROUTINE derivayu2n(a1,dimx,dimy,dimz,campo_saida)
! derivada por backward de segunda ordem na direção y (usado para upwind)

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=2,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=0.5*(a1(i,j-2,k)-4*a1(i,j-1,k)+3*a1(i,j,k))/dy
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	j = 1
!       !$omp parallel do collapse(2) private(i, k)
	do k=1,dimz
	 do i=1,dimx
	    campo_saida(i,j,k)=0.5*(a1(i,j+1,k)-a1(i,j-1,k))/dy
	 enddo
	enddo
!       !$omp end parallel do

END SUBROUTINE derivayu2n

!##############################################################

SUBROUTINE derivazu2n(a1,dimx,dimy,dimz,d,campo_saida)
! derivada por backward de segunda ordem na direção z (usado para upwind)

!! NOTE QUE MESMO PRA MAIOR ORDEM SE USOU APENAS O DZ DO PONTO CALCULADO

	USE disc
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(0:dimx+1,0:dimy+1,0:dimz+1)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1, d) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=2,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=0.5*(a1(i,j,k-2)-4*a1(i,j,k-1)+3*a1(i,j,k))/d(k)
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	k = 1
!       !$omp parallel do collapse(2) private(i, j)
	do j=1,dimy
	 do i=1,dimx
	    campo_saida(i,j,k)=0.5*(a1(i,j,k+1)-a1(i,j,k-1))/d(k)
	  enddo
	 enddo
!       !$omp end parallel do

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

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
        !a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy
	  do i=2,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i-1,j,k))*0.5
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	if (ccx0 .ne. 0) then

#ifdef USE_GPU
       !$omp target update from(campo_saida)
#endif
        !$omp parallel do private(j,k)
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
	
#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do private(j,k)
#else
        !$omp parallel do private(j,k)
#endif
	do k=1,dimz
	 do j=1,dimy
	   campo_saida(1,j,k)=(a1(1,j,k)+a1(dimx,j,k))*0.5
	   campo_saida(dimx+1,j,k)=(a1(1,j,k)+a1(dimx,j,k))*0.5
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

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

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy
	  do i=1,dimx-1
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i+1,j,k))*0.5
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE interpx_fc

!##############################################################

SUBROUTINE interpy_cf(a1,dimx,dimy,dimz,campo_saida)
! interpolação linear: cf --> centro pra face fc--> face pro centro

	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy+1,dimz) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=2,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j-1,k))*0.5
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	
#ifdef USE_GPU
        !$omp target update from(campo_saida)
#endif
        !$omp parallel do private(i,k)
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

#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz
	 do j=1,dimy-1
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)+a1(i,j+1,k))*0.5
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif


END SUBROUTINE interpy_fc

!##############################################################

SUBROUTINE interpz_cf(a1,dimx,dimy,dimz,d,campo_saida) 
! interpolação linear: cf --> centro pra face fc--> face pro centro

	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz+1) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1, d) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=2,dimz
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=(a1(i,j,k)*d(k)+a1(i,j,k-1)*d(k-1))/(d(k)+d(k-1))
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target update from(campo_saida)
#endif
        !$omp parallel do private(i,j)
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

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(0:dimz+1) :: d  !!!!! uma ideia !!!!!
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(dimx,dimy,dimz-1) :: campo_saida

#ifdef USE_GPU
        !$omp target map(to: a1, d) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	!a1=campo_entrada
	do k=1,dimz-1
	 do j=1,dimy
	  do i=1,dimx
	    campo_saida(i,j,k)=	(a1(i,j,k)*d(k)+a1(i,j,k+1)*d(k+1))/(d(k)+d(k+1))
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE interpz_fc

!##############################################################

SUBROUTINE interpxy_ff(u,nx,ny,nz,nx1,ny1,uy) !ff --> face pra face  !velocidade x na face y

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,nx1,ny,ny1,nz
	real(8),intent(IN),dimension(0:nx1+1,0:ny+1,0:nz+1)  :: u
	real(8),intent(OUT),dimension(nx,ny1,nz) :: uy

#ifdef USE_GPU
        !$omp target map(to: u) map(from: uy)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k=1,nz
	 do j=1,ny1
	  do i=1,nx
		uy(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * 0.25
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

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

#ifdef USE_GPU
        !$omp target map(to: u, d) map(from: uz)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k=1,nz1
	 do j=1,ny
	  do i=1,nx
		uz(i,j,k) = ((u(i,j,k) + u(i+1,j,k))*d(k) + (u(i,j,k-1) + u(i+1,j,k-1))*d(k-1))*0.5/(d(k)+d(k-1))
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE interpxz_ff

!##############################################################

SUBROUTINE interpyx_ff(v,nx,ny,nz,nx1,ny1,vx) !ff --> face pra face  !velocidade y na face x

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: nx,nx1,ny,ny1,nz
	real(8),intent(IN),dimension(0:nx+1,0:ny1+1,0:nz+1)  :: v
	real(8),intent(OUT),dimension(nx1,ny,nz) :: vx

#ifdef USE_GPU
        !$omp target map(to: v) map(from: vx)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k=1,nz
	 do j=1,ny
	  do i=1,nx1
		vx(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * 0.25
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

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

#ifdef USE_GPU
        !$omp target map(to: v, d) map(from: vz)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k=1,nz1
	 do j=1,ny
	  do i=1,nx
		vz(i,j,k) = ((v(i,j,k) + v(i,j+1,k))*d(k) + (v(i,j,k-1) + v(i,j+1,k-1))*d(k-1))*0.5/(d(k)+d(k-1))
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

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

#ifdef USE_GPU
        !$omp target map(to: w, d) map(from: wx)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k=1,nz
	 do j=1,ny
	  do i=1,nx1
		wx(i,j,k) = ((w(i,j,k) + w(i-1,j,k))*d(k) + (w(i,j,k+1)  + w(i-1,j,k+1))*d(k+1))*0.5/(d(k)+d(k+1))
	  enddo
	 enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

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

#ifdef USE_GPU
        !$omp target map(to: w, d) map(from: wy)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k=1,nz
        do j=1,ny1
        do i=1,nx
                wy(i,j,k) = ((w(i,j,k) + w(i,j-1,k))*d(k) + (w(i,j,k+1) + w(i,j-1,k+1))*d(k+1))*0.5/(d(k)+d(k+1))
        enddo
        enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE interpzy_ff

!##############################################################

SUBROUTINE fantasmas_cf(a1,dimx,dimy,dimz,campo_saida)
! interpolação linear: cf --> centro pra face fc--> face pro centro

	USE cond, only : ccx0
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz)  :: a1
	real(8),intent(OUT),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: campo_saida

	!a1=campo_entrada
#ifdef USE_GPU
        !$omp target map(to: a1) map(from: campo_saida)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k=1,dimz
        do j=1,dimy
        do i=1,dimx
                campo_saida(i,j,k)=a1(i,j,k)
        enddo
        enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

       !$omp parallel do collapse(2) private(i, k)
	do k=1,dimz
	  do i=1,dimx
	  j=0
	   campo_saida(i,j,k)=a1(i,j+1,k)
	  j=dimy+1
	   campo_saida(i,j,k)=a1(i,j-1,k)
	 enddo
	enddo
       !$omp end parallel do

       !$omp parallel do collapse(2) private(i, j)
	 do j=1,dimy
	  do i=1,dimx
	  k=0
	   campo_saida(i,j,k)=a1(i,j,k+1)
	  k=dimz+1
	   campo_saida(i,j,k)=a1(i,j,k-1)
	 enddo
	enddo
       !$omp end parallel do

	if (ccx0 == 0) then
	   campo_saida(0,:,:)=campo_saida(dimx,:,:)
	   campo_saida(dimx+1,:,:)=campo_saida(1,:,:)
	else
	   campo_saida(0,:,:)=campo_saida(1,:,:)
	   campo_saida(dimx+1,:,:)=campo_saida(dimx,:,:)
	endif

END SUBROUTINE fantasmas_cf

!##############################################################

SUBROUTINE wenox(a,dimx,dimy,dimz,dx1,dphidxp,dphidxn,ihs)

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i,j,k,kk,ii,ihs,dimx,dimy,dimz
	real(8),intent(in) :: dx1
	real(8),intent(in),dimension(dimx,dimy,dimz) :: a

	real(8),intent(out),dimension(dimx,dimy,dimz) :: dphidxp,dphidxn
	real(8),save,dimension(3) ::isup, isun, auxx
	real(8),save,dimension(3) ::alpup, omgup,alpun, omgun

	real(8),save :: mod_phi1,aux1,aux2,aux3,aux4,aux5,aux6,aux,aux11, aux12

	real(8),dimension(-2:dimx+3,dimy,dimz) :: phi1
	real(8),dimension(dimx+4,dimy,dimz)    :: un
	real(8),dimension(-3:dimx,dimy,dimz)   :: up


	aux1 = 13./12.
	aux2 = 1./4.
	aux3 = 1./6.
	auxx(1) = 0.1
	auxx(2) = 0.6
	auxx(3) = 0.3
	aux6 = 0.00000001

!================================================================
!	phi1(1:dimx,1:dimy,1:dimz) = a(1:dimx,1:dimy,1:dimz)
!================================================================

#ifdef USE_GPU
        !$omp target map(to: a) map(from: phi1)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz
        do j = 1, dimy
        do i = 1, dimx
                phi1(i,j,k) = a(i,j,k)
        enddo
        enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif


	if (ihs == 0) then !Contorno periódico
!               !$omp parallel do collapse(2) private(j, k)
		do k = 1, dimz
		do j = 1, dimy
			phi1(0,j,k)  = phi1(dimx-1,j,k)
			phi1(-1,j,k) = phi1(dimx-2,j,k)
			phi1(-2,j,k) = phi1(dimx-3,j,k)

			phi1(dimx+1,j,k) = phi1(2,j,k)
			phi1(dimx+2,j,k) = phi1(3,j,k)
			phi1(dimx+3,j,k) = phi1(4,j,k)

		enddo
		enddo
!		!$omp end parallel do

	elseif (ihs == 1) then !Distance extrapolation
!               !$omp parallel do collapse(2) private(j, k)
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
!		!$omp end parallel do

	elseif (ihs == 2) then !Derivative zero
!               !$omp parallel do collapse(2) private(j, k)
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
!               !$omp end parallel do
	endif

#ifdef USE_GPU
        !$omp target map(to: phi1) map(from: up)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz
	do j = 1, dimy
	do i =-3, dimx
	        up(i,j,k)=(phi1(i+3,j,k)-phi1(i+2,j,k))/dx1
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target map(to: phi1) map(from: un)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx+4
	        un(i,j,k)=(phi1(i-2,j,k)-phi1(i-3,j,k))/dx1
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif


#ifdef USE_GPU
        !$omp target map(to: un, up, auxx, aux1, aux2, aux3, aux6, dx1) map(from: dphidxp, dphidxn)
        !$omp teams distribute parallel do collapse(3) private(i, j, k, isup, isun, alpup, alpun, omgup, omgun)
#else
        !$omp parallel do collapse(3) private(i, j, k, isup, isun, alpup, alpun, omgup, omgun)
#endif
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
	isup(1) = aux1 * (up(i,j,k)-2*up(i-1,j,k)+up(i-2,j,k))*(up(i,j,k)-2*up(i-1,j,k)+up(i-2,j,k)) &
	+ aux2 * (up(i,j,k)-4*up(i-1,j,k)+3*up(i-2,j,k))*(up(i,j,k)-4*up(i-1,j,k)+3*up(i-2,j,k))
	isun(1) = aux1 * (un(i,j,k)-2*un(i+1,j,k)+un(i+2,j,k))*(un(i,j,k)-2*un(i+1,j,k)+un(i+2,j,k)) &
	+ aux2 * (un(i,j,k)-4*un(i+1,j,k)+3*un(i+2,j,k))*(un(i,j,k)-4*un(i+1,j,k)+3*un(i+2,j,k))

	isup(2) = aux1 * (up(i-1,j,k)-2*up(i-2,j,k)+up(i-3,j,k))*(up(i-1,j,k)-2*up(i-2,j,k)+up(i-3,j,k)) + aux2 * &
	(up(i-1,j,k)-up(i-3,j,k))*(up(i-1,j,k)-up(i-3,j,k))
	isun(2) = aux1 * (un(i+1,j,k)-2*un(i+2,j,k)+un(i+3,j,k))*(un(i+1,j,k)-2*un(i+2,j,k)+un(i+3,j,k)) + aux2 * &
	(un(i+1,j,k)-un(i+3,j,k))*(un(i+1,j,k)-un(i+3,j,k))

	isup(3) = aux1 * (up(i-2,j,k)-2*up(i-3,j,k)+up(i-4,j,k))*(up(i-2,j,k)-2*up(i-3,j,k)+up(i-4,j,k)) &
	+ aux2 * (3*up(i-2,j,k)-4*up(i-3,j,k)+up(i-4,j,k))*(3*up(i-2,j,k)-4*up(i-3,j,k)+up(i-4,j,k))
	isun(3) = aux1 * (un(i+2,j,k)-2*un(i+3,j,k)+un(i+4,j,k))*(un(i+2,j,k)-2*un(i+3,j,k)+un(i+4,j,k)) &
	+ aux2 * (3*un(i+2,j,k)-4*un(i+3,j,k)+un(i+4,j,k))*(3*un(i+2,j,k)-4*un(i+3,j,k)+un(i+4,j,k))

!================================================================
	do kk = 1, 3
	alpup(kk) = auxx(kk) / ((aux6 + isup(kk))*(aux6 + isup(kk)))
	alpun(kk) = auxx(kk) / ((aux6 + isun(kk))*(aux6 + isun(kk)))
	enddo
!================================================================

!	alpup(1) = auxx(1) / ((aux6 + isup(1))*(aux6 + isup(1)))
!	alpun(1) = auxx(1) / ((aux6 + isun(1))*(aux6 + isun(1)))

!	alpup(2) = auxx(2) / ((aux6 + isup(2))*(aux6 + isup(2)))
!	alpun(2) = auxx(2) / ((aux6 + isun(2))*(aux6 + isun(2)))

!	alpup(3) = auxx(3) / ((aux6 + isup(3))*(aux6 + isup(3)))
!	alpun(3) = auxx(3) / ((aux6 + isun(3))*(aux6 + isun(3)))
	
!================================================================
	do kk = 1, 3
	omgup(kk) = alpup(kk) / (alpup(1)+alpup(2)+alpup(3))
	omgun(kk) = alpun(kk) / (alpun(1)+alpun(2)+alpun(3))
	enddo
!================================================================

!	omgup(1) = alpup(1) / (alpup(1)+alpup(2)+alpup(3))
!	omgun(1) = alpun(1) / (alpun(1)+alpun(2)+alpun(3))

!	omgup(2) = alpup(2) / (alpup(1)+alpup(2)+alpup(3))
!	omgun(2) = alpun(2) / (alpun(1)+alpun(2)+alpun(3))

!	omgup(3) = alpup(3) / (alpup(1)+alpup(2)+alpup(3))
!	omgun(3) = alpun(3) / (alpun(1)+alpun(2)+alpun(3))
	
	dphidxp(i,j,k) = aux3* (omgup(1) * (2*up(i,j,k)-7*up(i-1,j,k)+11*up(i-2,j,k)) + &
	omgup(2) * (-up(i-1,j,k)+5*up(i-2,j,k)+2*up(i-3,j,k)) + omgup(3) * (2*up(i-2,j,k)+5*up(i-3,j,k)-up(i-4,j,k)))
	dphidxn(i,j,k) = aux3* (omgun(1) * (2*un(i,j,k)-7*un(i+1,j,k)+11*un(i+2,j,k)) + &
	omgun(2) * (-un(i+1,j,k)+5*un(i+2,j,k)+2*un(i+3,j,k)) + omgun(3) * (2*un(i+2,j,k)+5*un(i+3,j,k)-un(i+4,j,k)))
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE wenox

!##############################################################

SUBROUTINE wenoy(b,dimx,dimy,dimz,dy1,dphidyp,dphidyn,ihs)

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i,j,k,kk,ii,ihs,dimx,dimy,dimz
	real(8),intent(in) :: dy1
	real(8),intent(in),dimension(dimx,dimy,dimz) :: b
				
	real(8),intent(out),dimension(dimx,dimy,dimz) :: dphidyp,dphidyn
	real(8),save,dimension(3) ::isvp, isvn, auxx
	real(8),save,dimension(3) ::alpvp, omgvp,alpvn, omgvn

	real(8),save :: mod_phi1,aux1,aux2,aux3,aux4,aux5,aux6,aux,aux11,aux12
	
	real(8),dimension(dimx,-2:dimy+3,dimz) :: phi1
	real(8),dimension(dimx,dimy+4,dimz)    :: vn
	real(8),dimension(dimx,-3:dimy,dimz)   :: vp
	
	
	aux1 = 13./12.
	aux2 = 1./4.
	aux3 = 1./6.
	auxx(1) = 0.1
	auxx(2) = 0.6
	auxx(3) = 0.3
	aux6 = 0.00000001
	
!================================================================
!	phi1(1:dimx,1:dimy,1:dimz) = b(1:dimx,1:dimy,1:dimz)
!================================================================

#ifdef USE_GPU
        !$omp target map(to: b) map(from: phi1)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz
        do j = 1, dimy
        do i = 1, dimx
                phi1(i,j,k) = b(i,j,k)
        enddo
        enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif
	

	if (ihs == 0) then !Contorno periódico
!               !$omp parallel do collapse(2) private(i, k)
		do k = 1, dimz
		do i = 1, dimx
			phi1(i,0,k)  = phi1(i,dimy-1,k)
			phi1(i,-1,k) = phi1(i,dimy-2,k)
			phi1(i,-2,k) = phi1(i,dimy-3,k)

			phi1(i,dimy+1,k) = phi1(i,2,k)
			phi1(i,dimy+2,k) = phi1(i,3,k)
			phi1(i,dimy+3,k) = phi1(i,4,k)
			
		enddo
		enddo
!               !$omp end parallel do
		
	elseif (ihs == 1) then !Distance extrapolation
!               !$omp parallel do collapse(2) private(i, k)
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
!		!$omp end parallel do
		
	elseif (ihs == 2) then !Derivative zero
!               !$omp parallel do collapse(2) private(i, k)
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
!               !$omp end parallel do
	endif
    
#ifdef USE_GPU
        !$omp target map(to: phi1) map(from: vp)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz
	do j = -3, dimy
	do i = 1, dimx
                vp(i,j,k)=(phi1(i,j+3,k)-phi1(i,j+2,k))/dy1
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target map(to: phi1) map(from: vn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz
	do j = 1, dimy+4
	do i = 1, dimx
                vn(i,j,k)=(phi1(i,j-2,k)-phi1(i,j-3,k))/dy1
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target map(to: auxx, aux1, aux2, aux3, aux6, dy1) map(from: dphidyp, dphidyn)
        !$omp teams distribute parallel do collapse(3) private(i, j, k, isvp, isvn, alpvp, alpvn, omgvp, omgvn)
#else
        !$omp parallel do collapse(3) private(i, j, k, isvp, isvn, alpvp, alpvn, omgvp, omgvn)
#endif
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
	isvp(1) = aux1 * (vp(i,j,k)-2*vp(i,j-1,k)+vp(i,j-2,k))*(vp(i,j,k)-2*vp(i,j-1,k)+vp(i,j-2,k)) &
	+ aux2 * (vp(i,j,k)-4*vp(i,j-1,k)+3*vp(i,j-2,k))*(vp(i,j,k)-4*vp(i,j-1,k)+3*vp(i,j-2,k))
	isvn(1) = aux1 * (vn(i,j,k)-2*vn(i,j+1,k)+vn(i,j+2,k))*(vn(i,j,k)-2*vn(i,j+1,k)+vn(i,j+2,k)) &
	+ aux2 * (vn(i,j,k)-4*vn(i,j+1,k)+3*vn(i,j+2,k))*(vn(i,j,k)-4*vn(i,j+1,k)+3*vn(i,j+2,k))

	isvp(2) = aux1 * (vp(i,j-1,k)-2*vp(i,j-2,k)+vp(i,j-3,k))*(vp(i,j-1,k)-2*vp(i,j-2,k)+vp(i,j-3,k)) + aux2 * &
	(vp(i,j-1,k)-vp(i,j-3,k))*(vp(i,j-1,k)-vp(i,j-3,k))
	isvn(2) = aux1 * (vn(i,j+1,k)-2*vn(i,j+2,k)+vn(i,j+3,k))*(vn(i,j+1,k)-2*vn(i,j+2,k)+vn(i,j+3,k)) + aux2 * &
	(vn(i,j+1,k)-vn(i,j+3,k))*(vn(i,j+1,k)-vn(i,j+3,k))

	isvp(3) = aux1 * (vp(i,j-2,k)-2*vp(i,j-3,k)+vp(i,j-4,k))*(vp(i,j-2,k)-2*vp(i,j-3,k)+vp(i,j-4,k)) &
	+ aux2 * (3*vp(i,j-2,k)-4*vp(i,j-3,k)+vp(i,j-4,k))*(3*vp(i,j-2,k)-4*vp(i,j-3,k)+vp(i,j-4,k))
	isvn(3) = aux1 * (vn(i,j+2,k)-2*vn(i,j+3,k)+vn(i,j+4,k))*(vn(i,j+2,k)-2*vn(i,j+3,k)+vn(i,j+4,k)) &
	+ aux2 * (3*vn(i,j+2,k)-4*vn(i,j+3,k)+vn(i,j+4,k))*(3*vn(i,j+2,k)-4*vn(i,j+3,k)+vn(i,j+4,k))
	
!================================================================
	do kk = 1, 3
	alpvp(kk) = auxx(kk) / ((aux6 + isvp(kk))*(aux6 + isvp(kk)))
	alpvn(kk) = auxx(kk) / ((aux6 + isvn(kk))*(aux6 + isvn(kk)))
	enddo
!================================================================

!	alpvp(1) = auxx(1) / ((aux6 + isvp(1))*(aux6 + isvp(1)))
!	alpvn(1) = auxx(1) / ((aux6 + isvn(1))*(aux6 + isvn(1)))
	
!	alpvp(2) = auxx(2) / ((aux6 + isvp(2))*(aux6 + isvp(2)))
!	alpvn(2) = auxx(2) / ((aux6 + isvn(2))*(aux6 + isvn(2)))
	
!	alpvp(3) = auxx(3) / ((aux6 + isvp(3))*(aux6 + isvp(3)))
!	alpvn(3) = auxx(3) / ((aux6 + isvn(3))*(aux6 + isvn(3)))
	
!================================================================
	do kk = 1, 3
	omgvp(kk) = alpvp(kk) / (alpvp(1)+alpvp(2)+alpvp(3))
	omgvn(kk) = alpvn(kk) / (alpvn(1)+alpvn(2)+alpvn(3))
	enddo
!================================================================

!	omgvp(1) = alpvp(1) / (alpvp(1)+alpvp(2)+alpvp(3))
!	omgvn(1) = alpvn(1) / (alpvn(1)+alpvn(2)+alpvn(3))
	
!	omgvp(2) = alpvp(2) / (alpvp(1)+alpvp(2)+alpvp(3))
!	omgvn(2) = alpvn(2) / (alpvn(1)+alpvn(2)+alpvn(3))
	
!	omgvp(3) = alpvp(3) / (alpvp(1)+alpvp(2)+alpvp(3))
!	omgvn(3) = alpvn(3) / (alpvn(1)+alpvn(2)+alpvn(3))

	dphidyp(i,j,k) = aux3* (omgvp(1) * (2*vp(i,j,k)-7*vp(i,j-1,k)+11*vp(i,j-2,k)) + &
	omgvp(2) * (-vp(i,j-1,k)+5*vp(i,j-2,k)+2*vp(i,j-3,k)) + omgvp(3) * (2*vp(i,j-2,k)+5*vp(i,j-3,k)-vp(i,j-4,k)))
	dphidyn(i,j,k) = aux3* (omgvn(1) * (2*vn(i,j,k)-7*vn(i,j+1,k)+11*vn(i,j+2,k)) + &
	omgvn(2) * (-vn(i,j+1,k)+5*vn(i,j+2,k)+2*vn(i,j+3,k)) + omgvn(3) * (2*vn(i,j+2,k)+5*vn(i,j+3,k)-vn(i,j+4,k)))
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE wenoy

!##############################################################

SUBROUTINE wenoz(d,dimx,dimy,dimz,dz1,dphidzp,dphidzn,ihs)

        !$ USE omp_lib

	IMPLICIT NONE
	integer :: i,j,k,kk,ii,ihs,dimx,dimy,dimz
	real(8),intent(in),dimension(0:dimz+1) :: dz1 !ter certeza desses limites!
	real(8),intent(in),dimension(dimx,dimy,dimz) :: d
				
	real(8),intent(out),dimension(dimx,dimy,dimz) :: dphidzp,dphidzn
	real(8),save,dimension(3) ::iswp, iswn, auxx
	real(8),save,dimension(3) ::alpwp, omgwp,alpwn, omgwn

	real(8),save :: mod_phi1,aux1,aux2,aux3,aux4,aux5,aux6,aux,aux11,aux12

	real(8),dimension(-2:dimz+3) :: dz11	
	real(8),dimension(dimx,dimy,-2:dimz+3) :: phi1
	real(8),dimension(dimx,dimy,dimz+4)    :: wn
	real(8),dimension(dimx,dimy,-3:dimz)   :: wp
	
	
	aux1 = 13./12.
	aux2 = 1./4.
	aux3 = 1./6.
	auxx(1) = 0.1
	auxx(2) = 0.6
	auxx(3) = 0.3
	aux6 = 0.00000001

!================================================================
!	phi1(1:dimx,1:dimy,1:dimz) = d(1:dimx,1:dimy,1:dimz)
!================================================================

#ifdef USE_GPU
        !$omp target map(to: d) map(from: phi1)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
	        phi1(i,j,k) = d(i,j,k)
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	dz11(0:dimz+1) = dz1(0:dimz+1)
	
	
	
	if (ihs == 0) then !Contorno periódico
!               !$omp parallel do collapse(2) private(i, j)
		do j = 1, dimy
		do i = 1, dimx
			phi1(i,j,0)  = phi1(i,j,dimz-1)
			phi1(i,j,-1) = phi1(i,j,dimz-2)
			phi1(i,j,-2) = phi1(i,j,dimz-3)

			phi1(i,j,dimz+1) = phi1(i,j,2)
			phi1(i,j,dimz+2) = phi1(i,j,3)
			phi1(i,j,dimz+3) = phi1(i,j,4)
			
			
			dz11(0)  = dz11(dimz-1)
			dz11(-1) = dz11(dimz-2)
			dz11(-2) = dz11(dimz-3)

			dz11(dimz+1) = dz11(2)
			dz11(dimz+2) = dz11(3)
			dz11(dimz+3) = dz11(4)		
			
		enddo
		enddo
!               !$omp end parallel do
		
	elseif (ihs == 1) then !Velocidade zero com derivada de 1a ordem 
	! no-slip de velocidade = 0 ou free-slip para velocidade perpendicular
	        !!$omp parallel do collapse(2) private(i, j)
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
				

			dz11(0)  = 1./5. * (12.*dz11(1)  - 9.*dz11(2) + 2.*dz11(3) )
			dz11(-1) = 1./5. * (12.*dz11(0)  - 9.*dz11(1) + 2.*dz11(2) )
			dz11(-2) = 1./5. * (12.*dz11(-1) - 9.*dz11(0) + 2.*dz11(1) )

			dz11(dimz+1) = 1./11. * (18.*dz11(dimz)   - 9.*dz11(dimz-1) + 2.*dz11(dimz-2))
			dz11(dimz+2) = 1./11. * (18.*dz11(dimz+1) - 9.*dz11(dimz)   + 2.*dz11(dimz-1))
			dz11(dimz+3) = 1./11. * (18.*dz11(dimz+2) - 9.*dz11(dimz+1) + 2.*dz11(dimz)  )			
				
		enddo
		enddo
		!!$omp end parallel do
		
	elseif (ihs == 2) then !Derivative zero
!               !$omp parallel do collapse(2) private(i, j)
		do j = 1, dimy
		do i = 1, dimx
			
			phi1(i,j,0)  = 1./11. * (18.*phi1(i,j,1)  - 9.*phi1(i,j,2) + 2.*phi1(i,j,3) )
			phi1(i,j,-1) = 1./11. * (18.*phi1(i,j,0)  - 9.*phi1(i,j,1) + 2.*phi1(i,j,2) )
			phi1(i,j,-2) = 1./11. * (18.*phi1(i,j,-1) - 9.*phi1(i,j,0) + 2.*phi1(i,j,1) )

			phi1(i,j,dimz+1) = 1./11. * (18.*phi1(i,j,dimz)   - 9.*phi1(i,j,dimz-1) + 2.*phi1(i,j,dimz-2))
			phi1(i,j,dimz+2) = 1./11. * (18.*phi1(i,j,dimz+1) - 9.*phi1(i,j,dimz)   + 2.*phi1(i,j,dimz-1))
			phi1(i,j,dimz+3) = 1./11. * (18.*phi1(i,j,dimz+2) - 9.*phi1(i,j,dimz+1) + 2.*phi1(i,j,dimz)  )
			

			dz11(0)  = 1./11. * (18.*dz11(1)  - 9.*dz11(2) + 2.*dz11(3) )
			dz11(-1) = 1./11. * (18.*dz11(0)  - 9.*dz11(1) + 2.*dz11(2) )
			dz11(-2) = 1./11. * (18.*dz11(-1) - 9.*dz11(0) + 2.*dz11(1) )

			dz11(dimz+1) = 1./11. * (18.*dz11(dimz)   - 9.*dz11(dimz-1) + 2.*dz11(dimz-2))
			dz11(dimz+2) = 1./11. * (18.*dz11(dimz+1) - 9.*dz11(dimz)   + 2.*dz11(dimz-1))
			dz11(dimz+3) = 1./11. * (18.*dz11(dimz+2) - 9.*dz11(dimz+1) + 2.*dz11(dimz)  )
					
		enddo
		enddo
!               !$omp end parallel do
	endif
    
#ifdef USE_GPU
        !$omp target map(to: phi1, dz11) map(from: wp)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = -3, dimz
	do j = 1, dimy
	do i = 1, dimx
	wp(i,j,k)=(phi1(i,j,k+3)-phi1(i,j,k+2))/dz11(k+3)
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target map(to: phi1, dz11) map(from: wn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, dimz+4
	do j = 1, dimy
	do i = 1, dimx
	wn(i,j,k)=(phi1(i,j,k-2)-phi1(i,j,k-3))/dz11(k-2)
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

#ifdef USE_GPU
        !$omp target map(to: auxx, aux1, aux2, aux3, aux6, dz11) map(from: dphidzp, dphidzn)
        !$omp teams distribute parallel do collapse(3) private(i, j, k, iswp, iswn, alpwp, alpwn, omgwp, omgwn)
#else
        !$omp parallel do collapse(3) private(i, j, k, iswp, iswn, alpwp, alpwn, omgwp, omgwn)
#endif
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
	iswp(1) = aux1 * (wp(i,j,k)-2*wp(i,j,k-1)+wp(i,j,k-2))*(wp(i,j,k)-2*wp(i,j,k-1)+wp(i,j,k-2)) &
	+ aux2 * (wp(i,j,k)-4*wp(i,j,k-1)+3*wp(i,j,k-2))*(wp(i,j,k)-4*wp(i,j,k-1)+3*wp(i,j,k-2))
	iswn(1) = aux1 * (wn(i,j,k)-2*wn(i,j,k+1)+wn(i,j,k+2))*(wn(i,j,k)-2*wn(i,j,k+1)+wn(i,j,k+2)) &
	+ aux2 * (wn(i,j,k)-4*wn(i,j,k+1)+3*wn(i,j,k+2))*(wn(i,j,k)-4*wn(i,j,k+1)+3*wn(i,j,k+2))

	iswp(2) = aux1 * (wp(i,j,k-1)-2*wp(i,j,k-2)+wp(i,j,k-3))*(wp(i,j,k-1)-2*wp(i,j,k-2)+wp(i,j,k-3)) + aux2 * &
	(wp(i,j,k-1)-wp(i,j,k-3))*(wp(i,j,k-1)-wp(i,j,k-3))
	iswn(2) = aux1 * (wn(i,j,k+1)-2*wn(i,j,k+2)+wn(i,j,k+3))*(wn(i,j,k+1)-2*wn(i,j,k+2)+wn(i,j,k+3)) + aux2 * &
	(wn(i,j,k+1)-wn(i,j,k+3))*(wn(i,j,k+1)-wn(i,j,k+3))

	iswp(3) = aux1 * (wp(i,j,k-2)-2*wp(i,j,k-3)+wp(i,j,k-4))*(wp(i,j,k-2)-2*wp(i,j,k-3)+wp(i,j,k-4)) &
	+ aux2 * (3*wp(i,j,k-2)-4*wp(i,j,k-3)+wp(i,j,k-4))*(3*wp(i,j,k-2)-4*wp(i,j,k-3)+wp(i,j,k-4))
	iswn(3) = aux1 * (wn(i,j,k+2)-2*wn(i,j,k+3)+wn(i,j,k+4))*(wn(i,j,k+2)-2*wn(i,j,k+3)+wn(i,j,k+4)) &
	+ aux2 * (3*wn(i,j,k+2)-4*wn(i,j,k+3)+wn(i,j,k+4))*(3*wn(i,j,k+2)-4*wn(i,j,k+3)+wn(i,j,k+4))

!================================================================
	do kk = 1, 3
	alpwp(kk) = auxx(kk) / ((aux6 + iswp(kk))*(aux6 + iswp(kk)))
	alpwn(kk) = auxx(kk) / ((aux6 + iswn(kk))*(aux6 + iswn(kk)))
	enddo
!================================================================

!	alpwp(1) = auxx(1) / ((aux6 + iswp(1))*(aux6 + iswp(1)))
!	alpwn(1) = auxx(1) / ((aux6 + iswn(1))*(aux6 + iswn(1)))
	
!	alpwp(2) = auxx(2) / ((aux6 + iswp(2))*(aux6 + iswp(2)))
!	alpwn(2) = auxx(2) / ((aux6 + iswn(2))*(aux6 + iswn(2)))
	
!	alpwp(3) = auxx(3) / ((aux6 + iswp(3))*(aux6 + iswp(3)))
!	alpwn(3) = auxx(3) / ((aux6 + iswn(3))*(aux6 + iswn(3)))

!================================================================
	do kk = 1, 3
	omgwp(kk) = alpwp(kk) / (alpwp(1)+alpwp(2)+alpwp(3))
	omgwn(kk) = alpwn(kk) / (alpwn(1)+alpwn(2)+alpwn(3))
	enddo
!================================================================

!	omgwp(1) = alpwp(1) / (alpwp(1)+alpwp(2)+alpwp(3))
!	omgwn(1) = alpwn(1) / (alpwn(1)+alpwn(2)+alpwn(3))

!	omgwp(2) = alpwp(2) / (alpwp(1)+alpwp(2)+alpwp(3))
!	omgwn(2) = alpwn(2) / (alpwn(1)+alpwn(2)+alpwn(3))

!	omgwp(3) = alpwp(3) / (alpwp(1)+alpwp(2)+alpwp(3))
!	omgwn(3) = alpwn(3) / (alpwn(1)+alpwn(2)+alpwn(3))

	dphidzp(i,j,k) = aux3* (omgwp(1) * (2*wp(i,j,k)-7*wp(i,j,k-1)+11*wp(i,j,k-2)) + &
	omgwp(2) * (-wp(i,j,k-1)+5*wp(i,j,k-2)+2*wp(i,j,k-3)) + omgwp(3) * (2*wp(i,j,k-2)+5*wp(i,j,k-3)-wp(i,j,k-4)))
	dphidzn(i,j,k) = aux3* (omgwn(1) * (2*wn(i,j,k)-7*wn(i,j,k+1)+11*wn(i,j,k+2)) + &
	omgwn(2) * (-wn(i,j,k+1)+5*wn(i,j,k+2)+2*wn(i,j,k+3)) + omgwn(3) * (2*wn(i,j,k+2)+5*wn(i,j,k+3)-wn(i,j,k+4)))
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif
	
END SUBROUTINE wenoz
