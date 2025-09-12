

SUBROUTINE classico(uint,vint,wint)

	USE velpre
	USE param
	!$ USE omp_lib

	IMPLICIT NONE

	!===================================================================================================================
	!contadores
	integer :: i, j, k
	
	real(8),intent(out),dimension(nx1,ny,nz) :: uint
	real(8),intent(out),dimension(nx,ny1,nz) :: vint
	real(8),intent(out),dimension(nx,ny,nz1) :: wint

	real(8),dimension(nx1,ny,nz) :: dudx, dudy, dudz, bma, dma, dudxa, dudya, dudza
	real(8),dimension(nx,ny1,nz) :: dvdx, dvdy, dvdz, amb, dmb, dvdxa, dvdya, dvdza
	real(8),dimension(nx,ny,nz1) :: dwdx, dwdy, dwdz, amd, bmd, dwdxa, dwdya, dwdza
	real(8),dimension(nx,ny,nz)  :: aux

	!auxiliares
	real(8),save :: aux1, aux2, aux3, aux4, aux5, aux6, zero, one_quarter

	zero = 0.00000000
	one_quarter = 0.25000000

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then
	! upwind

#ifdef USE_GPU
        !$omp target map(to: u, v, w, dz, zero, one_quarter) map(from: uint, vint, wint)
#endif
        !$omp parallel private(i, j, k, aux1, aux2, aux3, aux4, aux5, aux6)
        
        !$omp do collapse(2)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * one_quarter
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * one_quarter

	aux1 = max(u(i,j,k),zero)
	aux2 = min(u(i,j,k),zero)
	dudx(i,j,k) = aux1*(u(i,j,k)-u(i-1,j,k))/dx + aux2*(u(i+1,j,k)-u(i,j,k))/dx

	aux1 = max(bma(i,j,k),zero)
	aux2 = min(bma(i,j,k),zero)
	dudy(i,j,k) = aux1*(u(i,j,k)-u(i,j-1,k))/dy + aux2*(u(i,j+1,k)-u(i,j,k))/dy

	aux1 = max(dma(i,j,k),zero)
	aux2 = min(dma(i,j,k),zero)
	dudz(i,j,k) = aux1*(u(i,j,k)-u(i,j,k-1))/dz(k) + aux2*(u(i,j,k+1)-u(i,j,k))/dz(k)

	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

	enddo
	enddo
	enddo
        !$omp end do nowait

        !$omp do collapse(2)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * one_quarter
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * one_quarter

	aux3 = max(amb(i,j,k),zero)
	aux4 = min(amb(i,j,k),zero)
	dvdx(i,j,k) = aux3*(v(i,j,k)-v(i-1,j,k))/dx + aux4*(v(i+1,j,k)-v(i,j,k))/dx

	aux3 = max(v(i,j,k),zero)
	aux4 = min(v(i,j,k),zero)
	dvdy(i,j,k) = aux3*(v(i,j,k)-v(i,j-1,k))/dy + aux4*(v(i,j+1,k)-v(i,j,k))/dy

	aux3 = max(dmb(i,j,k),zero)
	aux4 = min(dmb(i,j,k),zero)
	dvdz(i,j,k) = aux3*(v(i,j,k)-v(i,j,k-1))/dz(k) + aux4*(v(i,j,k+1)-v(i,j,k))/dz(k)

	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

	enddo
	enddo
	enddo
	!$omp end do nowait

        !$omp do collapse(2)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * one_quarter
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * one_quarter

	aux5 = max(amd(i,j,k),zero)
	aux6 = min(amd(i,j,k),zero)
	dwdx(i,j,k) = aux5*(w(i,j,k)-w(i-1,j,k))/dx + aux6*(w(i+1,j,k)-w(i,j,k))/dx

	aux5 = max(bmd(i,j,k),zero)
	aux6 = min(bmd(i,j,k),zero)
	dwdy(i,j,k) = aux5*(w(i,j,k)-w(i,j-1,k))/dy + aux6*(w(i,j+1,k)-w(i,j,k))/dy

	aux5 = max(w(i,j,k),zero)
	aux6 = min(w(i,j,k),zero)
	dwdz(i,j,k) = aux5*(w(i,j,k)-w(i,j,k-1))/dzz(k) + aux6*(w(i,j,k+1)-w(i,j,k))/dzz(k)

	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)

	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif


	elseif (der == 2) then
	! centrado de segunda ordem

	CALL derivax(u,nx1,ny,nz,dudx)
	CALL derivay(u,nx1,ny,nz,dudy)
	CALL derivaz(u,nx1,ny,nz,dz,dudz)

	CALL derivax(v,nx,ny1,nz,dvdx)
	CALL derivay(v,nx,ny1,nz,dvdy)
	CALL derivaz(v,nx,ny1,nz,dz,dvdz)

	CALL derivax(w,nx,ny,nz1,dwdx)
	CALL derivay(w,nx,ny,nz1,dwdy)
	CALL derivaz(w,nx,ny,nz1,dzz,dwdz)

	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,bma) !(nx1,ny,nz)

	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,dzz,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,dma) !(nx1,ny,nz)

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,amb) !(nx,ny1,nz)

	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,dzz,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,dmb) !(nx,ny1,nz)

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,dz,amd) !(nx,ny,nz1)

	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,dz,bmd) !(nx,ny,nz1)

#ifdef USE_GPU
        !$omp target map(to: u, v, w, dudx, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, bma, dma, amb, dmb, amd, bmd) map(from: uint, vint, wint)
#endif
        !$omp parallel
        
        !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
	uint(i,j,k) = u(i,j,k)*dudx(i,j,k) + bma(i,j,k)*dudy(i,j,k) + dma(i,j,k)*dudz(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do nowait

        !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
	vint(i,j,k) = amb(i,j,k)*dvdx(i,j,k) + v(i,j,k)*dvdy(i,j,k) + dmb(i,j,k)*dvdz(i,j,k)
	enddo
	enddo
	enddo
        !$omp end do nowait

        !$omp do collapse(3)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
	wint(i,j,k) = amd(i,j,k)*dwdx(i,j,k) + bmd(i,j,k)*dwdy(i,j,k) + w(i,j,k)*dwdz(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif

	elseif (der == 3) then
	! upwind 2nd order
	
!====================================================
!       uint(1:nx1,1:ny,1:nz) = u(1:nx1,1:ny,1:nz)
!====================================================

#ifdef USE_GPU
        !$omp target map(to: u) map(from: uint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx1
                uint(i,j,k) = u(i,j,k)
        enddo
        enddo
        enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	call interpx_fc(uint,nx1,ny,nz,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,amb) !(nx,ny1,nz)
	call interpx_fc(uint,nx1,ny,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,dz,amd) !(nx,ny,nz1)
	
!====================================================
!	vint(1:nx,1:ny1,1:nz) = v(1:nx,1:ny1,1:nz)
!====================================================

#ifdef USE_GPU
        !$omp target map(to: v) map(from: vint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
	        vint(i,j,k) = v(i,j,k)
	enddo
	enddo
	enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	call interpy_fc(vint,nx,ny1,nz,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,bma) !(nx1,ny,nz)
	call interpy_fc(vint,nx,ny1,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,dz,bmd) !(nx,ny,nz1)
	
!====================================================
!	wint(1:nx,1:ny,1:nz1) = w(1:nx,1:ny,1:nz1)
!====================================================

#ifdef USE_GPU
        !$omp target map(to: w) map(from: wint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
        do k = 1, nz1
        do j = 1, ny
        do i = 1, nx
	        wint(i,j,k) = w(i,j,k)
        enddo
        enddo
        enddo
#ifdef USE_GPU
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	call interpz_fc(wint,nx,ny,nz1,dzz,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,dma) !(nx1,ny,nz)
	call interpz_fc(wint,nx,ny,nz1,dzz,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,dmb) !(nx,ny1,nz)

	CALL derivaxu2n(u,nx1,ny,nz,dudxa)
	CALL derivaxu2p(u,nx1,ny,nz,dudx)
	CALL derivayu2n(u,nx1,ny,nz,dudya)
	CALL derivayu2p(u,nx1,ny,nz,dudy)
	CALL derivazu2n(u,nx1,ny,nz,dz,dudza)
	CALL derivazu2p(u,nx1,ny,nz,dz,dudz)

	CALL derivaxu2n(v,nx,ny1,nz,dvdxa)
	CALL derivaxu2p(v,nx,ny1,nz,dvdx)
	CALL derivayu2n(v,nx,ny1,nz,dvdya)
	CALL derivayu2p(v,nx,ny1,nz,dvdy)
	CALL derivazu2n(v,nx,ny1,nz,dz,dvdza)
	CALL derivazu2p(v,nx,ny1,nz,dz,dvdz)

	CALL derivaxu2n(w,nx,ny,nz1,dwdxa)
	CALL derivaxu2p(w,nx,ny,nz1,dwdx)
	CALL derivayu2n(w,nx,ny,nz1,dwdya)
	CALL derivayu2p(w,nx,ny,nz1,dwdy)
	CALL derivazu2n(w,nx,ny,nz1,dzz,dwdza)
	CALL derivazu2p(w,nx,ny,nz1,dzz,dwdz)

#ifdef USE_GPU
        !$omp target map(to: dudxa, dudya, dudza, dvdxa, dvdya, dvdza, dwdxa, dwdya, dwdza, dudx, dudz, dvdx, dvdy, dvdz, dwdx, dwdy, dwdz, zero) map(from: uint, vint, wint)
#endif
        !$omp parallel private(i, j, k, aux1, aux2, aux3, aux4, aux5, aux6)

        !$omp do collapse(2)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	aux1 = max(u(i,j,k),zero)
	aux2 = min(u(i,j,k),zero)
	dudx(i,j,k) = aux1*dudxa(i,j,k) + aux2*dudx(i,j,k)

	aux1 = max(bma(i,j,k),zero)
	aux2 = min(bma(i,j,k),zero)
	dudy(i,j,k) = aux1*dudya(i,j,k) + aux2*dudy(i,j,k)

	aux1 = max(dma(i,j,k),zero)
	aux2 = min(dma(i,j,k),zero)
	dudz(i,j,k) = aux1*dudza(i,j,k) + aux2*dudz(i,j,k)

	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

	enddo
	enddo
	enddo
	!$omp end do nowait

        !$omp do collapse(2)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	aux3 = max(amb(i,j,k),zero)
	aux4 = min(amb(i,j,k),zero)
	dvdx(i,j,k) = aux3*dvdxa(i,j,k) + aux4*dvdx(i,j,k)

	aux3 = max(v(i,j,k),zero)
	aux4 = min(v(i,j,k),zero)
	dvdy(i,j,k) = aux3*dvdya(i,j,k) + aux4*dvdy(i,j,k)

	aux3 = max(dmb(i,j,k),zero)
	aux4 = min(dmb(i,j,k),zero)
	dvdz(i,j,k) = aux3*dvdza(i,j,k) + aux4*dvdz(i,j,k)

	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

	enddo
	enddo
	enddo
	!$omp end do nowait

        !$omp do collapse(2)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	aux5 = max(amd(i,j,k),zero)
	aux6 = min(amd(i,j,k),zero)
	dwdx(i,j,k) = aux5*dwdxa(i,j,k) + aux6*dwdx(i,j,k)

	aux5 = max(bmd(i,j,k),zero)
	aux6 = min(bmd(i,j,k),zero)
	dwdy(i,j,k) = aux5*dwdya(i,j,k) + aux6*dwdy(i,j,k)

	aux5 = max(w(i,j,k),zero)
	aux6 = min(w(i,j,k),zero)
	dwdz(i,j,k) = aux5*dwdza(i,j,k) + aux6*dwdz(i,j,k)

	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)

	enddo
	enddo
	enddo
        !$omp end do
        !$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif

	endif

!==================================================================================================================
END SUBROUTINE classico


SUBROUTINE rotacional(uint,vint,wint)

	USE velpre
	USE param
	!$ USE omp_lib

	IMPLICIT NONE

	!===================================================================================================================
	real(8),intent(out),dimension(nx1,ny,nz) :: uint
	real(8),intent(out),dimension(nx,ny1,nz) :: vint
	real(8),intent(out),dimension(nx,ny,nz1) :: wint

	real(8),dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: ap, an, bma, dma
	real(8),dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: bp, bn, amb, dmb
	real(8),dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: dp, dn, amd, bmd

	!
	real(8),dimension(nx1,ny1,nz1) :: dudx, dvdx, dwdx
	real(8),dimension(nx1,ny1,nz1) :: dudy, dvdy, dwdy
	real(8),dimension(nx1,ny1,nz1) :: dudz, dvdz, dwdz

	!
	real(8),save :: aa, bb, dd

	!contadores
	integer :: i, j, k, ai, bi, di

	!plotagem
	real(8),save :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8),save :: aux1, aux2, zero, one_quarter
	
	zero = 0.00000000
	one_quarter = 0.25000000


	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then
	! upwind
        
#ifdef USE_GPU
        !$omp target map(to: u, v, w, zero, one_quarter) map(from: ap, an, bp, bn, dp, dn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 0, nz+1
	do j = 0, ny+1
	do i = 1, nx+1

	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * one_quarter
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * one_quarter

	ap(i,j,k) = max(u(i,j,k),zero)
	an(i,j,k) = min(u(i,j,k),zero)

	bp(i,j,k) = max(bma(i,j,k),zero)
	bn(i,j,k) = min(bma(i,j,k),zero)

	dp(i,j,k) = max(dma(i,j,k),zero)
	dn(i,j,k) = min(dma(i,j,k),zero)

	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

!====================================================================
!	ap(0,0:ny+1,0:nz+1)     = max(u(0,0:ny+1,0:nz+1),0.)
!	an(nx1+1,0:ny+1,0:nz+1) = min(u(nx1+1,0:ny+1,0:nz+1),0.)
!====================================================================

#ifdef USE_GPU
        !$omp target map(to: u, zero) map(from: ap, an)
#endif
        !$omp parallel private(j, k)

        !$omp do collapse(2) private(j, k)
        do k = 0, nz+1
        do j = 0, ny+1
            ap(0, j, k) = max(u(0, j, k), zero)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(j, k)
        do k = 0, nz+1
        do j = 0, ny+1
            an(nx1+1, j, k) = min(u(nx1+1, j, k), zero)
        enddo
        enddo
        !$omp end do
        !$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif


#ifdef USE_GPU
        !$omp target map(to: ap, an, bp, bp, bn, dp, dn, u, dx, dy, dz) map(from: uint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k,dudx,dudy,dudz)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	dudx(i,j,k) = (ap(i,j,k)*u(i,j,k)-ap(i-1,j,k)*u(i-1,j,k))/dx + (an(i+1,j,k)*u(i+1,j,k)-an(i,j,k)*u(i,j,k))/dx
	dudy(i,j,k) = (bp(i,j,k)*u(i,j,k)-bp(i,j-1,k)*u(i,j-1,k))/dy + (bn(i,j+1,k)*u(i,j+1,k)-bn(i,j,k)*u(i,j,k))/dy
	dudz(i,j,k) = (dp(i,j,k)*u(i,j,k)-dp(i,j,k-1)*u(i,j,k-1))/dz(k) + (dn(i,j,k+1)*u(i,j,k+1)-dn(i,j,k)*u(i,j,k))/dz(k)

	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

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
        !$omp target map(to: u, v, w, zero, one_quarter) map(from: ap, an, bp, bn, dp, dn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 0, nz+1
	do j = 1, ny+1
	do i = 0, nx+1

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * one_quarter
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * one_quarter

	ap(i,j,k) = max(amb(i,j,k),zero)
	an(i,j,k) = min(amb(i,j,k),zero)

	bp(i,j,k) = max(v(i,j,k),zero)
	bn(i,j,k) = min(v(i,j,k),zero)

	dp(i,j,k) = max(dmb(i,j,k),zero)
	dn(i,j,k) = min(dmb(i,j,k),zero)


	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

!====================================================================
!	bp(0:nx+1,0,0:nz+1)     = max(v(0:nx+1,0,0:nz+1),0.)
!	bn(0:nx+1,ny1+1,0:nz+1) = min(v(0:nx+1,ny1+1,0:nz+1),0.)
!====================================================================

#ifdef USE_GPU
        !$omp target map(to: v, zero) map(from: bp, bn)
#endif
        !$omp parallel private(i, k)
        !$omp do collapse(2)
        !private(i, k)
        do k = 0, nz+1
        do i = 0, nx+1
            bp(i, 0, k) = max(v(i, 0, k),zero)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2)
        !private(i, k)
        do k = 0, nz+1
        do i = 0, nx+1
            bn(i, ny1+1, k) = min(v(i, ny1+1, k),zero)
        enddo
        enddo
        !$omp end do
        !$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif


#ifdef USE_GPU
        !$omp target map(to: ap, an, bp, bp, bn, dp, dn, v, dx, dy, dz) map(from: vint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k,dvdx,dvdy,dvdz)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	dvdx(i,j,k) = (ap(i,j,k)*v(i,j,k)-ap(i-1,j,k)*v(i-1,j,k))/dx + (an(i+1,j,k)*v(i+1,j,k)-an(i,j,k)*v(i,j,k))/dx
	dvdy(i,j,k) = (bp(i,j,k)*v(i,j,k)-bp(i,j-1,k)*v(i,j-1,k))/dy + (bn(i,j+1,k)*v(i,j+1,k)-bn(i,j,k)*v(i,j,k))/dy
	dvdz(i,j,k) = (dp(i,j,k)*v(i,j,k)-dp(i,j,k-1)*v(i,j,k-1))/dz(k) + (dn(i,j,k+1)*v(i,j,k+1)-dn(i,j,k)*v(i,j,k))/dz(k)

	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

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
        !$omp target map(to: u, v, w, one_quarter) map(from: ap, an, bp, bn, dp, dn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(3) private(i,j,k)
#endif
	do k = 1, nz+1
	do j = 0, ny+1
	do i = 0, nx+1

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * one_quarter
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * one_quarter

	ap(i,j,k) = max(amd(i,j,k),zero)
	an(i,j,k) = min(amd(i,j,k),zero)

	bp(i,j,k) = max(bmd(i,j,k),zero)
	bn(i,j,k) = min(bmd(i,j,k),zero)

	dp(i,j,k) = max(w(i,j,k),zero)
	dn(i,j,k) = min(w(i,j,k),zero)

	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

!====================================================================
!	dp(0:nx+1,0:ny+1,0)     = max(w(0:nx+1,0:ny+1,0),0.)
!	dp(0:nx+1,0:ny+1,nz1+1) = min(w(0:nx+1,0:ny+1,nz1+1),0.)
!====================================================================

#ifdef USE_GPU
        !$omp target map(to: v, zero) map(from: bp, bn)
#endif
        !$omp parallel private(i, j)
        !$omp do collapse(2)
        !private(i, j)
        do j = 0, nz+1
        do i = 0, nx+1
                bp(i, j, 0) = max(v(i, j, 0),zero)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2)
        !private(i, j)
        do j = 0, nz+1
        do i = 0, nx+1
                bn(i, j, nz1+1) = min(v(i, j, nz1+1),zero)
        enddo
        enddo
        !$omp end do
        !$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif

#ifdef USE_GPU
        !$omp target map(to: ap, an, bp, bp, bn, dp, dn, w, dx, dy, dzz) map(from: wint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k,dwdx,dwdy,dwdz)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	dwdx(i,j,k) = (ap(i,j,k)*w(i,j,k)-ap(i-1,j,k)*w(i-1,j,k))/dx + (an(i+1,j,k)*w(i+1,j,k)-an(i,j,k)*w(i,j,k))/dx
	dwdy(i,j,k) = (bp(i,j,k)*w(i,j,k)-bp(i,j-1,k)*w(i,j-1,k))/dy + (bn(i,j+1,k)*w(i,j+1,k)-bn(i,j,k)*w(i,j,k))/dy
	dwdz(i,j,k) = (dp(i,j,k)*w(i,j,k)-dp(i,j,k-1)*w(i,j,k-1))/dzz(k) + (dn(i,j,k+1)*w(i,j,k+1)-dn(i,j,k)*w(i,j,k))/dzz(k)

	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)

	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	else

	write(*,*) "não possui este esquema para o rotacional"
	STOP

	endif


!==================================================================================================================
END SUBROUTINE rotacional


SUBROUTINE antissim(uint,vint,wint)

	USE velpre
	USE param
	!$ USE omp_lib

	IMPLICIT NONE

	!===================================================================================================================
	real(8),intent(out),dimension(nx1,ny,nz) :: uint
	real(8),intent(out),dimension(nx,ny1,nz) :: vint
	real(8),intent(out),dimension(nx,ny,nz1) :: wint

	real(8),dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: ap, an
	real(8),dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: bp, bn
	real(8),dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: dp, dn

	!!
	real(8),dimension(nx1,ny1,nz1) :: dudx, dvdx, dwdx
	real(8),dimension(nx1,ny1,nz1) :: dudy, dvdy, dwdy 
	real(8),dimension(nx1,ny1,nz1) :: dudz, dvdz, dwdz
	real(8),dimension(nx1,ny1,nz1) :: bma, dma,amb, dmb, amd, bmd

	!!
	real(8),save :: aa, bb, dd

	!contadores
	integer :: i, j, k, ai, bi, di

	!plotagem
	real(8),save :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8),save :: aux1, aux2, zero, one_quarter, half

        zero = 0.00000000
        one_quarter = 0.25000000
        half = 0.50000000

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then
	! upwind

#ifdef USE_GPU
        !$omp target map(to: u, v, w, zero, one_quarter) map(from: ap, an, bp, bn, dp, dn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, nz+1
	do j = 1, ny+1
	do i = 1, nx+1
	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * one_quarter
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * one_quarter

	ap(i,j,k) = max(u(i,j,k),zero)
	an(i,j,k) = min(u(i,j,k),zero)

	bp(i,j,k) = max(bma(i,j,k),zero)
	bn(i,j,k) = min(bma(i,j,k),zero)

	dp(i,j,k) = max(dma(i,j,k),zero)
	dn(i,j,k) = min(dma(i,j,k),zero)
	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

!====================================================================
!	ap(0,:,:) = ap(1,:,:)
!	bp(:,0,:) = bp(:,1,:)
!	dp(:,:,0) = dp(:,:,1)
!       an(nx1+1,:,:)=an(nx1,:,:)
!====================================================================

#ifdef USE_GPU
        !$omp target map(tofrom: ap, bp, dp, an)
#endif
        !$omp parallel

        !$omp do collapse(2) private(j, k)
        do k = 1, nz1+1
        do j = 1, ny1+1
            ap(0, j, k) = ap(1, j, k)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, k)
        do k = 1, nz1+1
        do i = 1, nx1+1
            bp(i, 0, k) = bp(i, 1, k)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, j)
        do j = 1, ny1+1
        do i = 1, nx1+1
            dp(i, j, 0) = dp(i, j, 1)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(j, k)
        do k = 1, nz1+1
        do j = 1, ny1+1
            an(nx1+1, j, k) = an(nx1, j, k)
        enddo
        enddo
        !$omp end do
        !$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif

#ifdef USE_GPU
        !$omp target map(to: ap, an, bp, bn, dp, dn, dx, dy, dz, u, half) map(from: uint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k,dudx,dudy,dudz)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	dudx(i,j,k) = ((ap(i,j,k)*u(i,j,k)-ap(i-1,j,k)*u(i-1,j,k))/dx + (an(i+1,j,k)*u(i+1,j,k)-an(i,j,k)*u(i,j,k))/dx & 
	+ ap(i,j,k)*(u(i,j,k)-u(i-1,j,k))/dx + an(i,j,k)*(u(i+1,j,k)-u(i,j,k))/dx) * half


	dudy(i,j,k) = ((bp(i,j,k)*u(i,j,k)-bp(i,j-1,k)*u(i,j-1,k))/dy + (bn(i,j+1,k)*u(i,j+1,k)-bn(i,j,k)*u(i,j,k))/dy &
	+ bp(i,j,k)*(u(i,j,k)-u(i,j-1,k))/dy + bn(i,j,k)*(u(i,j+1,k)-u(i,j,k))/dy) * half

	dudz(i,j,k) = ((dp(i,j,k)*u(i,j,k)-dp(i,j,k-1)*u(i,j,k-1))/dz(k)+(dn(i,j,k+1)*u(i,j,k+1)-dn(i,j,k)*u(i,j,k))/dz(k) &
	+ dp(i,j,k)*(u(i,j,k)-u(i,j,k-1))/dz(k) + dn(i,j,k)*(u(i,j,k+1)-u(i,j,k))/dz(k)) * half


	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

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
        !$omp target map(to: u, v, w, zero, one_quarter) map(tofrom: amb, dmb, ap, an, bp, bn, dp, dn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, nz+1
	do j = 1, ny+1
	do i = 1, nx+1

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * one_quarter
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * one_quarter

	ap(i,j,k) = max(amb(i,j,k),zero)
	an(i,j,k) = min(amb(i,j,k),zero)

	bp(i,j,k) = max(v(i,j,k),zero)
	bn(i,j,k) = min(v(i,j,k),zero)

	dp(i,j,k) = max(dmb(i,j,k),zero)
	dn(i,j,k) = min(dmb(i,j,k),zero)

	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

!====================================================================
!	ap(0,:,:) = ap(1,:,:)
!	bp(:,0,:) = bp(:,1,:)
!	dp(:,:,0) = dp(:,:,1)
!       bn(:,ny1+1,:)=bn(:,ny1,:)
!====================================================================

#ifdef USE_GPU
        !$omp target map(tofrom: ap, bp, dp, bn)
#endif
        !$omp parallel

        !$omp do collapse(2) private(j, k)
        do k = 1, nz1+1
        do j = 1, ny1+1
            ap(0, j, k) = ap(1, j, k)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, k)
        do k = 1, nz1+1
        do i = 1, nx1+1
            bp(i, 0, k) = bp(i, 1, k)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, j)
        do j = 1, ny1+1
        do i = 1, nx1+1
            dp(i, j, 0) = dp(i, j, 1)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, k)
        do k = 1, nz1+1
        do i = 1, ny1+1
            bn(i, ny1+1, k) = bn(i, ny1, k)
        enddo
        enddo
        !$omp end do
        !$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif


#ifdef USE_GPU
        !$omp target map(to: ap, an, bp, bn, dp, dn, dz, v, half) map(tofrom:dvdx, dvdy, dvdz, vint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(2) private(i, j, k)
#endif
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	dvdx(i,j,k) = ((ap(i,j,k)*v(i,j,k)-ap(i-1,j,k)*v(i-1,j,k))/dx + (an(i+1,j,k)*v(i+1,j,k)-an(i,j,k)*v(i,j,k))/dx &
	+ ap(i,j,k)*(v(i,j,k)-v(i-1,j,k))/dx + an(i,j,k)*(v(i+1,j,k)-v(i,j,k))/dx) * half


	dvdy(i,j,k) = ((bp(i,j,k)*v(i,j,k)-bp(i,j-1,k)*v(i,j-1,k))/dy + (bn(i,j+1,k)*v(i,j+1,k)-bn(i,j,k)*v(i,j,k))/dy &
	+ bp(i,j,k)*(v(i,j,k)-v(i,j-1,k))/dy + bn(i,j,k)*(v(i,j+1,k)-v(i,j,k))/dy) * half


	dvdz(i,j,k) = ((dp(i,j,k)*v(i,j,k)-dp(i,j,k-1)*v(i,j,k-1))/dz(k)+(dn(i,j,k+1)*v(i,j,k+1)-dn(i,j,k)*v(i,j,k))/dz(k) &
	+ dp(i,j,k)*(v(i,j,k)-v(i,j,k-1))/dz(k) + dn(i,j,k)*(v(i,j,k+1)-v(i,j,k))/dz(k)) * half


	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

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
        !$omp target map(to: u, v, w, zero, one_quarter) map(tofrom: amd, bmd, ap, an, bp, bn, dp, dn)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, nz+1
	do j = 1, ny+1
	do i = 1, nx+1

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * one_quarter
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * one_quarter

	ap(i,j,k) = max(amd(i,j,k),zero)
	an(i,j,k) = min(amd(i,j,k),zero)

	bp(i,j,k) = max(bmd(i,j,k),zero)
	bn(i,j,k) = min(bmd(i,j,k),zero)

	dp(i,j,k) = max(w(i,j,k),zero)
	dn(i,j,k) = min(w(i,j,k),zero)

	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

!====================================================================
!	ap(0,:,:) = ap(1,:,:)
!	bp(:,0,:) = bp(:,1,:)
!	dp(:,:,0) = dp(:,:,1)
!       dn(:,:,nz1+1)=dn(:,:,nz1)
!====================================================================

#ifdef USE_GPU
        !$omp target map(tofrom: ap, bp, dp, dn)
#endif
        !$omp parallel

        !$omp do collapse(2) private(j, k)
        do k = 1, nz1+1
        do j = 1, ny1+1
            ap(0,j,k) = ap(1,j,k)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, k)
        do k = 1, nz1+1
        do i = 1, nx1+1
            bp(i,0,k) = bp(i,1,k)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, j)
        do j = 1, ny1+1
        do i = 1, nx1+1
            dp(i,j,0) = dp(i,j,1)
        enddo
        enddo
        !$omp end do nowait

        !$omp do collapse(2) private(i, j)
        do j = 1, nz1+1
        do i = 1, ny1+1
            dn(i,j,nz1+1) = dn(i,j,nz1)
        enddo
        enddo
        !$omp end do
        !$omp end parallel
#ifdef USE_GPU
        !$omp end target
#endif


#ifdef USE_GPU
        !$omp target map(to: ap, an, bp, bn, dp, dn, dzz, w, half) map(tofrom:dwdx, dwdy, dwdz, wint)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(2) private(i, j, k)
#endif
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	dwdx(i,j,k) = ((ap(i,j,k)*w(i,j,k)-ap(i-1,j,k)*w(i-1,j,k))/dx + (an(i+1,j,k)*w(i+1,j,k)-an(i,j,k)*w(i,j,k))/dx &
	+ ap(i,j,k)*(w(i,j,k)-w(i-1,j,k))/dx + an(i,j,k)*(w(i+1,j,k)-w(i,j,k))/dx) * half

	dwdy(i,j,k) = ((bp(i,j,k)*w(i,j,k)-bp(i,j-1,k)*w(i,j-1,k))/dy + (bn(i,j+1,k)*w(i,j+1,k)-bn(i,j,k)*w(i,j,k))/dy &
	+ bp(i,j,k)*(w(i,j,k)-w(i,j-1,k))/dy + bn(i,j,k)*(w(i,j+1,k)-w(i,j,k))/dy) * half

	dwdz(i,j,k) = ((dp(i,j,k)*w(i,j,k)-dp(i,j,k-1)*w(i,j,k-1))/dzz(k)+(dn(i,j,k+1)*w(i,j,k+1)-dn(i,j,k)*w(i,j,k))/dzz(k) &
	+ dp(i,j,k)*(w(i,j,k)-w(i,j,k-1))/dzz(k) + dn(i,j,k)*(w(i,j,k+1)-w(i,j,k))/dzz(k)) * half


	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)


	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

	else
	
	write(*,*) "não possui esquema para o antissimétrico"
	STOP

	endif


!==================================================================================================================
END SUBROUTINE antissim

!##################################################################################################################

SUBROUTINE adv_esc(e,dimx,dimy,dimz,adve) ! Advectivo genérico para variáveis escalares centralizadas

	USE disc, only : nx,ny,nz,nx1,ny1,nz1,dz,dzz
	USE velpre, only : u, v, w
	USE param
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k,dimx,dimy,dimz ! contadores
	real(8),save :: aux1, aux2, zero ! auxiliares
	real(8),intent(out),dimension(dimx,dimy,dimz) :: adve        ! termo advectivo
	real(8),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: e           ! e = escalar
	real(8),dimension(nx,ny,nz)  :: uc, vc, wc             ! c = centro
	real(8),dimension(dimx,dimy,dimz)  :: udedx, dedxn, dedxp    ! n = negativo, p = positivo
	real(8),dimension(dimx,dimy,dimz)  :: vdedy, dedyn, dedyp    ! n = negativo, p = positivo
	real(8),dimension(dimx,dimy,dimz)  :: wdedz, dedzn, dedzp    ! n = negativo, p = positivo

        zero = 0.00000000
	! upwind 2nd order

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,uc) ! Interpolação de velocidades da face pro centro
	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,vc)
	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,dzz,wc)

	CALL derivaxu2n(e,dimx,dimy,dimz,dedxn) ! Derivadas do escalar conforme método
	CALL derivaxu2p(e,dimx,dimy,dimz,dedxp)
	CALL derivayu2n(e,dimx,dimy,dimz,dedyn)
	CALL derivayu2p(e,dimx,dimy,dimz,dedyp)
	CALL derivazu2n(e,dimx,dimy,dimz,dz,dedzn)
	CALL derivazu2p(e,dimx,dimy,dimz,dz,dedzp)

#ifdef USE_GPU
        !$omp target map(to: uc, vc, wc, dedxn, dedxp, dedyn, dedyp, dedzn, dedzp, zero) map(tofrom: udedx, vdedy, wdedz, adve)
        !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
        !$omp parallel do collapse(2) private(i,j,k)
#endif
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx

		aux1 = max(uc(i,j,k),zero)
		aux2 = min(uc(i,j,k),zero)
		udedx(i,j,k) = aux1*dedxn(i,j,k) + aux2*dedxp(i,j,k)

		aux1 = max(vc(i,j,k),zero)
		aux2 = min(vc(i,j,k),zero)
		vdedy(i,j,k) = aux1*dedyn(i,j,k) + aux2*dedyp(i,j,k)

		aux1 = max(wc(i,j,k),zero)
		aux2 = min(wc(i,j,k),zero)
		wdedz(i,j,k) = aux1*dedzn(i,j,k) + aux2*dedzp(i,j,k)

		adve(i,j,k) = udedx(i,j,k) + vdedy(i,j,k) + wdedz(i,j,k)

	enddo
	enddo
	enddo
#ifdef USE_GPU 
        !$omp end teams distribute parallel do
        !$omp end target
#else
        !$omp end parallel do
#endif

END SUBROUTINE adv_esc


