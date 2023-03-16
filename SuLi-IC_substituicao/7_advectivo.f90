

SUBROUTINE classico(uint,vint,wint)

	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	real(8), dimension(nx1,ny,nz) :: uint
	real(8), dimension(nx,ny1,nz) :: vint
	real(8), dimension(nx,ny,nz1) :: wint

	real(8), dimension(nx1,ny,nz) :: dudx, dudy, dudz, bma, dma, dudxa, dudya, dudza
	real(8), dimension(nx,ny1,nz) :: dvdx, dvdy, dvdz, amb, dmb, dvdxa, dvdya, dvdza
	real(8), dimension(nx,ny,nz1) :: dwdx, dwdy, dwdz, amd, bmd, dwdxa, dwdya, dwdza
	real(8), dimension(nx,ny,nz)  :: aux

	!contadores
	integer :: i, j, k

	!auxiliares
	real(8) :: aux1, aux2




	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then
	! upwind

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * 0.25
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * 0.25

	aux1 = max(u(i,j,k),0.)
	aux2 = min(u(i,j,k),0.)
	dudx(i,j,k) = aux1*(u(i,j,k)-u(i-1,j,k))/dx + aux2*(u(i+1,j,k)-u(i,j,k))/dx

	aux1 = max(bma(i,j,k),0.)
	aux2 = min(bma(i,j,k),0.)
	dudy(i,j,k) = aux1*(u(i,j,k)-u(i,j-1,k))/dy + aux2*(u(i,j+1,k)-u(i,j,k))/dy

	aux1 = max(dma(i,j,k),0.)
	aux2 = min(dma(i,j,k),0.)
	dudz(i,j,k) = aux1*(u(i,j,k)-u(i,j,k-1))/dz + aux2*(u(i,j,k+1)-u(i,j,k))/dz

	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

	enddo
	enddo
	enddo



	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * 0.25
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * 0.25

	aux1 = max(amb(i,j,k),0.)
	aux2 = min(amb(i,j,k),0.)
	dvdx(i,j,k) = aux1*(v(i,j,k)-v(i-1,j,k))/dx + aux2*(v(i+1,j,k)-v(i,j,k))/dx

	aux1 = max(v(i,j,k),0.)
	aux2 = min(v(i,j,k),0.)
	dvdy(i,j,k) = aux1*(v(i,j,k)-v(i,j-1,k))/dy + aux2*(v(i,j+1,k)-v(i,j,k))/dy

	aux1 = max(dmb(i,j,k),0.)
	aux2 = min(dmb(i,j,k),0.)
	dvdz(i,j,k) = aux1*(v(i,j,k)-v(i,j,k-1))/dz + aux2*(v(i,j,k+1)-v(i,j,k))/dz

	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

	enddo
	enddo
	enddo


	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * 0.25
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * 0.25

	aux1 = max(amd(i,j,k),0.)
	aux2 = min(amd(i,j,k),0.)
	dwdx(i,j,k) = aux1*(w(i,j,k)-w(i-1,j,k))/dx + aux2*(w(i+1,j,k)-w(i,j,k))/dx

	aux1 = max(bmd(i,j,k),0.)
	aux2 = min(bmd(i,j,k),0.)
	dwdy(i,j,k) = aux1*(w(i,j,k)-w(i,j-1,k))/dy + aux2*(w(i,j+1,k)-w(i,j,k))/dy

	aux1 = max(w(i,j,k),0.)
	aux2 = min(w(i,j,k),0.)
	dwdz(i,j,k) = aux1*(w(i,j,k)-w(i,j,k-1))/dz + aux2*(w(i,j,k+1)-w(i,j,k))/dz

	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)

	enddo
	enddo
	enddo


	elseif (der == 2) then


	CALL derivax(u,nx1,ny,nz,dudx)
	CALL derivay(u,nx1,ny,nz,dudy)
	CALL derivaz(u,nx1,ny,nz,dudz)

	CALL derivax(v,nx,ny1,nz,dvdx)
	CALL derivay(v,nx,ny1,nz,dvdy)
	CALL derivaz(v,nx,ny1,nz,dvdz)

	CALL derivax(w,nx,ny,nz1,dwdx)
	CALL derivay(w,nx,ny,nz1,dwdy)
	CALL derivaz(w,nx,ny,nz1,dwdz)

	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,bma) !(nx1,ny,nz)

	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,dma) !(nx1,ny,nz)

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,amb) !(nx,ny1,nz)

	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,dmb) !(nx,ny1,nz)

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,amd) !(nx,ny,nz1)

	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,bmd) !(nx,ny,nz1)

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
	uint(i,j,k) = u(i,j,k)*dudx(i,j,k) + bma(i,j,k)*dudy(i,j,k) + dma(i,j,k)*dudz(i,j,k)
	enddo
	enddo
	enddo

	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
	vint(i,j,k) = amb(i,j,k)*dvdx(i,j,k) + v(i,j,k)*dvdy(i,j,k) + dmb(i,j,k)*dvdz(i,j,k)
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
	wint(i,j,k) = amd(i,j,k)*dwdx(i,j,k) + bmd(i,j,k)*dwdy(i,j,k) + w(i,j,k)*dwdz(i,j,k)
	enddo
	enddo
	enddo



	elseif (der == 3) then
	! upwind 2nd order


	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,bma) !(nx1,ny,nz)

	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,aux) !(nx,ny,nz)
	call interpx_cf(aux,nx,ny,nz,dma) !(nx1,ny,nz)

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,amb) !(nx,ny1,nz)

	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,aux) !(nx,ny,nz)
	call interpy_cf(aux,nx,ny,nz,dmb) !(nx,ny1,nz)

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,amd) !(nx,ny,nz1)

	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,aux) !(nx,ny,nz)
	call interpz_cf(aux,nx,ny,nz,bmd) !(nx,ny,nz1)


	CALL derivaxu2n(u,nx1,ny,nz,dudxa)
	CALL derivaxu2p(u,nx1,ny,nz,dudx)
	CALL derivayu2n(u,nx1,ny,nz,dudya)
	CALL derivayu2p(u,nx1,ny,nz,dudy)
	CALL derivazu2n(u,nx1,ny,nz,dudza)
	CALL derivazu2p(u,nx1,ny,nz,dudz)

	CALL derivaxu2n(v,nx,ny1,nz,dvdxa)
	CALL derivaxu2p(v,nx,ny1,nz,dvdx)
	CALL derivayu2n(v,nx,ny1,nz,dvdya)
	CALL derivayu2p(v,nx,ny1,nz,dvdy)
	CALL derivazu2n(v,nx,ny1,nz,dvdza)
	CALL derivazu2p(v,nx,ny1,nz,dvdz)

	CALL derivaxu2n(w,nx,ny,nz1,dwdxa)
	CALL derivaxu2p(w,nx,ny,nz1,dwdx)
	CALL derivayu2n(w,nx,ny,nz1,dwdya)
	CALL derivayu2p(w,nx,ny,nz1,dwdy)
	CALL derivazu2n(w,nx,ny,nz1,dwdza)
	CALL derivazu2p(w,nx,ny,nz1,dwdz)


	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	aux1 = max(u(i,j,k),0.)
	aux2 = min(u(i,j,k),0.)
	dudx(i,j,k) = aux1*dudxa(i,j,k) + aux2*dudx(i,j,k)

	aux1 = max(bma(i,j,k),0.)
	aux2 = min(bma(i,j,k),0.)
	dudy(i,j,k) = aux1*dudya(i,j,k) + aux2*dudy(i,j,k)

	aux1 = max(dma(i,j,k),0.)
	aux2 = min(dma(i,j,k),0.)
	dudz(i,j,k) = aux1*dudza(i,j,k) + aux2*dudz(i,j,k)

	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

	enddo
	enddo
	enddo

	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	aux1 = max(amb(i,j,k),0.)
	aux2 = min(amb(i,j,k),0.)
	dvdx(i,j,k) = aux1*dvdxa(i,j,k) + aux2*dvdx(i,j,k)

	aux1 = max(v(i,j,k),0.)
	aux2 = min(v(i,j,k),0.)
	dvdy(i,j,k) = aux1*dvdya(i,j,k) + aux2*dvdy(i,j,k)

	aux1 = max(dmb(i,j,k),0.)
	aux2 = min(dmb(i,j,k),0.)
	dvdz(i,j,k) = aux1*dvdza(i,j,k) + aux2*dvdz(i,j,k)

	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	aux1 = max(amd(i,j,k),0.)
	aux2 = min(amd(i,j,k),0.)
	dwdx(i,j,k) = aux1*dwdxa(i,j,k) + aux2*dwdx(i,j,k)

	aux1 = max(bmd(i,j,k),0.)
	aux2 = min(bmd(i,j,k),0.)
	dwdy(i,j,k) = aux1*dwdya(i,j,k) + aux2*dwdy(i,j,k)

	aux1 = max(w(i,j,k),0.)
	aux2 = min(w(i,j,k),0.)
	dwdz(i,j,k) = aux1*dwdza(i,j,k) + aux2*dwdz(i,j,k)

	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)

	enddo
	enddo
	enddo



	endif

!==================================================================================================================
END SUBROUTINE classico


SUBROUTINE rotacional(uint,vint,wint)

	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	real(8), dimension(nx1,ny,nz) :: uint
	real(8), dimension(nx,ny1,nz) :: vint
	real(8), dimension(nx,ny,nz1) :: wint

	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: ap, an, bma, dma
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: bp, bn, amb, dmb
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: dp, dn, amd, bmd

	!
	real(8), dimension(nx1,ny1,nz1) :: dudx, dvdx, dwdx
	real(8), dimension(nx1,ny1,nz1) :: dudy, dvdy, dwdy
	real(8), dimension(nx1,ny1,nz1) :: dudz, dvdz, dwdz

	!
	real(8) :: aa, bb, dd

	!contadores
	integer :: i, j, k, ai, bi, di

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2


	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then
	! upwind

	do k = 0, nz+1
	do j = 0, ny+1
	do i = 1, nx+1

	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * 0.25
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * 0.25

	ap(i,j,k) = max(u(i,j,k),0.)
	an(i,j,k) = min(u(i,j,k),0.)

	bp(i,j,k) = max(bma(i,j,k),0.)
	bn(i,j,k) = min(bma(i,j,k),0.)

	dp(i,j,k) = max(dma(i,j,k),0.)
	dn(i,j,k) = min(dma(i,j,k),0.)

	enddo
	enddo
	enddo

	ap(0,0:ny+1,0:nz+1)     = max(u(0,0:ny+1,0:nz+1),0.)
	an(nx1+1,0:ny+1,0:nz+1) = min(u(nx1+1,0:ny+1,0:nz+1),0.)

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	dudx(i,j,k) = (ap(i,j,k)*u(i,j,k)-ap(i-1,j,k)*u(i-1,j,k))/dx + (an(i+1,j,k)*u(i+1,j,k)-an(i,j,k)*u(i,j,k))/dx
	dudy(i,j,k) = (bp(i,j,k)*u(i,j,k)-bp(i,j-1,k)*u(i,j-1,k))/dy + (bn(i,j+1,k)*u(i,j+1,k)-bn(i,j,k)*u(i,j,k))/dy
	dudz(i,j,k) = (dp(i,j,k)*u(i,j,k)-dp(i,j,k-1)*u(i,j,k-1))/dz + (dn(i,j,k+1)*u(i,j,k+1)-dn(i,j,k)*u(i,j,k))/dz

	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

	enddo
	enddo
	enddo


	do k = 0, nz+1
	do j = 1, ny+1
	do i = 0, nx+1

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * 0.25
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * 0.25

	ap(i,j,k) = max(amb(i,j,k),0.)
	an(i,j,k) = min(amb(i,j,k),0.)

	bp(i,j,k) = max(v(i,j,k),0.)
	bn(i,j,k) = min(v(i,j,k),0.)

	dp(i,j,k) = max(dmb(i,j,k),0.)
	dn(i,j,k) = min(dmb(i,j,k),0.)


	enddo
	enddo
	enddo

	bp(0:nx+1,0,0:nz+1)     = max(v(0:nx+1,0,0:nz+1),0.)
	bn(0:nx+1,ny1+1,0:nz+1) = min(v(0:nx+1,ny1+1,0:nz+1),0.)


	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	dvdx(i,j,k) = (ap(i,j,k)*v(i,j,k)-ap(i-1,j,k)*v(i-1,j,k))/dx + (an(i+1,j,k)*v(i+1,j,k)-an(i,j,k)*v(i,j,k))/dx
	dvdy(i,j,k) = (bp(i,j,k)*v(i,j,k)-bp(i,j-1,k)*v(i,j-1,k))/dy + (bn(i,j+1,k)*v(i,j+1,k)-bn(i,j,k)*v(i,j,k))/dy
	dvdz(i,j,k) = (dp(i,j,k)*v(i,j,k)-dp(i,j,k-1)*v(i,j,k-1))/dz + (dn(i,j,k+1)*v(i,j,k+1)-dn(i,j,k)*v(i,j,k))/dz

	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

	enddo
	enddo
	enddo

	do k = 1, nz+1
	do j = 0, ny+1
	do i = 0, nx+1

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * 0.25
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * 0.25

	ap(i,j,k) = max(amd(i,j,k),0.)
	an(i,j,k) = min(amd(i,j,k),0.)

	bp(i,j,k) = max(bmd(i,j,k),0.)
	bn(i,j,k) = min(bmd(i,j,k),0.)

	dp(i,j,k) = max(w(i,j,k),0.)
	dn(i,j,k) = min(w(i,j,k),0.)

	enddo
	enddo
	enddo

	dp(0:nx+1,0:ny+1,0)     = max(w(0:nx+1,0:ny+1,0),0.)
	dp(0:nx+1,0:ny+1,nz1+1) = min(w(0:nx+1,0:ny+1,nz1+1),0.)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	dwdx(i,j,k) = (ap(i,j,k)*w(i,j,k)-ap(i-1,j,k)*w(i-1,j,k))/dx + (an(i+1,j,k)*w(i+1,j,k)-an(i,j,k)*w(i,j,k))/dx
	dwdy(i,j,k) = (bp(i,j,k)*w(i,j,k)-bp(i,j-1,k)*w(i,j-1,k))/dy + (bn(i,j+1,k)*w(i,j+1,k)-bn(i,j,k)*w(i,j,k))/dy
	dwdz(i,j,k) = (dp(i,j,k)*w(i,j,k)-dp(i,j,k-1)*w(i,j,k-1))/dz + (dn(i,j,k+1)*w(i,j,k+1)-dn(i,j,k)*w(i,j,k))/dz

	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)

	enddo
	enddo
	enddo

	else

	write(*,*) "não possui este esquema para o rotacional"
	STOP

	endif


!==================================================================================================================
END SUBROUTINE rotacional


SUBROUTINE antissim(uint,vint,wint)

	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	real(8), dimension(nx1,ny,nz) :: uint
	real(8), dimension(nx,ny1,nz) :: vint
	real(8), dimension(nx,ny,nz1) :: wint

	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: ap, an
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: bp, bn
	real(8), dimension(0:nx1+1,0:ny1+1,0:nz1+1) :: dp, dn

	!
	real(8), dimension(nx1,ny1,nz1) :: dudx, dvdx, dwdx
	real(8), dimension(nx1,ny1,nz1) :: dudy, dvdy, dwdy 
	real(8), dimension(nx1,ny1,nz1) :: dudz, dvdz, dwdz
	real(8), dimension(nx1,ny1,nz1) :: bma, dma,amb, dmb, amd, bmd

	!
	real(8) :: aa, bb, dd

	!contadores
	integer :: i, j, k, ai, bi, di

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2


	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then
	! upwind

	do k = 1, nz+1
	do j = 1, ny+1
	do i = 1, nx+1
	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * 0.25
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * 0.25

	ap(i,j,k) = max(u(i,j,k),0.)
	an(i,j,k) = min(u(i,j,k),0.)

	bp(i,j,k) = max(bma(i,j,k),0.)
	bn(i,j,k) = min(bma(i,j,k),0.)

	dp(i,j,k) = max(dma(i,j,k),0.)
	dn(i,j,k) = min(dma(i,j,k),0.)
	enddo
	enddo
	enddo

	ap(0,:,:) = ap(1,:,:)
	bp(:,0,:) = bp(:,1,:)
	dp(:,:,0) = dp(:,:,1)
        an(nx1+1,:,:)=an(nx1,:,:)


	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	dudx(i,j,k) = ((ap(i,j,k)*u(i,j,k)-ap(i-1,j,k)*u(i-1,j,k))/dx + (an(i+1,j,k)*u(i+1,j,k)-an(i,j,k)*u(i,j,k))/dx & 
	+ ap(i,j,k)*(u(i,j,k)-u(i-1,j,k))/dx + an(i,j,k)*(u(i+1,j,k)-u(i,j,k))/dx) * 0.5


	dudy(i,j,k) = ((bp(i,j,k)*u(i,j,k)-bp(i,j-1,k)*u(i,j-1,k))/dy + (bn(i,j+1,k)*u(i,j+1,k)-bn(i,j,k)*u(i,j,k))/dy &
	+ bp(i,j,k)*(u(i,j,k)-u(i,j-1,k))/dy + bn(i,j,k)*(u(i,j+1,k)-u(i,j,k))/dy) * 0.5

	dudz(i,j,k) = ((dp(i,j,k)*u(i,j,k)-dp(i,j,k-1)*u(i,j,k-1))/dz+(dn(i,j,k+1)*u(i,j,k+1)-dn(i,j,k)*u(i,j,k))/dz &
	+ dp(i,j,k)*(u(i,j,k)-u(i,j,k-1))/dz + dn(i,j,k)*(u(i,j,k+1)-u(i,j,k))/dz) * 0.5


	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

	enddo
	enddo
	enddo


	do k = 1, nz+1
	do j = 1, ny+1
	do i = 1, nx+1

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * 0.25
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * 0.25

	ap(i,j,k) = max(amb(i,j,k),0.)
	an(i,j,k) = min(amb(i,j,k),0.)

	bp(i,j,k) = max(v(i,j,k),0.)
	bn(i,j,k) = min(v(i,j,k),0.)

	dp(i,j,k) = max(dmb(i,j,k),0.)
	dn(i,j,k) = min(dmb(i,j,k),0.)


	enddo
	enddo
	enddo


	ap(0,:,:) = ap(1,:,:)
	bp(:,0,:) = bp(:,1,:)
	dp(:,:,0) = dp(:,:,1)
        bn(:,ny1+1,:)=bn(:,ny1,:)

	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	dvdx(i,j,k) = ((ap(i,j,k)*v(i,j,k)-ap(i-1,j,k)*v(i-1,j,k))/dx + (an(i+1,j,k)*v(i+1,j,k)-an(i,j,k)*v(i,j,k))/dx &
	+ ap(i,j,k)*(v(i,j,k)-v(i-1,j,k))/dx + an(i,j,k)*(v(i+1,j,k)-v(i,j,k))/dx) * 0.5


	dvdy(i,j,k) = ((bp(i,j,k)*v(i,j,k)-bp(i,j-1,k)*v(i,j-1,k))/dy + (bn(i,j+1,k)*v(i,j+1,k)-bn(i,j,k)*v(i,j,k))/dy &
	+ bp(i,j,k)*(v(i,j,k)-v(i,j-1,k))/dy + bn(i,j,k)*(v(i,j+1,k)-v(i,j,k))/dy) * 0.5


	dvdz(i,j,k) = ((dp(i,j,k)*v(i,j,k)-dp(i,j,k-1)*v(i,j,k-1))/dz+(dn(i,j,k+1)*v(i,j,k+1)-dn(i,j,k)*v(i,j,k))/dz &
	+ dp(i,j,k)*(v(i,j,k)-v(i,j,k-1))/dz + dn(i,j,k)*(v(i,j,k+1)-v(i,j,k))/dz) * 0.5


	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

	enddo
	enddo
	enddo

	do k = 1, nz+1
	do j = 1, ny+1
	do i = 1, nx+1

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * 0.25
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * 0.25

	ap(i,j,k) = max(amd(i,j,k),0.)
	an(i,j,k) = min(amd(i,j,k),0.)

	bp(i,j,k) = max(bmd(i,j,k),0.)
	bn(i,j,k) = min(bmd(i,j,k),0.)

	dp(i,j,k) = max(w(i,j,k),0.)
	dn(i,j,k) = min(w(i,j,k),0.)

	enddo
	enddo
	enddo

	ap(0,:,:) = ap(1,:,:)
	bp(:,0,:) = bp(:,1,:)
	dp(:,:,0) = dp(:,:,1)
        dn(:,:,nz1+1)=dn(:,:,nz1)

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	dwdx(i,j,k) = ((ap(i,j,k)*w(i,j,k)-ap(i-1,j,k)*w(i-1,j,k))/dx + (an(i+1,j,k)*w(i+1,j,k)-an(i,j,k)*w(i,j,k))/dx &
	+ ap(i,j,k)*(w(i,j,k)-w(i-1,j,k))/dx + an(i,j,k)*(w(i+1,j,k)-w(i,j,k))/dx) * 0.5

	dwdy(i,j,k) = ((bp(i,j,k)*w(i,j,k)-bp(i,j-1,k)*w(i,j-1,k))/dy + (bn(i,j+1,k)*w(i,j+1,k)-bn(i,j,k)*w(i,j,k))/dy &
	+ bp(i,j,k)*(w(i,j,k)-w(i,j-1,k))/dy + bn(i,j,k)*(w(i,j+1,k)-w(i,j,k))/dy) * 0.5

	dwdz(i,j,k) = ((dp(i,j,k)*w(i,j,k)-dp(i,j,k-1)*w(i,j,k-1))/dz+(dn(i,j,k+1)*w(i,j,k+1)-dn(i,j,k)*w(i,j,k))/dz &
	+ dp(i,j,k)*(w(i,j,k)-w(i,j,k-1))/dz + dn(i,j,k)*(w(i,j,k+1)-w(i,j,k))/dz) * 0.5


	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)


	enddo
	enddo
	enddo

	else

	write(*,*) "não possui esquema para o antissimétrico"
	STOP

	endif


!==================================================================================================================
END SUBROUTINE antissim
