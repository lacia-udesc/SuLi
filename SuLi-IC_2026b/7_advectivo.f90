

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

	real(8),dimension(nx1,ny,nz) :: dudx, dudy, dudz, bma, dma, dudxa, dudya, dudza, dudxp, dudyp, dudzp
	real(8),dimension(nx,ny1,nz) :: dvdx, dvdy, dvdz, amb, dmb, dvdxa, dvdya, dvdza, dvdxp, dvdyp, dvdzp
	real(8),dimension(nx,ny,nz1) :: dwdx, dwdy, dwdz, amd, bmd, dwdxa, dwdya, dwdza, dwdxp, dwdyp, dwdzp
	real(8),dimension(nx,ny,nz)  :: aux

	!auxiliares
	real(8) :: aux1, aux2, aux3, aux4, aux5, aux6, zero, one_quarter
	real(8) :: cdu, upx, upy, upz ! termos centrado e upwind do esquema hibrido

	zero = 0.0d0
	one_quarter = 0.25d0

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then
	! upwind


    !$omp parallel private(aux1,aux2,aux3,aux4,aux5,aux6)
    !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * one_quarter
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * one_quarter

	aux1 = max(u(i,j,k),zero)
	aux2 = min(u(i,j,k),zero)
	dudx(i,j,k) = aux1*(u(i,j,k)-u(i-1,j,k))*inv_dx + aux2*(u(i+1,j,k)-u(i,j,k))*inv_dx

	aux1 = max(bma(i,j,k),zero)
	aux2 = min(bma(i,j,k),zero)
	dudy(i,j,k) = aux1*(u(i,j,k)-u(i,j-1,k))*inv_dy + aux2*(u(i,j+1,k)-u(i,j,k))*inv_dy

	aux1 = max(dma(i,j,k),zero)
	aux2 = min(dma(i,j,k),zero)
	dudz(i,j,k) = aux1*(u(i,j,k)-u(i,j,k-1))*inv_dz(k) + aux2*(u(i,j,k+1)-u(i,j,k))*inv_dz(k)

	uint(i,j,k) = dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k)

	enddo
	enddo
	enddo
    !$omp end do nowait

    !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * one_quarter
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * one_quarter

	aux3 = max(amb(i,j,k),zero)
	aux4 = min(amb(i,j,k),zero)
	dvdx(i,j,k) = aux3*(v(i,j,k)-v(i-1,j,k))*inv_dx + aux4*(v(i+1,j,k)-v(i,j,k))*inv_dx

	aux3 = max(v(i,j,k),zero)
	aux4 = min(v(i,j,k),zero)
	dvdy(i,j,k) = aux3*(v(i,j,k)-v(i,j-1,k))*inv_dy + aux4*(v(i,j+1,k)-v(i,j,k))*inv_dy

	aux3 = max(dmb(i,j,k),zero)
	aux4 = min(dmb(i,j,k),zero)
	dvdz(i,j,k) = aux3*(v(i,j,k)-v(i,j,k-1))*inv_dz(k) + aux4*(v(i,j,k+1)-v(i,j,k))*inv_dz(k)

	vint(i,j,k) = dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k)

	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do collapse(3)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * one_quarter
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * one_quarter

	aux5 = max(amd(i,j,k),zero)
	aux6 = min(amd(i,j,k),zero)
	dwdx(i,j,k) = aux5*(w(i,j,k)-w(i-1,j,k))*inv_dx + aux6*(w(i+1,j,k)-w(i,j,k))*inv_dx

	aux5 = max(bmd(i,j,k),zero)
	aux6 = min(bmd(i,j,k),zero)
	dwdy(i,j,k) = aux5*(w(i,j,k)-w(i,j-1,k))*inv_dy + aux6*(w(i,j+1,k)-w(i,j,k))*inv_dy

	aux5 = max(w(i,j,k),zero)
	aux6 = min(w(i,j,k),zero)
	dwdz(i,j,k) = aux5*(w(i,j,k)-w(i,j,k-1))*inv_dzz(k) + aux6*(w(i,j,k+1)-w(i,j,k))*inv_dzz(k)

	wint(i,j,k) = dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k)

	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel


	elseif (der == 2) then
	! centrado de segunda ordem

	CALL derivax(u,nx1,ny,nz,dudx)
	CALL derivay(u,nx1,ny,nz,dudy)
	CALL derivaz(u,nx1,ny,nz,inv_2dz,dudz)

	CALL derivax(v,nx,ny1,nz,dvdx)
	CALL derivay(v,nx,ny1,nz,dvdy)
	CALL derivaz(v,nx,ny1,nz,inv_2dz,dvdz)

	CALL derivax(w,nx,ny,nz1,dwdx)
	CALL derivay(w,nx,ny,nz1,dwdy)
	CALL derivaz(w,nx,ny,nz1,inv_2dzz,dwdz)

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

	elseif (der == 3) then
	! upwind 2nd order

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


	CALL derivaxu2n(u,nx1,ny,nz,dudxa)
	CALL derivaxu2p(u,nx1,ny,nz,dudx)
	CALL derivayu2n(u,nx1,ny,nz,dudya)
	CALL derivayu2p(u,nx1,ny,nz,dudy)
	CALL derivazu2n(u,nx1,ny,nz,inv_2dz,dudza)
	CALL derivazu2p(u,nx1,ny,nz,inv_2dz,dudz)

	CALL derivaxu2n(v,nx,ny1,nz,dvdxa)
	CALL derivaxu2p(v,nx,ny1,nz,dvdx)
	CALL derivayu2n(v,nx,ny1,nz,dvdya)
	CALL derivayu2p(v,nx,ny1,nz,dvdy)
	CALL derivazu2n(v,nx,ny1,nz,inv_2dz,dvdza)
	CALL derivazu2p(v,nx,ny1,nz,inv_2dz,dvdz)

	CALL derivaxu2n(w,nx,ny,nz1,dwdxa)
	CALL derivaxu2p(w,nx,ny,nz1,dwdx)
	CALL derivayu2n(w,nx,ny,nz1,dwdya)
	CALL derivayu2p(w,nx,ny,nz1,dwdy)
	CALL derivazu2n(w,nx,ny,nz1,inv_2dzz,dwdza)
	CALL derivazu2p(w,nx,ny,nz1,inv_2dzz,dwdz)


    !$omp parallel private(aux1,aux2,aux3,aux4,aux5,aux6)
    !$omp do collapse(3)
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

    !$omp do collapse(3)
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

    !$omp do collapse(3)
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

	elseif (der == 4) then
	! hibrido: centrado de 2a ordem + upwind de 2a ordem (fator de mistura beta_cu)


	! cálculo do centrado
	CALL derivax(u,nx1,ny,nz,dudx)
	CALL derivay(u,nx1,ny,nz,dudy)
	CALL derivaz(u,nx1,ny,nz,inv_2dz,dudz)

	CALL derivax(v,nx,ny1,nz,dvdx)
	CALL derivay(v,nx,ny1,nz,dvdy)
	CALL derivaz(v,nx,ny1,nz,inv_2dz,dvdz)

	CALL derivax(w,nx,ny,nz1,dwdx)
	CALL derivay(w,nx,ny,nz1,dwdy)
	CALL derivaz(w,nx,ny,nz1,inv_2dzz,dwdz)

	!cálculo do upwind
	CALL derivaxu2n(u,nx1,ny,nz,dudxa)
	CALL derivaxu2p(u,nx1,ny,nz,dudxp)
	CALL derivayu2n(u,nx1,ny,nz,dudya)
	CALL derivayu2p(u,nx1,ny,nz,dudyp)
	CALL derivazu2n(u,nx1,ny,nz,inv_2dz,dudza)
	CALL derivazu2p(u,nx1,ny,nz,inv_2dz,dudzp)

	CALL derivaxu2n(v,nx,ny1,nz,dvdxa)
	CALL derivaxu2p(v,nx,ny1,nz,dvdxp)
	CALL derivayu2n(v,nx,ny1,nz,dvdya)
	CALL derivayu2p(v,nx,ny1,nz,dvdyp)
	CALL derivazu2n(v,nx,ny1,nz,inv_2dz,dvdza)
	CALL derivazu2p(v,nx,ny1,nz,inv_2dz,dvdzp)

	CALL derivaxu2n(w,nx,ny,nz1,dwdxa)
	CALL derivaxu2p(w,nx,ny,nz1,dwdxp)
	CALL derivayu2n(w,nx,ny,nz1,dwdya)
	CALL derivayu2p(w,nx,ny,nz1,dwdyp)
	CALL derivazu2n(w,nx,ny,nz1,inv_2dzz,dwdza)
	CALL derivazu2p(w,nx,ny,nz1,inv_2dzz,dwdzp)

	!alteração da posição do vetor velociade para o cálculo da advecção
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

    !$omp parallel private(aux1,aux2,aux3,aux4,aux5,aux6,cdu,upx,upy,upz)
    !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1

	!cálculo centrado
	cdu = u(i,j,k)*dudx(i,j,k) + bma(i,j,k)*dudy(i,j,k) + dma(i,j,k)*dudz(i,j,k)

	!cálculo upwind
	aux1 = max(u(i,j,k),zero)
	aux2 = min(u(i,j,k),zero)
	upx = aux1*dudxa(i,j,k) + aux2*dudxp(i,j,k)

	aux1 = max(bma(i,j,k),zero)
	aux2 = min(bma(i,j,k),zero)
	upy = aux1*dudya(i,j,k) + aux2*dudyp(i,j,k)

	aux1 = max(dma(i,j,k),zero)
	aux2 = min(dma(i,j,k),zero)
	upz = aux1*dudza(i,j,k) + aux2*dudzp(i,j,k)

	!junção
	uint(i,j,k) = (1.0d0 - beta_cu)*cdu + beta_cu*(upx + upy + upz)

	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx

	!cálculo centrado
	cdu = amb(i,j,k)*dvdx(i,j,k) + v(i,j,k)*dvdy(i,j,k) + dmb(i,j,k)*dvdz(i,j,k)

	aux3 = max(amb(i,j,k),zero)
	aux4 = min(amb(i,j,k),zero)
	upx = aux3*dvdxa(i,j,k) + aux4*dvdxp(i,j,k)

	!cálculo upwind
	aux3 = max(v(i,j,k),zero)
	aux4 = min(v(i,j,k),zero)
	upy = aux3*dvdya(i,j,k) + aux4*dvdyp(i,j,k)

	aux3 = max(dmb(i,j,k),zero)
	aux4 = min(dmb(i,j,k),zero)
	upz = aux3*dvdza(i,j,k) + aux4*dvdzp(i,j,k)

	!junção
	vint(i,j,k) = (1.0d0 - beta_cu)*cdu + beta_cu*(upx + upy + upz)

	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do collapse(3)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx

	!cálculo centrado
	cdu = amd(i,j,k)*dwdx(i,j,k) + bmd(i,j,k)*dwdy(i,j,k) + w(i,j,k)*dwdz(i,j,k)

	aux5 = max(amd(i,j,k),zero)
	aux6 = min(amd(i,j,k),zero)
	upx = aux5*dwdxa(i,j,k) + aux6*dwdxp(i,j,k)

	!cálculo upwind
	aux5 = max(bmd(i,j,k),zero)
	aux6 = min(bmd(i,j,k),zero)
	upy = aux5*dwdya(i,j,k) + aux6*dwdyp(i,j,k)

	!junção
	aux5 = max(w(i,j,k),zero)
	aux6 = min(w(i,j,k),zero)
	upz = aux5*dwdza(i,j,k) + aux6*dwdzp(i,j,k)

	wint(i,j,k) = (1.0d0 - beta_cu)*cdu + beta_cu*(upx + upy + upz)

	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel

	endif

!==================================================================================================================
END SUBROUTINE classico


!##################################################################################################################

SUBROUTINE adv_esc(e,dimx,dimy,dimz,adve) ! Advectivo genérico para variáveis escalares centralizadas

	USE disc, only : nx,ny,nz,nx1,ny1,nz1,dz,dzz,inv_2dz
	USE velpre, only : u, v, w
	USE param
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k,dimx,dimy,dimz ! contadores
	real(8) :: aux1, aux2, zero ! auxiliares
	real(8),intent(out),dimension(dimx,dimy,dimz) :: adve        ! termo advectivo
	real(8),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: e           ! e = escalar
	real(8),dimension(nx,ny,nz)  :: uc, vc, wc             ! c = centro
	real(8),dimension(dimx,dimy,dimz)  :: udedx, dedxn, dedxp    ! n = negativo, p = positivo
	real(8),dimension(dimx,dimy,dimz)  :: vdedy, dedyn, dedyp    ! n = negativo, p = positivo
	real(8),dimension(dimx,dimy,dimz)  :: wdedz, dedzn, dedzp    ! n = negativo, p = positivo

    zero = 0.0d0
	! upwind 2nd order

	call interpx_fc(u(1:nx1,1:ny,1:nz),nx1,ny,nz,uc) ! Interpolação de velocidades da face pro centro
	call interpy_fc(v(1:nx,1:ny1,1:nz),nx,ny1,nz,vc)
	call interpz_fc(w(1:nx,1:ny,1:nz1),nx,ny,nz1,dzz,wc)

	CALL derivaxu2n(e,dimx,dimy,dimz,dedxn) ! Derivadas do escalar conforme método
	CALL derivaxu2p(e,dimx,dimy,dimz,dedxp)
	CALL derivayu2n(e,dimx,dimy,dimz,dedyn)
	CALL derivayu2p(e,dimx,dimy,dimz,dedyp)
	CALL derivazu2n(e,dimx,dimy,dimz,inv_2dz,dedzn)
	CALL derivazu2p(e,dimx,dimy,dimz,inv_2dz,dedzp)


    !$omp parallel do collapse(3) private(aux1,aux2)
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
    !$omp end parallel do


END SUBROUTINE adv_esc

