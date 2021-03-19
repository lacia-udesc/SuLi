	!Subrotina para calcular o Fu, utilizando o ponto de vista lagrangiano
	!Referencia: Casulli (1990, 1992)

	!!! Implementação 15/04/2014
	! Leonardo Romero Monteiro

	!!! Modificações
	! Leonardo Romero Monteiro - 13/01/2015

	!!! Observações
	! no lagrangiano as velocidades no tempo n = t serão calculadas para o tempo n = t-1

	! são calculadas as Fu, Fv e Fw
	SUBROUTINE lagr(uint,vint,wint)

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	! auxiliares de velocidades: velocidades lagrangianas
	real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: uint
	real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: vint
	real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: wint

	! números de Courant e auxiliares de velocidades: velocidades eulerianas com contornos adicionais
	real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: a, uu
	real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: b, vv
	real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: d, ww

	!
	real(8), dimension(nx1,ny,nz) :: bma, dma
	real(8), dimension(nx,ny1,nz) :: amb, dmb
	real(8), dimension(nx,ny,nz1) :: amd, bmd

	!
	real(8) :: aa, bb, dd

	!contadores
	integer :: i, j, k, ai, bi, di
	integer :: ntal, tt
	real(8)    :: tal

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2

	!===================================================================================================================
	

	!reinicializando a subrotina
	a = 0.
	b = 0.
	d = 0.
	aa = 0.
	bb = 0.
	dd = 0.
	ai = 0
	bi = 0
	di = 0
	acont = 0.
	bcont = 0.
	dcont = 0.


	!as velocidades iniciais do tempo n = t
	uint = u
	vint = v
	wint = w
	uu = u
	vv = v
	ww = w

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	!%%%!--VELOCIDADES LAGRANGIANAS--!%%%!

	! verifica se o número de courant será mais do que 1. Se for, será fracionado o tempo em "tal" partes (multi-step backward Euler - MSE: ver dissertação)
	if (maxval(abs(v))/dy > maxval(abs(u))/dx .or. (maxval(abs(v))/dy > maxval(abs(w)/dzz))) then
		ntal = ceiling((maxval(abs(v))/dy)*dt+0.1)
	elseif (maxval(abs(u))/dx > maxval(abs(w)/dzz)) then
		ntal = ceiling((maxval(abs(u))/dx)*dt+0.1)
	else
		ntal = ceiling((maxval(abs(w)/dzz))*dt+0.1)
	endif


	tal = dt/ntal

	! cálculo dos números de courant
	do tt = 1, ntal

	aux1 = tal / dx
	aux2 = tal / dy



	a = a + uint*aux1
	b = b + vint*aux2
	d = d + wint*tal/dzz
	
	! cálculo das velocidades lagrangianas (t = n-1 - backward step)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		! para deixar os números de grid courant em cima de suas respectivas velocidades
		bma(i,j,k) = (b(i,j,k) + b(i,j+1,k) + b(i-1,j,k) + b(i-1,j+1,k) ) * 0.25
		dma(i,j,k) = (d(i,j,k) + d(i,j,k+1) + d(i-1,j,k) + d(i-1,j,k+1) ) * 0.25
		! como se está voltando no tempo as velocidades positivas serão vistas em distância "-1" e as negativas serão vistas em "1"
		if (a(i,j,k) >= 0 ) then
			ai = 1
		else
			ai = -1
		endif
		if (bma(i,j,k) >= 0 ) then
			bi = 1
		else
			bi = -1
		endif
		if (dma(i,j,k) >= 0 ) then
			di = 1
		else
			di = -1
		endif

		! como já foi definido o sinal anteriormente aqui devemos abandonar o sinal
		aa = abs(a(i,j,k))
		bb = abs(bma(i,j,k))
		dd = abs(dma(i,j,k))
						
		uint(i,j,k) = (1-dd) * ( (1-aa) * ( (1-bb) * uu(i,j,k) + bb * uu(i,j-bi,k) ) &
		+ aa * ( (1-bb) * uu(i-ai,j,k) + bb * uu(i-ai,j-bi,k))) &
		+ dd * ( (1-aa) * ( (1-bb) * uu(i,j,k-di) + bb * uu(i,j-bi,k-di) ) &
		+ aa * ( (1-bb) * uu(i-ai,j,k-di) + bb * uu(i-ai,j-bi,k-di) ) )
	enddo
	enddo

	do j = 1, ny1
	do i = 1, nx
		amb(i,j,k) =  (a(i,j,k) + a(i+1,j,k) + a(i,j-1,k) + a(i+1,j-1,k) ) * 0.25
		dmb(i,j,k) =  (d(i,j,k) + d(i,j,k+1) + d(i,j-1,k) + d(i,j-1,k+1) ) * 0.25

		if (amb(i,j,k) >= 0 ) then
			ai = 1
		else
			ai = -1
		endif

		if (b(i,j,k) >= 0 ) then
			bi = 1
		else
			bi = -1
		endif

		if (dmb(i,j,k) >= 0 ) then
			di = 1
		else
			di = -1
		endif

		aa = abs(amb(i,j,k))
		bb = abs(b(i,j,k))
		dd = abs(dmb(i,j,k))

		vint(i,j,k) = (1-dd) * ( (1-aa) * ( (1-bb) * vv(i,j,k) + bb * vv(i,j-bi,k) ) &
		+ aa * ( (1-bb) * vv(i-ai,j,k) + bb * vv(i-ai,j-bi,k))) &
		+ dd * ( (1-aa) * ( (1-bb) * vv(i,j,k-di) + bb * vv(i,j-bi,k-di) ) &
		+ aa * ( (1-bb) * vv(i-ai,j,k-di) + bb * vv(i-ai,j-bi,k-di) ) )
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		
		amd(i,j,k) =  (a(i,j,k) + a(i+1,j,k) + a(i,j,k-1) + a(i+1,j,k-1) ) * 0.25
		bmd(i,j,k) =  (b(i,j,k) + b(i,j+1,k) + b(i,j,k-1) + b(i,j+1,k-1) ) * 0.25

		if (amd(i,j,k) >= 0 ) then
			ai = 1
		else
			ai = -1
		endif

		if (bmd(i,j,k) >= 0 ) then
			bi = 1
		else
			bi = -1
		endif

		if (d(i,j,k) >= 0 ) then
			di = 1
		else
			di = -1
		endif

		aa = abs(amd(i,j,k))
		bb = abs(bmd(i,j,k))
		dd = abs(d(i,j,k))

		wint(i,j,k) = (1-dd) * ( (1-aa) * ( (1-bb) * ww(i,j,k) + bb * ww(i,j-bi,k) ) &
		+ aa * ( (1-bb) * ww(i-ai,j,k) + bb * ww(i-ai,j-bi,k))) &
		+ dd * ( (1-aa) * ( (1-bb) * ww(i,j,k-di) + bb * ww(i,j-bi,k-di) ) &
		+ aa * ( (1-bb) * ww(i-ai,j,k-di) + bb * ww(i-ai,j-bi,k-di) ) )
	enddo
	enddo
	enddo

	enddo !fechou o tempo


	!===============================================================================================================
	! Plotagem dos Números de Courant
	acont = maxval(abs(a))
	bcont = maxval(abs(b))
	dcont = maxval(abs(d))

	loca = maxloc(abs(a))
	locb = maxloc(abs(b))
	locd = maxloc(abs(d))

	!Números de Courant **************************************************
	if (it == 1) open (unit=9999991, action= 'write', file= 'dados//courant.txt', status= 'unknown')
	write(9999991,*) it*dt, ntal, tal, acont, bcont, dcont, loca, locb, locd
	!*********************************************************************

	! interrompe simulalação por desrespeitar o CFL
	if (acont > 1. .or. bcont > 1 .or. dcont > 1) then
		write(*,*) "Pelo menos um dos valores de CFL ultrafassou o limite de 1"
		write(*,*) "acont", acont, "bcont", bcont, "dcont", dcont
		stop
	endif

	if (der .ne. 1) then
		write(*,*) "A variação de esquema não altera este método"
	endif
	
END SUBROUTINE lagr

!==================================================================================================================

SUBROUTINE classico(uint,vint,wint)

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	! auxiliares de velocidades: velocidades lagrangianas
	real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: uint
	real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: vint
	real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: wint

	! números de Courant e auxiliares de velocidades: velocidades eulerianas com contornos adicionais
	real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: a, uu
	real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: b, vv
	real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: d, ww

	!
	real(8), dimension(nx1,ny,nz) :: dudx, dudy, dudz, bma, dma
	real(8), dimension(nx,ny1,nz) :: dvdx, dvdy, dvdz, amb, dmb
	real(8), dimension(nx,ny,nz1) :: dwdx, dwdy, dwdz, amd, bmd

	!
	real(8) :: aa, bb, dd

	!contadores
	integer :: i, j, k, ai, bi, di
	integer :: ntal, tt
	real(8)    :: tal

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2


	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then ! upwind 1st

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
	dudz(i,j,k) = aux1*(u(i,j,k)-u(i,j,k-1))/dzx(i,j,k) + aux2*(u(i,j,k+1)-u(i,j,k))/dzx(i,j,k)

	uint(i,j,k) = u(i,j,k) - (dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k))*dt

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
	dvdz(i,j,k) = aux1*(v(i,j,k)-v(i,j,k-1))/dzy(i,j,k) + aux2*(v(i,j,k+1)-v(i,j,k))/dzy(i,j,k)

	vint(i,j,k) = v(i,j,k) - (dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k))*dt

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
	dwdz(i,j,k) = aux1*(w(i,j,k)-w(i,j,k-1))/dzz(i,j,k) + aux2*(w(i,j,k+1)-w(i,j,k))/dzz(i,j,k)

	wint(i,j,k) = w(i,j,k) - (dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k))*dt

	enddo
	enddo
	enddo


	elseif (der == 2) then !centrado 2nd

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx

	!defirenças centras
	dudx(i,j,k) = (u(i+1,j,k) - u(i-1,j,k)) / (2*dx)
	dudy(i,j,k) = (u(i,j+1,k) - u(i,j-1,k)) / (2*dy)
	dudz(i,j,k) = (u(i,j,k+1) - u(i,j,k-1)) / (2*dzx(i,j,k))

	dvdx(i,j,k) = (v(i+1,j,k) - v(i-1,j,k)) / (2*dx)
	dvdy(i,j,k) = (v(i,j+1,k) - v(i,j-1,k)) / (2*dy)
	dvdz(i,j,k) = (v(i,j,k+1) - v(i,j,k-1)) / (2*dzy(i,j,k))

	dwdx(i,j,k) = (w(i+1,j,k) - w(i-1,j,k)) / (2*dx)
	dwdy(i,j,k) = (w(i,j+1,k) - w(i,j-1,k)) / (2*dy)
	dwdz(i,j,k) = (w(i,j,k+1) - w(i,j,k-1)) / (2*dzz(i,j,k))

	enddo
	enddo
	enddo

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	! interpolações
	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * 0.25
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * 0.25

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * 0.25
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * 0.25

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * 0.25
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * 0.25

	uint(i,j,k) = u(i,j,k) - (u(i,j,k)*dudx(i,j,k) + bma(i,j,k)*dudy(i,j,k) + dma(i,j,k)*dudz(i,j,k))*dt
	vint(i,j,k) = v(i,j,k) - (amb(i,j,k)*dvdx(i,j,k) + v(i,j,k)*dvdy(i,j,k) + dmb(i,j,k)*dvdz(i,j,k))*dt
	wint(i,j,k) = w(i,j,k) - (amd(i,j,k)*dwdx(i,j,k) + bmd(i,j,k)*dwdy(i,j,k) + w(i,j,k)*dwdz(i,j,k))*dt

	enddo
	enddo
	enddo

	i = nx1
	do k = 1, nz
	do j = 1, ny

	dudx(i,j,k) = (u(i+1,j,k) - u(i-1,j,k)) / (2*dx)
	dudy(i,j,k) = (u(i,j+1,k) - u(i,j-1,k)) / (2*dy)
	dudz(i,j,k) = (u(i,j,k+1) - u(i,j,k-1)) / (2*dzx(i,j,k))

	bma(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j+1,k) + v(i-1,j+1,k)) * 0.25
	dma(i,j,k) = (w(i,j,k) + w(i-1,j,k) + w(i,j,k+1) + w(i-1,j,k+1)) * 0.25

	uint(i,j,k) = u(i,j,k) - (u(i,j,k)*dudx(i,j,k) + bma(i,j,k)*dudy(i,j,k) + dma(i,j,k)*dudz(i,j,k))*dt

	enddo
	enddo


	j = ny1
	do k = 1, nz
	do i = 1, nx

	dvdx(i,j,k) = (v(i+1,j,k) - v(i-1,j,k)) / (2*dx)
	dvdy(i,j,k) = (v(i,j+1,k) - v(i,j-1,k)) / (2*dy)
	dvdz(i,j,k) = (v(i,j,k+1) - v(i,j,k-1)) / (2*dzy(i,j,k))

	amb(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j-1,k) + u(i+1,j-1,k)) * 0.25
	dmb(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i,j,k+1) + w(i,j-1,k+1)) * 0.25

	vint(i,j,k) = v(i,j,k) - (amb(i,j,k)*dvdx(i,j,k) + v(i,j,k)*dvdy(i,j,k) + dmb(i,j,k)*dvdz(i,j,k))*dt

	enddo
	enddo


	k = nz1
	do j = 1, ny
	do i = 1, nx

	dwdx(i,j,k) = (w(i+1,j,k) - w(i-1,j,k)) / (2*dx)
	dwdy(i,j,k) = (w(i,j+1,k) - w(i,j-1,k)) / (2*dy)
	dwdz(i,j,k) = (w(i,j,k+1) - w(i,j,k-1)) / (2*dzz(i,j,k))

	amd(i,j,k) = (u(i,j,k) + u(i+1,j,k) + u(i,j,k-1) + u(i+1,j,k-1)) * 0.25
	bmd(i,j,k) = (v(i,j,k) + v(i,j+1,k) + v(i,j,k-1) + v(i,j+1,k-1)) * 0.25

	wint(i,j,k) = w(i,j,k) - (amd(i,j,k)*dwdx(i,j,k) + bmd(i,j,k)*dwdy(i,j,k) + w(i,j,k)*dwdz(i,j,k))*dt

	enddo
	enddo

	endif


END SUBROUTINE classico



!==================================================================================================================


SUBROUTINE rotacional(uint,vint,wint)

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	! auxiliares de velocidades: velocidades lagrangianas
	real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: uint
	real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: vint
	real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: wint

	! números de Courant e auxiliares de velocidades: velocidades eulerianas com contornos adicionais
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
	integer :: ntal, tt
	real(8)    :: tal

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2


	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then ! upwind 1st

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
	dudz(i,j,k) = (dp(i,j,k)*u(i,j,k)-dp(i,j,k-1)*u(i,j,k-1))/dzx(i,j,k)+(dn(i,j,k+1)*u(i,j,k+1)-dn(i,j,k)*u(i,j,k))/dzx(i,j,k)

	uint(i,j,k) = u(i,j,k) - (dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k))*dt

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
	dvdz(i,j,k) = (dp(i,j,k)*v(i,j,k)-dp(i,j,k-1)*v(i,j,k-1))/dzy(i,j,k)+(dn(i,j,k+1)*v(i,j,k+1)-dn(i,j,k)*v(i,j,k))/dzy(i,j,k)

	vint(i,j,k) = v(i,j,k) - (dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k))*dt

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
	dwdz(i,j,k) = (dp(i,j,k)*w(i,j,k)-dp(i,j,k-1)*w(i,j,k-1))/dzz(i,j,k)+(dn(i,j,k+1)*w(i,j,k+1)-dn(i,j,k)*w(i,j,k))/dzz(i,j,k)

	wint(i,j,k) = w(i,j,k) - (dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k))*dt

	enddo
	enddo
	enddo

	else

	write(*,*) "não possui outro esquema que não seja upwind 1st"
	STOP

	endif

END SUBROUTINE rotacional


!==================================================================================================================


SUBROUTINE antissim(uint,vint,wint)

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	! auxiliares de velocidades: velocidades lagrangianas
	real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: uint
	real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: vint
	real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: wint

	! números de Courant e auxiliares de velocidades: velocidades eulerianas com contornos adicionais
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
	integer :: ntal, tt
	real(8)    :: tal

	!plotagem
	real(8) :: acont, bcont, dcont 
	integer :: loca(3), locb(3), locd(3)

	!auxiliares
	real(8) :: aux1, aux2


	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	if (der == 1) then ! upwind 1st

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

	dudz(i,j,k) = ((dp(i,j,k)*u(i,j,k)-dp(i,j,k-1)*u(i,j,k-1))/dzx(i,j,k)+(dn(i,j,k+1)*u(i,j,k+1)-dn(i,j,k)*u(i,j,k))/dzx(i,j,k) &
	+ dp(i,j,k)*(u(i,j,k)-u(i,j,k-1))/dzx(i,j,k) + dn(i,j,k)*(u(i,j,k+1)-u(i,j,k))/dzx(i,j,k)) * 0.5


	uint(i,j,k) = u(i,j,k) - (dudx(i,j,k) + dudy(i,j,k) + dudz(i,j,k))*dt

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


	dvdz(i,j,k) = ((dp(i,j,k)*v(i,j,k)-dp(i,j,k-1)*v(i,j,k-1))/dzy(i,j,k)+(dn(i,j,k+1)*v(i,j,k+1)-dn(i,j,k)*v(i,j,k))/dzy(i,j,k) &
	+ dp(i,j,k)*(v(i,j,k)-v(i,j,k-1))/dzy(i,j,k) + dn(i,j,k)*(v(i,j,k+1)-v(i,j,k))/dzy(i,j,k)) * 0.5


	vint(i,j,k) = v(i,j,k) - (dvdx(i,j,k) + dvdy(i,j,k) + dvdz(i,j,k))*dt

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

	dwdz(i,j,k) = ((dp(i,j,k)*w(i,j,k)-dp(i,j,k-1)*w(i,j,k-1))/dzz(i,j,k)+(dn(i,j,k+1)*w(i,j,k+1)-dn(i,j,k)*w(i,j,k))/dzz(i,j,k) &
	+ dp(i,j,k)*(w(i,j,k)-w(i,j,k-1))/dzz(i,j,k) + dn(i,j,k)*(w(i,j,k+1)-w(i,j,k))/dzz(i,j,k)) * 0.5


	wint(i,j,k) = w(i,j,k) - (dwdx(i,j,k) + dwdy(i,j,k) + dwdz(i,j,k))*dt


	enddo
	enddo
	enddo

	else

	write(*,*) "não possui outro esquema que não seja upwind 1st para antissimétrico"
	STOP

	endif

END SUBROUTINE antissim
