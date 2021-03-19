!Subrotina definir as condições de contorno das velocidades e suas influências
! Referência: Gotoh, 2013

!!! Implementação 15/04/2014
! Leonardo Romero Monteiro

!!! Últimas Modificações
! Leonardo Romero Monteiro - 10/08/2020


!subrotina para plotagem inicial
SUBROUTINE plot_i()

	USE dzs
	USE velpre
	USE tempo
	USE smag
	USE obst
	USE mms_m
	USE cond

	IMPLICIT NONE
	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA

	real(8), dimension(nx1,ny1,nz1) :: uaux, vaux, waux, uaux1, vaux1, waux1, x1, y1, z1
	real(8), dimension(nx1,ny,nz)    :: dudy, dudz
	real(8), dimension(nx,ny1,nz)    :: dvdx, dvdz
	real(8), dimension(nx,ny,nz1)    :: dwdx, dwdy
	real(8), dimension(nx,ny,nz)    :: prdaux, kaux, vorti, vortj, vortk
	real(8), dimension(nx,ny,nz)    :: xnuta,ynuta,znuta
	real(8), dimension(nx1,ny,nz1)  :: auxy
	real(8), dimension(nx1,ny1,nz)  :: auxz
	real(8), dimension(nx1,ny1)     :: etaaux
	integer :: i, j, k


!###! TERMINAL INICIAL !###!

write(*,*) "nx = ", nx,"ny = ", ny, "nz = ", nz
write(*,*) "dx = ", dx,"dy = ", dy, "dz = ", dz1
write(*,*) "dt = ", dt,"dt_frame = ", dt_frame, "ts = ", ts
write(*,*) "der = ", der,"advectivo = ", adv_type, "modelo de turb. = ", m_turb
write(*,*) "obst_t = ", obst_t,"wave_t = ", wave_t, "mms_t  = ", mms_t 
write(*,*) "ccx0 = ", ccx0 ,"ccxf = ", ccxf, "ccy0 = ", ccy0, "ccyf = ", ccyf, "ccz0 = ", ccz0, "cczf = ", cczf


!###! PLOTAGENS CONDIÇÂO INICIAL !###!
call sonda()


!###! PLOTAGENS !###!
!conservação de massa ************************************************
open (unit=100002, action= 'write', file= 'dados//conservacao_massa.txt', status= 'unknown')
write(100002,*) "it*dt", " ",  "consmax (%) ", " ", "div"
!*********************************************************************

!contagem temporal ***************************************************
open (unit=200000, file='dados//contagem')
write(200000,*) "it", " ", "it*dt", " ", "hoje", " ", "agora"
call date_and_time(values = agora)
agora1 = agora
write(200000,*) it, it*dt, agora
!*********************************************************************

! PLOTAGENS DE PLANOS (ESCOLHER UM POR SIMULAÇÃO - FAZER ALTERAÇÕES DEPENDENDO DO PLANO A SER EXTRAÍDO)
cont = 0

! interpolação das variáveis para eles ficare nos vértices dos elementos.
call interp2d_fp(eta1,nx,ny,etaaux)
call interpx_fp(u,nx1,ny,nz,uaux)
call interpy_fp(v,nx,ny1,nz,vaux)
call interpz_fp(w,nx,ny,nz1,waux)


do j = 1, ny1
do i = 1, nx1
	x1(i,j,:) = dx*(i-1)
	y1(i,j,:) = dy*(j-1)
enddo
enddo

z1(:,:,1) = 0.

do k = 2, nz1
do j = 1, ny1
do i = 1, nx1
	z1(i,j,k) = z1(i,j,k-1) + (dzx(i,j-1,k-1)+dzx(i,j,k-1))*0.5
enddo
enddo
enddo

! definição do objeti de fundo
kaux = 0
if (obst_t .ne. 0) then
	do j = 1, ny
	do i = 1, nx
	do k = 1, kw(i,j)-1 !-1 prq o kw está posicionado na aresta inferior e não no centro da célula
         !print *, ku(i,j)
         kaux(i,j,k)=1.
	enddo
	enddo
	enddo
endif


	! cálculo da vorticidade
	CALL derivax(v,nx,ny1,nz,dvdx)
	CALL derivax(w,nx,ny,nz1,dwdx)

	CALL derivay(u,nx1,ny,nz,dudy)
	CALL derivay(w,nx,ny,nz1,dwdy)

	CALL derivaz_x(u,nx1,ny,nz,dudz)
	CALL derivaz_y(v,nx,ny1,nz,dvdz)



	CALL interpy_fc(dvdz,nx,ny1,nz,ynuta)
	CALL interpz_fc(dwdy,nx,ny,nz1,znuta)
	vorti = znuta - ynuta
	
	CALL interpx_fc(dudz,nx1,ny,nz,xnuta)
	CALL interpz_fc(dwdx,nx,ny,nz1,znuta)
	vortj = xnuta - znuta
	
	
	CALL interpx_fc(dudy,nx1,ny,nz,xnuta)
	CALL interpy_fc(dvdx,nx,ny1,nz,ynuta)
	vortk = ynuta - xnuta

CALL interpx_fc(xnut,nx1,ny,nz,xnuta)
CALL interpy_fc(ynut,nx,ny1,nz,ynuta)
CALL interpz_fc(znut,nx,ny,nz1,znuta)

! apenas para arrumar a variável para plotagem
prdaux(1:nx,1:ny,1:nz) = prd1(1:nx,1:ny,1:nz)

! a plotagem para o MMS se altera
if (mms_t > 0) then
	CALL interpy_cf(erro_u,nx1,ny,nz,auxz)
	CALL interpz_cf(auxz,nx1,ny1,nz,uaux1)

	CALL interpx_cf(erro_v,nx,ny1,nz,auxz)
	CALL interpz_cf(auxz,nx1,ny1,nz,vaux1)

	CALL interpx_cf(erro_w,nx,ny,nz1,auxy)
	CALL interpy_cf(auxy,nx1,ny,nz1,waux1)
	
!*********************************************************************
open (unit=10, file='arquivos//campos_00000',form='unformatted',status='unknown')
write(10) real(x1,4),real(y1,4),real(z1,4),real(uaux,4), & 
real(vaux,4),real(waux,4),real(uaux1,4),real(vaux1,4),real(waux1,4)
 close (unit=10)
!*********************************************************************

else

!*********************************************************************
open (unit=10, file='arquivos//campos_00000',form='unformatted',status='unknown')
write(10) real(x1,4),real(y1,4),real(z1,4),real(uaux,4), & 
real(vaux,4),real(waux,4),real(prdaux,4), &
real(vorti,4),real(vortj,4),real(vortk,4),real(kaux,4),real(xnuta,4),real(ynuta,4),real(znuta,4)
 close (unit=10)
!*********************************************************************

endif

!! chama paraview
CALL visu ()
!!

END SUBROUTINE plot_i


!###################################################################################
!###################################################################################


!subrotina para plotagem ao longo do tempo
SUBROUTINE plot_f()

	USE dzs
	USE velpre
	USE tempo
	USE smag
	USE obst
	USE analis
	USE mms_m
	
	IMPLICIT NONE
	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA

	real(8), dimension(nx1,ny1,nz1) :: uaux, vaux, waux, x1, y1, z1, uaux1, vaux1, waux1
	
	
	real(8), dimension(nx1,ny,nz)    :: dudy, dudz
	real(8), dimension(nx,ny1,nz)    :: dvdx, dvdz
	real(8), dimension(nx,ny,nz1)    :: dwdx, dwdy
	
	real(8), dimension(nx,ny,nz)    :: prdaux, kaux, vorti, vortj, vortk
	real(8), dimension(nx,ny,nz)    :: xnuta,ynuta,znuta
	real(8), dimension(nx1,ny1)     :: etaaux
	real(8), dimension(nx1,ny,nz1)  :: auxy
	real(8), dimension(nx1,ny1,nz)  :: auxz
	integer ifile, nfil, i, j, k
	real(8) :: u_m, v_m, w_m, q
	! número do arquivo de saída
	integer dig1, dig2, dig3, dig4, dig5

	! nome do arquivo de saída
	 character(5) chits


! Cálculo para o a estimativa do tempo restante
ciclo = (agora(5)-agora1(5)) * 60 * 60 + (agora(6)-agora1(6)) * 60 + (agora(7)-agora1(7)) + real(agora(8)-agora1(8))/1000
prev = (prev*6 + (ts-it)*ciclo*1./(60.*60.))/(7.)
agora1 = agora
call date_and_time(values = agora)

call sonda()

nfil=it
if(mod(it, ceiling(dt_frame/dt)).eq.0) then

	cont = cont + 1

	! interpolação das variáveis para eles ficare nos vértices dos elementos.
	call interp2d_fp(eta1,nx,ny,etaaux)
	call interpx_fp(u,nx1,ny,nz,uaux)
	call interpy_fp(v,nx,ny1,nz,vaux)
	call interpz_fp(w,nx,ny,nz1,waux)


	do j = 1, ny1
		y1(:,j,:) = dy*(j-1)
	enddo

	do i = 1, nx1
		x1(i,:,:) = dx*(i-1)
	enddo
	
	z1(:,:,1) = 0.

	do k = 2, nz1
	do j = 1, ny1
	do i = 1, nx1
		z1(i,j,k) = z1(i,j,k-1) + (dzx(i,j-1,k-1)+dzx(i,j,k-1))*0.5
	enddo
	enddo
	enddo

	! definição do objeti de fundo
	kaux = 0
	if (obst_t .ne. 0) then
		do j = 1, ny
		do i = 1, nx
		do k = 1, kw(i,j)-1!-1 prq o kw está posicionado na aresta inferior e não no centro da célula
			!print *, ku(i,j)
			kaux(i,j,k)=1.
		enddo
		enddo
		enddo
	endif


	! cálculo da vorticidade
	CALL derivax(v,nx,ny1,nz,dvdx)
	CALL derivax(w,nx,ny,nz1,dwdx)

	CALL derivay(u,nx1,ny,nz,dudy)
	CALL derivay(w,nx,ny,nz1,dwdy)

	CALL derivaz_x(u,nx1,ny,nz,dudz)
	CALL derivaz_y(v,nx,ny1,nz,dvdz)



	CALL interpy_fc(dvdz,nx,ny1,nz,ynuta)
	CALL interpz_fc(dwdy,nx,ny,nz1,znuta)
	vorti = znuta - ynuta
	
	CALL interpx_fc(dudz,nx1,ny,nz,xnuta)
	CALL interpz_fc(dwdx,nx,ny,nz1,znuta)
	vortj = xnuta - znuta
	
	CALL interpx_fc(dudy,nx1,ny,nz,xnuta)
	CALL interpy_fc(dvdx,nx,ny1,nz,ynuta)
	vortk = ynuta - xnuta

	CALL interpx_fc(xnut,nx1,ny,nz,xnuta)
	CALL interpy_fc(ynut,nx,ny1,nz,ynuta)
	CALL interpz_fc(znut,nx,ny,nz1,znuta)


	! definição do número que vai na frente do arquivo
	do ifile = 1, cont
		dig1 =    ifile/10000 + 48
		dig2 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
		dig3 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
		dig4 = ( ifile - 100*( ifile/100 ) )/10 + 48
		dig5 = ( ifile - 10*( ifile/10 ) )/1 + 48
		chits(1:5) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//char(dig5)
	enddo

	! a plotagem para o MMS se altera
	if (mms_t > 0) then
		CALL interpy_cf(erro_u,nx1,ny,nz,auxz)
		CALL interpz_cf(auxz,nx1,ny1,nz,uaux1)

		CALL interpx_cf(erro_v,nx,ny1,nz,auxz)
		CALL interpz_cf(auxz,nx1,ny1,nz,vaux1)

		CALL interpx_cf(erro_w,nx,ny,nz1,auxy)
		CALL interpy_cf(auxy,nx1,ny,nz1,waux1)
			
		!*********************************************************************
		open (unit=10, action= 'write', file= 'arquivos//campos_'//chits,form='unformatted',status='unknown')
		write(10) real(x1,4),real(y1,4),real(z1,4),real(uaux,4), & 
		real(vaux,4),real(waux,4),real(uaux1,4),real(vaux1,4),real(waux1,4)
		 close (unit=10)
		!*********************************************************************

	!prdaux = sqrt(erro_p)
	else

		prdaux(1:nx,1:ny,1:nz) = prd1(1:nx,1:ny,1:nz)
		
		open (unit=10, action= 'write', file= 'arquivos//campos_'//chits,form='unformatted',status='unknown')
		write(10) real(x1,4),real(y1,4),real(z1,4),real(uaux,4), & 
		real(vaux,4),real(waux,4),real(prdaux,4), &
		real(vorti,4),real(vortj,4),real(vortk,4),real(kaux,4),real(xnuta,4),real(ynuta,4),real(znuta,4)
		close (unit=10)
	endif


	!contagem temporal ***************************************************
	write(*,*) "it,", " ", "it*dt,", " ", "ciclo,", " ","tempo restante aproximado (horas),", " ","tempo restante aproximado (min)."
	write(*,*) it, it*dt, ciclo, prev, prev*60.
	write(200000,*) it, it*dt, agora
	!*********************************************************************
endif


CALL est()

END SUBROUTINE plot_f

!*************************************************************************************
!=====================================================================================

!subrotina que apresenta dados finais da simulação (Daniel)
SUBROUTINE plot_atrib()

	USE velpre
	USE obst

	open(unit=9, action= 'write', file= 'dados//atributos.txt', status= 'unknown')
	write(9,*) "d_min = ", d_min
	write(9,*) "d_max = ", d_max
	write(9,*) "u inicial = ", uinicial
	write(9,*) "dx = ", dx, ", dy = ", dy, ", dz1 = ", dz1, "dz topo = ", dztopo
	write(9,*) "nx = ", nx, "ny = ", ny, "nz = ", nz
	write(9,*) "dt = ", dt, "ts = ", ts
	write(9,*) "amp = ", amp, "comp = ", comp
	write(9,*) "Lz dep= ", (nz-1-elev)*dz1+dztopo
	write(9,*) "elev = ", elev
	close (unit=9)

END SUBROUTINE plot_atrib



!###################################################################################

! subrotina que calcula estabilidade e conservação de volume
SUBROUTINE est()

	use dzs
	USE velpre

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA
	real(8), dimension(nx,ny,nz) :: div

	!contadores
	integer :: i, j, k

	!!! Teste de Estabilidade !!!
	! estabilidade: conservação de massa perante o desnível, conservação de massa máxima, divergente máximo (s-1) e volume total (m³)
	real(8) consmax, divmax, vol

	! estabilidade: divergente, conservação de massa e conservação de massa acumulado 
	real(8), dimension(nx,ny,nz) :: consmas, consmasac 

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	! Cálculo do divergente
	! falta produzir
	!div = ...
	div = 0.
	
	!variação do volume mássico por iteração
	consmas = 0.
	consmasac = 0.

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		consmas(i,j,k) = (u(i+1,j,k)*dy*dzx(i+1,j,k)-u(i,j,k)*dy*dzx(i,j,k) &
		+ v(i,j+1,k)*dx*dzy(i,j+1,k)-v(i,j,k)*dx*dzy(i,j,k) + (w(i,j,k+1)-w(i,j,k))*dx*dy)*dt*998.
		consmasac(i,j,k) = consmasac(i,j,k) + consmas(i,j,k)
		!volume mássico inicial
		vol = dx*dy*dz(i,j,k)*998. + vol
	enddo
	enddo
	enddo
	if (consmax < maxval(abs(consmasac))/vol) then	
		consmax = maxval(abs(consmasac))/vol
	endif

	if (divmax < maxval(div)) then
		divmax = maxval(div)
	endif

	!===============================================================================================================
	! Plotagem

	!conservação de massa ************************************************
	!open (unit=100002, action= 'write', file= 'conservacao_massa.txt', status= 'unknown')
	write(100002,*) it*dt, consmax, divmax

	if(mod(it, ceiling(dt_frame/dt)).eq.0) then
		write(*,*) " "
 		write(*,*) "**div. max", divmax, "***dz min.", minval(dz)
		write(*,*) " "
	endif

END SUBROUTINE est


!==================================================================================================================
!==================================================================================================================


! subrotina para a criação de sondas 
SUBROUTINE sonda()

	use dzs
	USE velpre

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA

	real(8), dimension(nx,ny) :: eta_ana 
	real(8) :: error_eta, x, y
	!contadores
	integer :: i, j, k
!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================
	error_eta = 0.
	
if (it == 0) then 
	!medição da superfície livre ************************************************
	open (unit=100003, action= 'write', file= 'dados//sonda1.txt', status= 'unknown')
	write(100003,*) "it*dt", " ",  "eta (m) "
	i = 1; j = 1
	write(100003,*) it*dt, eta1(i,j)
	!*********************************************************************

else

	!medição da superfície livre ************************************************
	i = 1; j= 1
	write(100003,*) it*dt, eta1(i,j)
	!*********************************************************************
endif

if (it == ts) then
	do j = 1, ny; y = (real(j)-0.5)*dy
	do i = 1, nx; x = (real(i)-0.5)*dx
		eta_ana(i,j) = 0.1*cos(x*pi/(dx*nx))*cos(y*pi/(dy*ny))*cos(2.*pi*it*dt/3.01)
		error_eta = error_eta + abs(eta_ana(i,j) - eta1(i,j))
	enddo
	enddo

	error_eta = error_eta/(nx*ny)

	write(*,*) "error_eta", error_eta
endif

!==================================================================================================================
END SUBROUTINE sonda






