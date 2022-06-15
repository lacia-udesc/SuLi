!Subrotina que implementa o método level set
!Referenciado por seção

!Implementação 06/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro 

SUBROUTINE level_set_ini()

USE ls_param
IMPLICIT NONE

!Declarado também no programa
integer :: i, j, k, tipo, ii, jj, ilamb
real(8), dimension(nx,ny,nz) :: dist_sign, distx, disty, distz
real(8), dimension(nx,ny) :: b
real(8), dimension(nx) :: x
real(8), dimension(ny) :: y
real(8), dimension(nz) :: z
real(8) :: dist, ampl, lambdax, lambday, prof, lsaux, m, aux1, xaux1, xaux2, xaux3, xaux4, erro1, erro2, erro3, erro4

!Inicialização das variáveis
alpha1 = 1.5 !Número de células que varão parte da espessura, pois é uma função suave
rho_f1 = 1.204!kg/m³ ar (ls negativo) 20°C
mi_f1 = 0.000018253 !Pa/s  !0.00001516!m²/s ar (ls negativo) 20°C !no incompact3d tá como NI
rho_f2 =998.0!kg/m³ água saturada (ls positivo) 20°C.
mi_f2 = 0.00100798  !Pa/s  !0.00000101!m²/s água saturada (ls positivo) 20°C !no incompact3d tá como NI
sigma = 0.0728!N/m tensão superficial da água 20°C
rho_m = abs(rho_f2-rho_f1)*0.5

!Coeficientes de integração RK3 TVD
adtl(1)=1.
bdtl(1)=0.
gdtl(1)=1.

adtl(2)=3./4.
bdtl(2)=1./4.
gdtl(2)=1./4.

adtl(3)=1./3.
bdtl(3)=2./3.
gdtl(3)=2./3.

!Resolução do problema
if (t_hs == 0) then
  ls_m = 0.
else
  ls_m = 0.481803535*dx1*alpha1
endif

do i = 1, nx
x(i) = (i-0.5) * dx
enddo

do j = 1, ny
y(j) = (j-0.5) * dy
enddo

do k = 1, nz
z(k) = (k-0.5) * dz
enddo

!Tipos
!1 = Onda
!2 = Onda 2 
!3 = Barragem
!4 = Gota
!5 = fennema1990 ou aureli2008
!6 = MMS
!7 = Koshizuka et al., 1995

!Definindo o tipo:
tipo = 1

if (tipo == 1) then
	!CALL waves_coef()
	!Condição inicial de onda
	ampl = 0. 	 !Amplitude da onda
	lambdax = 2. !Comprimento da onda na direção x
	lambday = 2. !Wave length
	prof = 0.5   !Profundidade do escoamento sem a onda

	!Calcula efetivamente o ls
	if (lambday .ne. 0.) then
	    do k = 1, nz
	    do j = 1, ny
	    do i = 1, nx
		distz(i,j,k) = cos(2.*pi/lambday *y(j)) !Variação em y
	    enddo
	    enddo
	    enddo
    else
		distz = 1.
    endif

    do k = 1, nz
    do j = 1, ny
    do i = 1, nx
	distz(i,j,k) = -ampl * distz(i,j,k) * cos(2.*pi/lambdax *x(i)) +  prof !Variação em x
    enddo
    enddo
    enddo

	!Melhorar
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
    ls(i,j,k) =  distz(i,j,k) - z(k)
    !dist_sign(i,j,k) =  distz1(i,j,k) - z1(i,j,k) !function sign
	enddo
	enddo
	enddo
elseif (tipo == 2) then
	CALL waves_coef()
	!CALL waves()
elseif (tipo == 3) then
	!Condição inicial de barragem
	prof    = 0.0 !5715 !dam depth
	lambdax = 0.0 !5715 !dam x-length
	lambday = 0.0 	    !dam y-length
	ampl = 0.05715
	m = 30
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	lsaux = (z(k)-(prof))**m  + (x(i) - (lambdax))**m + (y(j) - lambday )**m/10**m
	ls(i,j,k) = -lsaux**(1./m)
	ls(i,j,k) = ls(i,j,k) + ampl 
	enddo
	enddo
	enddo
elseif (tipo == 4) then
	!Condição inicial de bolha quadrada/elipsoidal
	m = 2
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	lsaux = (x(i)-0.03)**m/2. + (y(j)-0.03)**m/2.  + (z(k)-0.01)**m/2.
	ls(i,j,k) = -lsaux**(1./m)
	ls(i,j,k) = ls(i,j,k) + 0.005 
	enddo
	enddo
	enddo
elseif (tipo == 5) then
	!Condição inicial de barragem
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	if ((x(i) >= 0.8) .and. (x(i) <= 0.825)) then
		b(i,j) = 0.15*sin((x(i)-0.8)*pi/0.05)  !5.*sin((x(i)-90)*pi/20.) !
	!if ((x(i) >= 90.) .and. (x(i) <= 100.)) then
	!b(i,j) = 10.*sin((x(i)-90)*pi/20.)  !5.*sin((x(i)-90)*pi/20.) !
	else
		b(i,j) = 0.
	endif
	enddo
	enddo
	enddo
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	if (x(i) <= 0.825) then
		ls(i,j,k) = 0.15 -b(i,j) - z(k)
	else
		ls(i,j,k) = 0. -z(k) !0.01 - z(k) !
	endif
	!if (x(i) <= 100) then
		!ls(i,j,k) = 10.-b(i,j) - z(k)
	!else
		!ls(i,j,k) =  0. - z(k) !5. - z(k) !
	!endif
	enddo
	enddo
	enddo
elseif (tipo == 6) then !MMS
	CALL mms_i()
elseif (tipo == 7) then
	!Condição inicial de barragem
	prof    = 0.0!5715 ! dam depth
	lambdax = 3.22!5715 ! dam x-length
	lambday = 0.0 ! dam y-length
	ampl = 1.22
	m = 30
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	lsaux = (z(k)-prof)**m/(0.450819672**m)  + (-x(i) +lambdax)**m
	ls(i,j,k) = -lsaux**(1./m)
	ls(i,j,k) = ls(i,j,k) + ampl 
	enddo
	enddo
	enddo
endif

do k = 1, nz
do j = 1, ny
do i = 1, nx
ls(i,j,k) = ls(i,j,k) + ls_m !corrigir por causa na nova HS
enddo
enddo
enddo

CALL heaviside()


vol_ini = 0.
do k = 1, nz
do j = 1, ny
do i = 1, nx
vol_ini = vol_ini + (1.-(1.-hs(i,j,k)))*dx*dy*dz
enddo
enddo
enddo

CALL mod_ls1()

END SUBROUTINE level_set_ini

SUBROUTINE level_set()

USE ls_param
USE velpre

IMPLICIT NONE
real(8),dimension(nx,ny,nz) :: sy7_ls,gx_ls,ta1_ls,sy7_ls1,gx_ls1,ta1_ls1
integer :: i, j, k, itrl
real(8) :: aux1, aux2, dtaux

dtaux = dt1
dt1 = dt

do itrl=1,3
CALL conv_weno(sy7_ls)
CALL intt_ls(sy7_ls,gx_ls,ta1_ls,itrl,ls)
enddo

dt1 = dtaux

CALL reinic_weno(sy7_ls1,gx_ls1,ta1_ls1)

CALL mod_ls1()
CALL heaviside()

vol_ins = 0.
do k = 1, nz
do j = 1, ny
do i = 1, nx
vol_ins = vol_ins + (1.-(1.-hs(i,j,k)))*dx*dy*dz
enddo
enddo
enddo

!if ((it == 1))) vol_ini = vol_ins

do while ((abs(vol_ins-vol_ini)/vol_ins > 0.1) .and. (mms_t .ne. 2)) !Erro aceitável de 1 % para conservação de volume e se não tiver obstáculos

do k = 1, nz
do j = 1, ny
do i = 1, nx
ls(i,j,k) = ls(i,j,k) - dt1 * (vol_ins-vol_ini)*mod_ls(i,j,k)/vol_ini
enddo
enddo
enddo

CALL mod_ls1()
CALL heaviside()

vol_ins = 0.
do k = 1, nz
do j = 1, ny
do i = 1, nx
vol_ins = vol_ins + (1.-(1.-hs(i,j,k)))*dx*dy*dz
enddo
enddo
enddo
enddo

CALL mod_ls1() !Função para plotagem e cálculo da curvatura

CALL heaviside()


END SUBROUTINE level_set

SUBROUTINE intt_ls(hx,gx,ta1,itrl,ls1)

USE ls_param

implicit none

integer :: i,j,k,itrl
real(8),dimension(nx,ny,nz) :: hx,gx,ta1,ls1

!RK3 TVD
if (adtl(itrl)==1.) then
	!Precisa ser ls0 para todos os subtempos
	ta1 = ls1
endif

gx = ls1 !ls1 e ls2 para os subtempos
ls1 = adtl(itrl)*ta1+bdtl(itrl)*gx+gdtl(itrl)*dt1*hx


END SUBROUTINE intt_ls

SUBROUTINE conv_weno(sy7)

USE ls_param
USE velpre

implicit none

integer :: i,j,k,ihs

real(8),dimension(nx,ny,nz) :: ta1,tb1,tc1,td1,te1,tf1,sy7
real(8) :: apos, aneg, bpos, bneg, cpos, cneg


!ihs = 1
CALL der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs)

do k = 1, nz
do j = 1, ny
do i = 1, nx

   apos = max((u(i,j,k)+u(i+1,j,k))*0.5,0.)
   aneg = min((u(i,j,k)+u(i+1,j,k))*0.5,0.)
   bpos = max((v(i,j,k)+v(i,j+1,k))*0.5,0.)
   bneg = min((v(i,j,k)+v(i,j+1,k))*0.5,0.)
   cpos = max((w(i,j,k)+w(i,j,k+1))*0.5,0.)
   cneg = min((w(i,j,k)+w(i,j,k+1))*0.5,0.)

   sy7(i,j,k) = -(apos*td1(i,j,k) + aneg*ta1(i,j,k) + bpos*te1(i,j,k) + bneg*tb1(i,j,k) + cpos*tf1(i,j,k) + cneg*tc1(i,j,k))

enddo
enddo
enddo

end subroutine conv_weno

SUBROUTINE der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs)

USE disc
USE cond

implicit none

integer :: i,j,k,ihs
real(8),dimension(nx,ny,nz) :: ls,ta1,tb1,tc1,td1,te1,tf1

if (ccx0 == 0) then
	ihs = 0
else
	ihs = 2
endif

do k = 1, nz
do j = 1, ny
call weno1(ta1(:,j,k),td1(:,j,k),nx,dx,ls(:,j,k),ihs)
enddo
enddo

if (ccy0 == 0) then
	ihs = 0
else
	ihs = 2
endif

do k = 1, nz
do i = 1, nx
call weno1(tb1(i,:,k),te1(i,:,k),ny,dy,ls(i,:,k),ihs)
enddo
enddo

ihs = 1
do j = 1, ny
do i = 1, nx
call weno1(tc1(i,j,:),tf1(i,j,:),nz,dz,ls(i,j,:),ihs)
enddo
enddo

END SUBROUTINE der_weno

SUBROUTINE reinic_weno(sy7_ls1,gx_ls1,ta1_ls1)

USE ls_param
IMPLICIT NONE
real(8),dimension(nx,ny,nz) :: sy1,sy4,func_s,ddd,ta1,tb1,tc1,td1,te1,tf1
real(8),dimension(nx,ny,nz) :: sy7_ls1,gx_ls1,ta1_ls1,lsaux,ls0
real(8) :: error
integer :: i, j, k, l, il, nr,ihs,itrl
real(8) :: mod_ls1, aux1, aux2

ls0 = ls
l = 3

!ihs = 1
CALL der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs)

do k = 1, nz
do j = 1, ny
do i = 1, nx

if (ls(i,j,k) > 0. ) then
	aux1 = max(td1(i,j,k),0.00000001)
	aux2 = -min(ta1(i,j,k),0.00000001)
	ta1(i,j,k) = max(aux1, aux2)

	aux1 = max(te1(i,j,k),0.00000001)
	aux2 = -min(tb1(i,j,k),0.00000001)
	tb1(i,j,k) = max(aux1, aux2)

	aux1 = max(tf1(i,j,k),0.00000001)
	aux2 = -min(tc1(i,j,k),0.00000001)
	tc1(i,j,k) = max(aux1, aux2)
else
	aux1 = max(ta1(i,j,k),0.00000001)
	aux2 = -min(td1(i,j,k),0.00000001)
	ta1(i,j,k) = max(aux1, aux2)

	aux1 = max(tb1(i,j,k),0.00000001)
	aux2 = -min(te1(i,j,k),0.00000001)
	tb1(i,j,k) = max(aux1, aux2)

	aux1 = max(tc1(i,j,k),0.00000001)
	aux2 = -min(tf1(i,j,k),0.00000001)
	tc1(i,j,k) = max(aux1, aux2)
endif
mod_ls1 = sqrt(ta1(i,j,k)*ta1(i,j,k) + tb1(i,j,k)*tb1(i,j,k) + tc1(i,j,k)*tc1(i,j,k))
func_s(i,j,k) = ls(i,j,k) / sqrt(ls(i,j,k)*ls(i,j,k) + mod_ls1*mod_ls1*dx1*dx1)
enddo
enddo
enddo

il = 0
error = 999.

do while (error > dt1*dx1*dx1 )
il = il + 1

do itrl=1,3

CALL der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs)

do k = 1, nz
do j = 1, ny
do i = 1, nx

if (ls(i,j,k) > 0. ) then
	aux1 = max(td1(i,j,k),0.00000001)
	aux2 = -min(ta1(i,j,k),0.00000001)
	ta1(i,j,k) = max(aux1, aux2)

	aux1 = max(te1(i,j,k),0.00000001)
	aux2 = -min(tb1(i,j,k),0.00000001)
	tb1(i,j,k) = max(aux1, aux2)

	aux1 = max(tf1(i,j,k),0.00000001)
	aux2 = -min(tc1(i,j,k),0.00000001)
	tc1(i,j,k) = max(aux1, aux2)
else
	aux1 = max(ta1(i,j,k),0.00000001)
	aux2 = -min(td1(i,j,k),0.00000001)
	ta1(i,j,k) = max(aux1, aux2)

	aux1 = max(tb1(i,j,k),0.00000001)
	aux2 = -min(te1(i,j,k),0.00000001)
	tb1(i,j,k) = max(aux1, aux2)

	aux1 = max(tc1(i,j,k),0.00000001)
	aux2 = -min(tf1(i,j,k),0.00000001)
	tc1(i,j,k) = max(aux1, aux2)
endif
mod_ls1      = sqrt(ta1(i,j,k)*ta1(i,j,k) + tb1(i,j,k)*tb1(i,j,k) + tc1(i,j,k)*tc1(i,j,k))
sy7_ls1(i,j,k) = func_s(i,j,k) * (1.-mod_ls1 )
lsaux(i,j,k)  = ls(i,j,k)

enddo
enddo
enddo
CALL intt_ls(sy7_ls1,gx_ls1,ta1_ls1,itrl,ls)
enddo

if (il == l ) then !Número máximo de iterações
	error = 0.
elseif (il < 1) then !Número mínimo de iterações
	error = 999.
else
	error = 0.
	nr = 0
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	if (abs(ls(i,j,k)) <= 2.*alpha1 * dx1) then
		error = error + ls(i,j,k) - lsaux(i,j,k)
		nr = nr + 1
	endif
	enddo
	enddo
	enddo
	error = error / nr
endif
enddo

!write(*,*) il

END SUBROUTINE reinic_weno

SUBROUTINE weno1(dphidxp,dphidxn,nx1,dx1,phi0,ihs)

IMPLICIT NONE
real(8),dimension(nx1) :: phiaux,dphidxp,dphidxn,phi0
real(8),dimension(3) ::isup, isun, auxx
real(8),dimension(3) ::alpup, omgup,alpun, omgun
integer :: i,kk, ii,nx1, ihs
real(8) :: dx1, mod_phi1,aux1,aux2,aux3,aux4,aux5,aux6,aux,aux11, aux12
real(8),dimension(-2:nx1+3) :: phi1
real(8),dimension(nx1+4)  :: un
real(8),dimension(-3:nx1) :: up

aux1 = 13./12.
aux2 = 1./4.
aux3 = 1./6.
auxx(1) = 0.1
auxx(2) = 0.6
auxx(3) = 0.3
aux6 = 0.00000001
phi1(1:nx1) = phi0(1:nx1)

if (ihs == 0) then !Contorno periódico
	phi1(0)  = phi1(nx1-1)
	phi1(-1) = phi1(nx1-2)
	phi1(-2) = phi1(nx1-3)


	phi1(nx1+1) = phi1(2)
	phi1(nx1+2) = phi1(3)
	phi1(nx1+3) = phi1(4)
elseif (ihs == 1) then !Distance extrapolation
	!phi1(0)     = 2*phi1(1)     - phi1(2)
	!phi1(-1)    = 2*phi1(0)     - phi1(1)
	!phi1(-2)    = 2*phi1(-1)    - phi1(0)
	!phi1(nx1+1) = 2*phi1(nx1)   - phi1(nx1-1)
	!phi1(nx1+2) = 2*phi1(nx1+1) - phi1(nx1)
	!phi1(nx1+3) = 2*phi1(nx1+2) - phi1(nx1+1)
	phi1(0)  = 1./5. * (12.*phi1(1)  - 9.*phi1(2) + 2.*phi1(3) )
	phi1(-1) = 1./5. * (12.*phi1(0)  - 9.*phi1(1) + 2.*phi1(2) )
	phi1(-2) = 1./5. * (12.*phi1(-1) - 9.*phi1(0) + 2.*phi1(1) )
	!phi1(nx1+1) = 1./5. * (12.*phi1(nx1)   - 9.*phi1(nx1-1) + 2.*phi1(nx1-2))
	!phi1(nx1+2) = 1./5. * (12.*phi1(nx1+1) - 9.*phi1(nx1)   + 2.*phi1(nx1-1))
	!phi1(nx1+3) = 1./5. * (12.*phi1(nx1+2) - 9.*phi1(nx1+1) + 2.*phi1(nx1)  )
	phi1(nx1+1) = 1./11. * (18.*phi1(nx1)   - 9.*phi1(nx1-1) + 2.*phi1(nx1-2))
	phi1(nx1+2) = 1./11. * (18.*phi1(nx1+1) - 9.*phi1(nx1)   + 2.*phi1(nx1-1))
	phi1(nx1+3) = 1./11. * (18.*phi1(nx1+2) - 9.*phi1(nx1+1) + 2.*phi1(nx1)  )
elseif (ihs == 2) then !Dderivative zero
	!phi1(0)     = phi1(1)
	!phi1(nx1+1) = phi1(nx1)
	!phi1(-1)    = phi1(2)
	!phi1(nx1+2) = phi1(nx1-1)
	!phi1(-2)    = phi1(3)
	!phi1(nx1+3) = phi1(nx1-2)
	phi1(0)  = 1./11. * (18.*phi1(1)  - 9.*phi1(2) + 2.*phi1(3) )
	phi1(-1) = 1./11. * (18.*phi1(0)  - 9.*phi1(1) + 2.*phi1(2) )
	phi1(-2) = 1./11. * (18.*phi1(-1) - 9.*phi1(0) + 2.*phi1(1) )
	phi1(nx1+1) = 1./11. * (18.*phi1(nx1)   - 9.*phi1(nx1-1) + 2.*phi1(nx1-2))
	phi1(nx1+2) = 1./11. * (18.*phi1(nx1+1) - 9.*phi1(nx1)   + 2.*phi1(nx1-1))
	phi1(nx1+3) = 1./11. * (18.*phi1(nx1+2) - 9.*phi1(nx1+1) + 2.*phi1(nx1)  )
endif

do i=-3,nx1
up(i)=(phi1(i+3)-phi1(i+2))/dx1
enddo

do i=1,nx1+4
un(i)=(phi1(i-2)-phi1(i-3))/dx1
enddo

do i=1,nx1
isup(1) = aux1 * (up(i)-2*up(i-1)+up(i-2))*(up(i)-2*up(i-1)+up(i-2)) &
+ aux2 * (up(i)-4*up(i-1)+3*up(i-2))*(up(i)-4*up(i-1)+3*up(i-2))
isun(1) = aux1 * (un(i)-2*un(i+1)+un(i+2))*(un(i)-2*un(i+1)+un(i+2)) &
+ aux2 * (un(i)-4*un(i+1)+3*un(i+2))*(un(i)-4*un(i+1)+3*un(i+2))

isup(2) = aux1 * (up(i-1)-2*up(i-2)+up(i-3))*(up(i-1)-2*up(i-2)+up(i-3)) + aux2 * (up(i-1)-up(i-3))*(up(i-1)-up(i-3))
isun(2) = aux1 * (un(i+1)-2*un(i+2)+un(i+3))*(un(i+1)-2*un(i+2)+un(i+3)) + aux2 * (un(i+1)-un(i+3))*(un(i+1)-un(i+3))

isup(3) = aux1 * (up(i-2)-2*up(i-3)+up(i-4))*(up(i-2)-2*up(i-3)+up(i-4)) &
+ aux2 * (3*up(i-2)-4*up(i-3)+up(i-4))*(3*up(i-2)-4*up(i-3)+up(i-4))
isun(3) = aux1 * (un(i+2)-2*un(i+3)+un(i+4))*(un(i+2)-2*un(i+3)+un(i+4)) &
+ aux2 * (3*un(i+2)-4*un(i+3)+un(i+4))*(3*un(i+2)-4*un(i+3)+un(i+4))

do kk = 1, 3
alpup(kk) = auxx(kk) / ((aux6 + isup(kk))*(aux6 + isup(kk)))
alpun(kk) = auxx(kk) / ((aux6 + isun(kk))*(aux6 + isun(kk)))
enddo

do kk = 1, 3
omgup(kk) = alpup(kk) / (alpup(1)+alpup(2)+alpup(3))
omgun(kk) = alpun(kk) / (alpun(1)+alpun(2)+alpun(3))
enddo

dphidxp(i) = aux3* (omgup(1) * (2*up(i)-7*up(i-1)+11*up(i-2)) + &
omgup(2) * (-up(i-1)+5*up(i-2)+2*up(i-3)) + omgup(3) * (2*up(i-2)+5*up(i-3)-up(i-4)) )
dphidxn(i) = aux3* (omgun(1) * (2*un(i)-7*un(i+1)+11*un(i+2)) + &
omgun(2) * (-un(i+1)+5*un(i+2)+2*un(i+3)) + omgun(3) * (2*un(i+2)+5*un(i+3)-un(i+4)) )
enddo

END SUBROUTINE weno1

SUBROUTINE heaviside()
USE ls_param
USE velpre
IMPLICIT NONE

integer :: i, j, k, coefa1,ihs
real(8) :: aux1, aux2, aux3, aux4
real(8), dimension(nx,ny,nz) :: sy60, sy61,ta1,tb1,tc1,td1,te1,tf1

!ihs = 2
CALL der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs)

do k = 1, nz
do j = 1, ny
do i = 1, nx
if (abs(ta1(i,j,k)) > abs(td1(i,j,k))) then
	ta1(i,j,k) = ta1(i,j,k)
else
	ta1(i,j,k) = td1(i,j,k)
endif

if (abs(tb1(i,j,k)) > abs(te1(i,j,k))) then
	tb1(i,j,k) = tb1(i,j,k)
else
	tb1(i,j,k) = te1(i,j,k)
endif

if (abs(tc1(i,j,k)) > abs(tf1(i,j,k))) then
	tc1(i,j,k) = tc1(i,j,k)
else
	tc1(i,j,k) = tf1(i,j,k)
endif
enddo
enddo
enddo

!drhodx = 0.
!drhody = 0.
!drhodz = 0.

!dmidx = 0.
!dmidy = 0.
!dmidz = 0.

aux1 = mi_f1 * 10000.
aux2 = mi_f2 * 10000.

do k = 1, nz
do j = 1, ny
do i = 1, nx
if (ls(i,j,k) < -alpha1 * dx1) then !rever este dx
	hs(i,j,k) = 0.

	hsx(i,j,k) = 0.	
	hsy(i,j,k) = 0.	
	hsz(i,j,k) = 0.	
elseif (ls(i,j,k) > alpha1 * dx1) then
	hs(i,j,k) = 1.
	hsx(i,j,k) = 0.	
	hsy(i,j,k) = 0.	
	hsz(i,j,k) = 0.	
else
	hs(i,j,k) = 0.5*(1.+ls(i,j,k)/(alpha1*dx1) + 1./pi * sin(pi*ls(i,j,k)/(alpha1*dx1)))
	hsx(i,j,k) = 0.5*ta1(i,j,k)*(1. + cos(pi*ls(i,j,k)/(alpha1*dx1)))/(alpha1*dx1)
	hsy(i,j,k) = 0.5*tb1(i,j,k)*(1. + cos(pi*ls(i,j,k)/(alpha1*dx1)))/(alpha1*dx1)
	hsz(i,j,k) = 0.5*tc1(i,j,k)*(1. + cos(pi*ls(i,j,k)/(alpha1*dx1)))/(alpha1*dx1)
	!if (t_hs == 0) then
		!drhodx(i,j,1) =  0.5* (rho_f2-rho_f1)*sy60(i,j,1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) / (alpha1*dx1) 
		!drhody(i,j,1) =  0.5* (rho_f2-rho_f1)*sy61(i,j,1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) / (alpha1*dx1)
		!dmidx(i,j,1) =  0.5* (mi_f2-mi_f1)*sy60(i,j,1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) / (alpha1*dx1) 
		!dmidy(i,j,1) =  0.5* (mi_f2-mi_f1)*sy61(i,j,1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) / (alpha1*dx1)
	!else
		!drhodx(i,j,1) = (-rho_f1**(1.-hs(i,j,1))*log(rho_f1) + rho_f2**hs(i,j,1)*log(rho_f2)) * &
		!	        sy60(i,j,1)/(alpha1*dx1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) 
		!drhody(i,j,1) = (-rho_f1**(1.-hs(i,j,1))*log(rho_f1) + rho_f2**hs(i,j,1)*log(rho_f2)) * &
		!	        sy61(i,j,1)/(alpha1*dx1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) 

		!dmidx(i,j,1) = (-mi_f1**(1.-hs(i,j,1))*log(mi_f1) + mi_f2**hs(i,j,1)*log(mi_f2)) * &
		!	        sy60(i,j,1)/(alpha1*dx1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) 
		!dmidy(i,j,1) = (-mi_f1**(1.-hs(i,j,1))*log(mi_f1) + mi_f2**hs(i,j,1)*log(mi_f2)) * &
		!	        sy61(i,j,1)/(alpha1*dx1) * ( 1. + cos(pi*phi(i,j,1)/(alpha1*dx1)) ) 
	!endif
endif

if (t_hs == 0) then
	rho(i,j,k) = rho_f1 * (1.-hs(i,j,k)) + rho_f2 * hs(i,j,k)
	ls_nu(i,j,k) = mi_f1  * (1.-hs(i,j,k)) + mi_f2  * hs(i,j,k)
else
	rho(i,j,k) = rho_f1**(1.-hs(i,j,k)) + rho_f2**hs(i,j,k) -1.
	ls_nu(i,j,k) = (aux1** (1.-hs(i,j,k)) + aux2**hs(i,j,k) -1.)/10000.
endif
enddo
enddo
enddo

END SUBROUTINE heaviside

SUBROUTINE mod_ls1()
USE ls_param
IMPLICIT NONE

integer :: i, j, k,ihs
real(8) :: aux1, aux2
real(8), dimension(nx,ny,nz) :: ta1,tb1,tc1,td1,te1,tf1,dlsdxa,dlsdya,dlsdza
!ihs = 1
CALL der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs)

do k = 1, nz
do j = 1, ny
do i = 1, nx
if (abs(ta1(i,j,k)) > abs(td1(i,j,k))) then
	ta1(i,j,k) = ta1(i,j,k)
else
	ta1(i,j,k) = td1(i,j,k)
endif

if (abs(tb1(i,j,k)) > abs(te1(i,j,k))) then
	tb1(i,j,k) = tb1(i,j,k)
else
	tb1(i,j,k) = te1(i,j,k)
endif

if (abs(tc1(i,j,k)) > abs(tf1(i,j,k))) then
	tc1(i,j,k) = tc1(i,j,k)
else
	tc1(i,j,k) = tf1(i,j,k)
endif
enddo
enddo
enddo

do k = 1, nz
do j = 1, ny
do i = 1, nx
mod_ls(i,j,k) = sqrt(ta1(i,j,k)*ta1(i,j,k) + tb1(i,j,k)*tb1(i,j,k) + tc1(i,j,k)*tc1(i,j,k))

if (mod_ls(i,j,k) == 0) then 
	dlsdxa(i,j,k) = 0.
	dlsdya(i,j,k) = 0.
	dlsdza(i,j,k) = 0.
else
	dlsdxa(i,j,k) = ta1(i,j,k)/mod_ls(i,j,k)
	dlsdya(i,j,k) = tb1(i,j,k)/mod_ls(i,j,k)
	dlsdza(i,j,k) = tc1(i,j,k)/mod_ls(i,j,k)
endif
enddo
enddo
enddo

do k = 1, nz
do j = 1, ny
call weno1(ta1(:,j,k),td1(:,j,k),nx,dx,dlsdxa(:,j,k),ihs)
enddo
enddo

do k = 1, nz
do i = 1, nx
call weno1(tb1(i,:,k),te1(i,:,k),ny,dy,dlsdya(i,:,k),ihs)
enddo
enddo

do j = 1, ny
do i = 1, nx
call weno1(tc1(i,j,:),tf1(i,j,:),nz,dz,dlsdza(i,j,:),ihs)
enddo
enddo

do k = 1, nz
do j = 1, ny
do i = 1, nx

if (abs(td1(i,j,k)) > abs(ta1(i,j,k))) then
	ddlsdx(i,j,k) = td1(i,j,k)
else
	ddlsdx(i,j,k) = ta1(i,j,k)
endif

if (abs(te1(i,j,k)) > abs(tb1(i,j,k))) then
	ddlsdy(i,j,k) = te1(i,j,k)
else
	ddlsdy(i,j,k) = tb1(i,j,k)
endif

if (abs(tf1(i,j,k)) > abs(tc1(i,j,k))) then
	ddlsdz(i,j,k) = tf1(i,j,k)
else
	ddlsdz(i,j,k) = tc1(i,j,k)
endif
kurv(i,j,k) = ddlsdx(i,j,k) + ddlsdy(i,j,k) + ddlsdz(i,j,k)
enddo
enddo
enddo

END SUBROUTINE mod_ls1
