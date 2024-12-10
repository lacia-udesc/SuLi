!Subrotina que implementa o método level set
!Referenciado por seção

!Implementação 06/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro 

SUBROUTINE level_set_ini()

	USE ls_param
	USE cond, only: ccx0
	USE velpre, only: blx1
	
	IMPLICIT NONE

	!Declarado também no programa
	integer :: i, j, k, ii, jj, ilamb
	real(8),dimension(nx,ny,nz) :: dlsdxa,dlsdya,dlsdza
	real(8),dimension(nx,ny,nz) :: dist_sign, distx, disty, distz
	real(8),dimension(nx,ny) :: b
	real(8),save :: dist, lsaux, aux1, xaux1, xaux2, xaux3, xaux4, erro1, erro2, erro3, erro4

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

	dx1 =  max(dx,dy,dza) !(dx+dy+dz)/3. !
	dt1 =0.1 * dx1
	
	if (tipo == 1) then
		!CALL waves_coef()
		!Condição inicial de onda
		!ampl = 0.    !Amplitude da onda
		!lambdax = 2. !Comprimento da onda na direção x
		!lambday = 2. !Wave length
		!prof = 0.127 !Profundidade do escoamento sem a onda 

		!Calcula efetivamente o ls
		if (lambday .ne. 0.) then
		    do k = 1, nz
		    do j = 1, ny
		    do i = 1, nx
			distz(i,j,k) = cos(2.*pi/lambday *y(j)) !Variação em y
		    enddo
		    enddo
		    enddo
	   	elseif (lambday == 0 .and. lambdax == 0) then
			distz = 0.
		elseif (lambday == 0 .and. lambdax .ne. 0) then
			distz = 1.		
	    	endif

		if (lambdax .ne. 0.) then
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				distz(i,j,k) = -ampl * distz(i,j,k) * cos(2.*pi/lambdax *x(i)) !Variação em x
			enddo
			enddo
			enddo
	    	endif
	    
		!Melhorar
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
	  		ls(i,j,k) =  distz(i,j,k) - z(k) + prof !Variação em x
			!dist_sign(i,j,k) =  distz1(i,j,k) - z1(i,j,k) !function sign
		enddo
		enddo
		enddo
	elseif (tipo == 2) then
		CALL waves_coef()
		!CALL waves()
	elseif (tipo == 3) then
		!Condição inicial de barragem
		!prof    !dam z-length
		!lambdax !dam x-length
		!lambday !dam y-length
		!ampl    !fator de adição do comprimento em todas as direções (m) - ex: aresta do cubo
		!m       ! Curvatura do chanfro entre planos na barragem [adimensional]
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
		!ampl    ! raio da gota (ou metade de uma aresta) (m) 
		!m       ! Curvatura do chanfro entre planos da bolha (2 para redondo) [adimensional]
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
		lsaux = (x(i)-lambdax)**m/2. + (y(j)-lambday)**m/2.  + (z(k)-prof)**m/2.
		ls(i,j,k) = -lsaux**(1./m)
		ls(i,j,k) = ls(i,j,k) + ampl  
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
	!	CALL mms_i() !desatualizado
	elseif (tipo == 7) then
		!Condição inicial de barragem
		prof    = 0.!5715 ! dam depth
		lambdax = 3.22!5715 ! dam x-length
		lambday = 0.0 ! dam y-length
		ampl = 1.22
		
		m = 30
		
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
		
		lsaux = (z(k)-prof)**m/(0.450819672**m)  + (x(i) -lambdax)**m
		ls(i,j,k) = -lsaux**(1./m)
		ls(i,j,k) = ls(i,j,k) + ampl 
		
		enddo
		enddo
		enddo
	endif

	CALL heaviside()

	vol_ini = 0.
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	vol_ini = vol_ini + (1.-(1.-hs(i,j,k)))*dx*dy*dz(k)
	enddo
	enddo
	enddo

	if (ccx0==3) then
		blx1(:,:) = ls(1,:,:)
	endif

END SUBROUTINE level_set_ini

!!!####################################################################################

SUBROUTINE level_set()

	USE ls_param
	USE velpre

	IMPLICIT NONE
	
	real(8),dimension(nx,ny,nz) :: dlsdxa,dlsdya,dlsdza
	real(8),dimension(nx,ny,nz) :: sy7_ls,gx_ls,ta1_ls,sy7_ls1,gx_ls1,ta1_ls1,mod_ls
	integer :: i, j, k, itrl
	real(8),save :: aux1, aux2, dtaux

	dtaux = dt1
	dt1 = dt

	! cálculo da advecção
	do itrl=1,3
		CALL conv_weno(sy7_ls)
		CALL intt_ls(sy7_ls,gx_ls,ta1_ls,itrl,ls,nx,ny,nz)
	enddo

	dt1 = dtaux
	
	! cálculo da reinicialização
	CALL reinic_weno(ls,nx,ny,nz)

	CALL heaviside()

	! cálculo do volume "inicial - ins" antes da correção
	vol_ins = 0.
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		vol_ins = vol_ins + (1.-(1.-hs(i,j,k)))*dx*dy*dz(k)
	enddo
	enddo
	enddo
	
	! verifica se o volume está variando de acordo com o primeiro volume de todos. Tomar cuidado se volume for acrescentado no modelo, como uma adição de onda!
	do while ((abs(vol_ins-vol_ini)/vol_ins > 0.1) .and. (mms_t .ne. 2)) !Erro aceitável de 1 % para conservação de volume e se não tiver obstáculos
	
		CALL mod_ls1(ls,mod_ls,dlsdxa,dlsdya,dlsdza,nx,ny,nz) !Função para plotagem e cálculo da curvatura
		
		! correção do volume, evolui por euler explícito
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			ls(i,j,k) = ls(i,j,k) - dt1 * (vol_ins-vol_ini)*mod_ls(i,j,k)/vol_ini
		enddo
		enddo
		enddo

		CALL heaviside()
	
		! cálculo do volume "inicial - ins" antes da correção
		vol_ins = 0.
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
		vol_ins = vol_ins + (1.-(1.-hs(i,j,k)))*dx*dy*dz(k)
		enddo
		enddo
		enddo
	enddo

	CALL mod_ls1(ls,mod_ls,dlsdxa,dlsdya,dlsdza,nx,ny,nz) !Função para plotagem e cálculo da curvatura
	if (t_tens == 1) CALL curv_ls1(dlsdxa,dlsdya,dlsdza)
	
	
END SUBROUTINE level_set

!!!####################################################################################

SUBROUTINE intt_ls(hx,gx,ta1,itrl,ls1,dimx,dimy,dimz)
!rotina para fazer a integração temporal do Level-Set

	USE ls_param, only : dt1, adtl, bdtl, gdtl

	implicit none

	integer :: i,j,k,itrl,dimx,dimy,dimz
	real(8),intent(inout),dimension(dimx,dimy,dimz) :: hx,gx,ta1,ls1

	!RK3 TVD
	if (adtl(itrl)==1.) then
		!Precisa ser ls0 para todos os subtempos
		ta1 = ls1
	endif

	gx = ls1 !ls1 e ls2 para os subtempos
	ls1 = adtl(itrl)*ta1+bdtl(itrl)*gx+gdtl(itrl)*dt1*hx


END SUBROUTINE intt_ls

!!!####################################################################################

SUBROUTINE conv_weno(sy7)
!rotina para cálculo dos termos convectivos do Level-Set

	USE ls_param
	USE velpre

	implicit none

	integer :: i,j,k,ihs

	real(8),intent(out),dimension(nx,ny,nz) :: sy7
	real(8),dimension(nx,ny,nz) :: ta1,tb1,tc1,td1,te1,tf1
	real(8),save :: apos, aneg, bpos, bneg, cpos, cneg

	CALL der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs,nx,ny,nz)

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

!!!####################################################################################

SUBROUTINE der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs,dimx,dimy,dimz)
!cálculo da derivada de WENO

	USE disc, only : dx,dy,dz
	USE cond

	implicit none

	integer :: i,j,k,ihs,dimx,dimy,dimz
	real(8),intent(inout),dimension(dimx,dimy,dimz) :: ls,ta1,tb1,tc1,td1,te1,tf1

	if (ccx0 == 0) then
		ihs = 0
	else
		ihs = 2
	endif

	call wenox(ls,dimx,dimy,dimz,dx,ta1,td1,ihs)

	if (ccy0 == 0) then
		ihs = 0
	else
		ihs = 2
	endif

	call wenoy(ls,dimx,dimy,dimz,dy,tb1,te1,ihs)

    	! em cima em baixo vai ser sempre ihs 1!
	ihs = 2
	
	call wenoz(ls,dimx,dimy,dimz,dz,tc1,tf1,ihs)

END SUBROUTINE der_weno

!!!####################################################################################

SUBROUTINE reinic_weno(ls1,dimx,dimy,dimz)
!cálculo da reinicialização da função distância

	USE ls_param, only : dx1, alpha1, dt1

	IMPLICIT NONE
    
	integer :: i, j, k, l, il, nr,ihs,itrl,dimx,dimy,dimz	
	real(8),dimension(dimx,dimy,dimz) :: sy1,sy4,func_s,ddd,ta1,tb1,tc1,td1,te1,tf1,lsaux,ls0
	real(8),intent(inout),dimension(dimx,dimy,dimz) :: ls1
	real(8),dimension(dimx,dimy,dimz) :: sy7_ls1,gx_ls1,ta1_ls1
	real(8) :: error
	real(8),save :: mod_ls, aux1, aux2

	ls0 = ls1
	l = 3 ! número máximo de iterações

	CALL der_weno(ls1,ta1,tb1,tc1,td1,te1,tf1,ihs,dimx,dimy,dimz)

	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx

		if (ls1(i,j,k) > 0. ) then
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
			mod_ls = sqrt(ta1(i,j,k)*ta1(i,j,k) + tb1(i,j,k)*tb1(i,j,k) + tc1(i,j,k)*tc1(i,j,k))
			func_s(i,j,k) = ls1(i,j,k) / sqrt(ls1(i,j,k)*ls1(i,j,k) + mod_ls*mod_ls*dx1*dx1)
	enddo
	enddo
	enddo

	il = 0
	error = 999.

	do while (error > dt1*dx1*dx1 )
	il = il + 1

	! RK3 - TVD
	do itrl=1,3

		CALL der_weno(ls1,ta1,tb1,tc1,td1,te1,tf1,ihs,dimx,dimy,dimz)

		do k = 1, dimz
		do j = 1, dimy
		do i = 1, dimx

			if (ls1(i,j,k) > 0. ) then
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
    
			mod_ls      = sqrt(ta1(i,j,k)*ta1(i,j,k) + tb1(i,j,k)*tb1(i,j,k) + tc1(i,j,k)*tc1(i,j,k))
			sy7_ls1(i,j,k) = func_s(i,j,k) * (1.-mod_ls )
			lsaux(i,j,k)  = ls1(i,j,k)

		enddo
		enddo
		enddo
    
		CALL intt_ls(sy7_ls1,gx_ls1,ta1_ls1,itrl,ls1,dimx,dimy,dimz)
    
	enddo

	if (il == l ) then !Número máximo de iterações
		error = 0.
	elseif (il < 1) then !Número mínimo de iterações
		error = 999.
	else
		error = 0.
		nr = 0
		do k = 1, dimz
		do j = 1, dimy
		do i = 1, dimx
		if (abs(ls1(i,j,k)) <= 2.*alpha1 * dx1) then
			error = error + ls1(i,j,k) - lsaux(i,j,k)
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

!!!####################################################################################

SUBROUTINE heaviside()
! cálculo da função heaviside e atualização das propriedades físicas

	USE ls_param
	USE velpre
	IMPLICIT NONE

	integer :: i, j, k, coefa1,ihs
	real(8),dimension(nx,ny,nz) :: sy60, sy61,ta1,tb1,tc1,td1,te1,tf1


	if (t_tens == 1) then 
		CALL der_weno(ls,ta1,tb1,tc1,td1,te1,tf1,ihs,nx,ny,nz)

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

		hsx(i,j,k) = 0.	
		hsy(i,j,k) = 0.	
		hsz(i,j,k) = 0.	

		enddo
		enddo
		enddo
	endif

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	if (ls(i,j,k) < -alpha1 * dx1) then !rever este dx
		hs(i,j,k) = 0.

	elseif (ls(i,j,k) > alpha1 * dx1) then
		hs(i,j,k) = 1.
	else
		hs(i,j,k) = 0.5*(1.+ls(i,j,k)/(alpha1*dx1) + 1./pi * sin(pi*ls(i,j,k)/(alpha1*dx1)))
		
		if (t_tens == 1) then 
			hsx(i,j,k) = 0.5*ta1(i,j,k)*(1. + cos(pi*ls(i,j,k)/(alpha1*dx1)))/(alpha1*dx1)
			hsy(i,j,k) = 0.5*tb1(i,j,k)*(1. + cos(pi*ls(i,j,k)/(alpha1*dx1)))/(alpha1*dx1)
			hsz(i,j,k) = 0.5*tc1(i,j,k)*(1. + cos(pi*ls(i,j,k)/(alpha1*dx1)))/(alpha1*dx1)
		endif
		
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

    rho(i,j,k) = rho_f1 * (1.-hs(i,j,k)) + rho_f2 * hs(i,j,k)
    ls_mu(i,j,k) = mu_f1  * (1.-hs(i,j,k)) + mu_f2  * hs(i,j,k)

	enddo
	enddo
	enddo

END SUBROUTINE heaviside

!!!####################################################################################

SUBROUTINE mod_ls1(ls1,mod_ls,dlsdxa,dlsdya,dlsdza,dimx,dimy,dimz)
! cálculo do vetor normal e da curvatura
 	
	IMPLICIT NONE

	integer :: i, j, k,ihs,dimx,dimy,dimz
	real(8),save :: aux1, aux2
	real(8),dimension(dimx,dimy,dimz) :: ta1,tb1,tc1,td1,te1,tf1,dlsdxa,dlsdya,dlsdza,ls1,mod_ls

	!cálculo das derivadas
	CALL der_weno(ls1,ta1,tb1,tc1,td1,te1,tf1,ihs,dimx,dimy,dimz)

	!utilizar o maior valor em absoluto como derivada representativa (o WENO tem o positivo e o negativo)
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
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

	!cálculo do vetor normal dlsdxa, dlsdya e dlsdza
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx
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

	
END SUBROUTINE mod_ls1


!!!####################################################################################

SUBROUTINE curv_ls1(dlsdxa,dlsdya,dlsdza)
! cálculo do vetor normal e da curvatura
 
	USE ls_param
	
	IMPLICIT NONE

	integer :: i, j, k,ihs
	real(8),save :: aux1, aux2
	real(8),dimension(nx,ny,nz) :: ta1,tb1,tc1,td1,te1,tf1,dlsdxa,dlsdya,dlsdza,ddlsdx,ddlsdy,ddlsdz

	ihs = 1

	call wenox(dlsdxa,nx,ny,nz,dx,ta1,td1,ihs)
	call wenoy(dlsdya,nx,ny,nz,dy,tb1,te1,ihs)
	call wenoz(dlsdza,nx,ny,nz,dz,tc1,tf1,ihs)


	!cálculo da curvatura em si
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

END SUBROUTINE curv_ls1
