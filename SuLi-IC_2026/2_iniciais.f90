!Subrotina: definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!!!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!!!Modificações
!Leonardo Romero Monteiro em 13/01/2022
!Bruna Fernanda Soares em 23/02/2024

SUBROUTINE parametros()

	USE disc
	USE restart
	USE cond
	USE param
	USE les
	USE ls_param

implicit none
	
character :: a*80

1000 format(a,80x) 
open(1,file='00_suli3d.prm',status='unknown',form='formatted') 
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) dx
read (1,*) dy
read (1,*) dzmin
read (1,*) dzmax
read (1,*) dz_t
read (1,*) nx
read (1,*) ny
read (1,*) nz
read (1,*) dt0 
read (1,*) dt_frame
read (1,*) t_s 
read (1,*) uinicial
read (1,*) iturb
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) t_plot 
read (1,*) t_sonda
read (1,*) t_tempo 
read (1,*) t_tempo_var 
read (1,*) der 
read (1,*) adv_type 
read (1,*) ibm_t 
read (1,*) obst_t 
read (1,*) m_turb 
read (1,*) ct_var
read (1,*) esp_type
read (1,*) wave_t 
read (1,*) t_press
read (1,*) pcg_prec_t
read (1,*) mms_t 
read (1,*) t_tens
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) irest
read (1,*) interv_rest
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) ccx0
read (1,*) ccxf
read (1,*) ccy0
read (1,*) ccyf
read (1,*) ccz0
read (1,*) cczf
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) chezy
read (1,*) decliv
read (1,*) z0
read (1,*) cka
read (1,*) csmag
read (1,*) cmu
read (1,*) cdes
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) alpha1
read (1,*) rho_f1
read (1,*) mu_f1
read (1,*) rho_f2
read (1,*) mu_f2
read (1,*) sigma
read (1,*) tipo 
read (1,*) ampl
read (1,*) lambdax
read (1,*) lambday
read (1,*) prof
read (1,*) m
read (1,1000) a 
read (1,*) numb_threads
read (1,*) omp_t

 close(1) 

nx1=nx+1
ny1=ny+1
nz1=nz+1
nxx = nx1*2
nyy = ny1*2
nzz = nz1*2

CALL init_disc()
CALL dz_dist()

dt = dt0
lx = dx*nx
ly = dy*ny
lz = sum(dz(1:nz)) !configurado como valor de entrada na verdade temos mais valores de entrada do que precisa

ts = ceiling(t_s/dt0)+1
gx = 9.80665 * dsin(datan(decliv)); gz = 9.80665 * dcos(datan(decliv))

delta = max(dx,dy,dza)
deltai(0:nz+1) = (dx*dy*dz(0:nz+1))**(1./3.)
 
CALL init_variables1()
CALL init_variables2()
CALL init_variables3()
CALL init_variables4()
CALL init_variables5()
CALL init_variables6()
CALL init_variables7()

END SUBROUTINE parametros

!##################################################################################################################

SUBROUTINE iniciais()


	USE velpre
	USE obst
	USE tempo_hora_data
	USE cond
	USE mms_m
	USE les
	USE ls_param
	
	IMPLICIT NONE

	integer :: i, j, k, ik,iik
	real(8) :: umed, aux1
	real(8),dimension(nx1,ny,nz) :: ls_x
	real(8),dimension(nx,ny,nz)  :: uc
	real(8),dimension(3):: rr
	
	write(*,'(A,I10)') "Tamanho do domínio:", nx*ny*nz
		
	if (nx*ny*nz > 30000000) then
		write(*,*) "Verifique se o seu computador tem capacidade para a simulação, se sim, cancele esta condicional no código."
		STOP
	endif
	
	if (ccx0.eq.4) then 
		!Para validacao
		!u(:,:,11+1) = 0.15946
		!u(:,:,12+1) = 0.2873
		!u(:,:,13+1) = 0.328
		!u(:,:,14+1) = 0.36376
		!u(:,:,15+1) = 0.39226
		!u(:,:,16+1) = 0.41742
		!u(:,:,17+1) = 0.44166
		!u(:,:,18+1) = 0.46318
		!u(:,:,19+1) = 0.48141
		!u(:,:,20+1) = 0.4867
	endif

	t = 0.
	dt = dt0
	
	contp = 0
	v = 0.
	w = 0.
	prd0 = 0.
	prd1 = 0.
	px = 0.
	py = 0.
	pz = 0.
	ub = 0.
	vb = 0.
	wb = 0.
	d_max = 0.
	d_min = 0.

    	hsx = 0
    	hsy = 0
    	hsz = 0

	CALL level_set_ini()
	CALL interpx_cf(ls,nx,ny,nz,ls_x)
		
	if (ccx0.eq.0 .or. ccx0.eq.3) then
		!Velocidade de entrada com incremento de turbulência 
			
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			if (ls_x(i,j,k)>=-2*dz(k)) then
				call random_number(rr(1))
				rr(1) = 2.*(rr(1)-0.5)
				u(i,j,k) = uinicial*(1.+rr(1)*iturb)				
				umed = u(i,j,k) + umed
				ik   = 1 + ik
				iik  = k		
			else
				u(i,j,k) = u(i,j,iik)*(1.+dtanh(pi*(iik-k-1)/2.))
			endif
		enddo
		enddo
		enddo
		
	else
		u = uinicial
	endif

	!Utilizado para o modelo de turbulência Smag
	if (m_turb == 1) then
	
		aux1 = iturb*iturb*(uinicial*uinicial+10.**(-6))*1.5 ! ka = turbulence kinetic energy ! ruído (10-6) adicionado para casos sem velocidade inicial.
		

		do k = 1, nz
		do j = 1, ny
		do i = 1, nx

			if (ls(i,j,k)>=-2*dz(k)) then
				
			call random_number(rr)

			
			u(i,j,k) = u(i,j,k)+ aux1**0.5*rr(1)
			v(i,j,k) = v(i,j,k)+ aux1**0.5*rr(2)
			w(i,j,k) = w(i,j,k)+ aux1**0.5*rr(3)		
	
			endif	
			
		enddo
		enddo
		enddo
		

	
	!Utilizado para o modelo de turbulência DES
	elseif (m_turb .ge. 2) then
	
		aux1 = iturb*iturb*(uinicial*uinicial+10.**(-6))*1.5 ! ka = turbulence kinetic energy ! ruído (10-6) adicionado para casos sem velocidade inicial.
		

		
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
		
			!if (ls(i,j,k)>=-2*dz(k)) then
				
			ka(i,j,k) = aux1
			
			!call random_number(rr)

			
			!u(i,j,k) = u(i,j,k)+ aux1**0.5*rr(1)
			!v(i,j,k) = v(i,j,k)+ aux1**0.5*rr(2)
			!w(i,j,k) = w(i,j,k)+ aux1**0.5*rr(3)	
			
			!endif	
			
		enddo
		enddo
		enddo
		

		
		call contorno_les()
	endif

	!Utilizado para definir obstáculos de fundo
	CALL obstaculo()

	!Utilizado para definir onda
	if (wave_t > 0) then
		write(*,*) "Com entrada de onda."
		CALL waves_coef()
		call boundary_waves()
	else
		write(*,*) "Sem entrada de onda."
		bxx1(:,:) = uinicial
		bxx0(:,:) = uinicial
	
	! condição antiga que não lembro mais pra que servia.	
	!	if (obst_t .ne. 0) then
	!		do j = 0, ny+1
	!		bxx1(:,0:ku(1,j)) = 0.
	!		bxx0(:,0:ku(0,j)) = 0.
	!		enddo
	!	endif	
		bxy0 = 0.
		bxz0 = 0.
	endif

	if (mms_t == 0) then
	   tf_p = 0.
	   tf_u = 0.
	   tf_v = 0.
	   tf_w = 0.
	endif

	!inicialização do código
	it = 0
	cont = 10

	call coef_tempo()

END SUBROUTINE iniciais

!##################################################################################################################

SUBROUTINE coef_tempo()

	USE disc

	IMPLICIT NONE
	
	!Método da integração temporal para velocidades
	if (t_tempo == 0) then
	if (it < 1) write(*,*) "Esquema temporal: Euler."
	if (it < 1) call print_pcg_prec()
		ntt = 1
		a_dt = 1.0*dt
		
	elseif (t_tempo == 1) then
	if (it < 1) write(*,*) "Esquema temporal: RK2."
	if (it < 1) call print_pcg_prec()
		ntt = 2
		a_dt(1) = 0.5*dt
		a_dt(2) = 1.0*dt
		a_dt(3) = 1.0*dt
		
	elseif (t_tempo == 2) then
	if (it < 1) write(*,*) "Esquema temporal: RK3."
	if (it < 1) call print_pcg_prec()
		ntt = 3
		a_dt(1) = 0.5*dt
		a_dt(2) = 1.0*dt
		a_dt(3) = 1.0*dt

	elseif (t_tempo == 3) then
	if (it < 1) write(*,*) "Esquema temporal: AB2."
	if (it < 1) call print_pcg_prec()
		ntt = 1
		a_dt = 1.0*dt

    	elseif (t_tempo == 4) then
	if (it < 1) write(*,*) "Esquema temporal: AB3."
	if (it < 1) call print_pcg_prec()
		ntt = 1
		a_dt = 1.0*dt
	endif

CONTAINS

	SUBROUTINE print_pcg_prec()
		USE disc, only: pcg_prec_t, t_press
		IMPLICIT NONE

		if (t_press == 2) then
			write(*,*) "Poisson (t_press=2): FFTW 3D (sem CG/PCG). pcg_prec_t ignorado."
			return
		endif

		select case (pcg_prec_t)
		case (0,1)
			write(*,'(A,I2,A)') "PCG precondicionador: nenhum (pcg_prec_t=", pcg_prec_t, ")."
		case (2)
			write(*,*) "PCG precondicionador: line-z (tridiagonal) (pcg_prec_t=2)."
		case (3)
			write(*,*) "PCG precondicionador: SGS red-black (pcg_prec_t=3)."
		case default
			write(*,'(A,I6,A)') "PCG precondicionador: valor invalido (pcg_prec_t=", pcg_prec_t, ") -> default=2."
		end select
	END SUBROUTINE print_pcg_prec

END SUBROUTINE coef_tempo

!##################################################################################################################
! calcula como será a variação do dz ao longo de y
SUBROUTINE dz_dist()

	USE disc
	USE ls_param
	
	IMPLICIT NONE
	
	integer :: i, j, k
	real(8) :: dzsum, dist, coef
	real(8),dimension(0:nz+1)  :: hss

dza = dzmin !aqui estou considerando a referência como o dz mínimo	

if (dz_t == 0 ) then ! sem variação
	do k = 0, nz+1
		dz(k) = dza!(dz(k)+dz(k-1))*0.5 
	enddo

elseif (dz_t == 1) then!variação linear
	do k = 0, nz+1
		dz(k) = dzmin + (dzmax-dzmin)*k/nz
	enddo


elseif (dz_t == 2) then!variação logarítmica
	do k = 1, nz+1
		dz(k) = dzmin + log(real(k))*dzmax
		if (dz(k) > 0.005) dz(k) = 0.005 
	enddo
	dz(0) = dz(1) 
	
		
	
elseif (dz_t == 3) then!variação logarítmica


	dist = 2
	coef = 4 !10
	do k = 0, nz+1
	
		hss(k) = 0.5 * (1.+ dtanh((-dist + real(k)/coef)))
		dz(k) = dzmin * (1.-hss(k)) + dzmax * hss(k)
		
	enddo




elseif (dz_t == 4) then !variação a partir da superfície livre

	dzsum = 0.
	do k = 1, nz
		if (prof > dzsum)  then
			dz(k) = dzmin
			dzsum = dzsum + dz(k)
		else
			dz(k) = dz(k-1) + 0.001*dzmin 
			if (dz(k) > dzmax) dz(k) = dzmax
			dzsum = dzsum + dz(k)
		endif
	enddo
	dz(0) = dz(1) 
	dz(nz+1) = dz(nz) 	

endif


!variação a partir por Heaviside


!! cálculo do dz referente a w
do k = 1, nz+1
	dzz(k) = (dz(k)+dz(k-1))*0.5 
enddo
dzz(0) = dzz(1)
dzz(nz+1) = dzz(nz)
dzz(nz1+1) = dzz(nz+1)

!cálculo de distâncias

! estas vão até o centro da célula
do i = 1,nx
	x(i)  = (i-0.5)*dx
enddo

do j = 1,ny
	y(j)  = (j-0.5)*dy
enddo

z(1) = dz(1)*0.5

do k = 2, nz
	z(k) = z(k-1) + (dz(k-1)+dz(k))*0.5
enddo

! estas passam pelas as arestas e centro das células
do i = -1, nxx+2
	xx(i)  = (i-1.)*dx*0.5
enddo
		
do j = -1, nyy+2
	yy(j)  = (j-1.)*dy*0.5
enddo		


zz(-1) = -dz(1)      !! condição de contorno
zz(0)  = -dz(1)*0.5  !! condição de contorno
zz(1)  = 0.  
do k = 2, nzz+1
	zz(k) = zz(k-1) + dz(floor(k*0.5))*0.5
enddo
zz(nzz+2) = zz(nzz+1) + dz(nz+1)*0.5


END SUBROUTINE dz_dist


