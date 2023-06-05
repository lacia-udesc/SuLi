!Subrotina: definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 13/01/2022

SUBROUTINE parametros()

	USE disc
	USE restart
	USE cond
	USE param
	USE smag
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
read (1,*) dz
read (1,*) nx
read (1,*) ny
read (1,*) nz
read (1,*) dt0 
read (1,*) dt_frame
read (1,*) t_s 
read (1,*) uinicial
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) t_plot 
read (1,*) t_tempo 
read (1,*) t_tempo_var 
read (1,*) der 
read (1,*) adv_type 
read (1,*) obst_t 
read (1,*) m_turb 
read (1,*) esp_type
read (1,*) wave_t 
read (1,*) t_press
read (1,*) mms_t 
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
read (1,*) csmag
read (1,1000) a 
read (1,1000) a 
read (1,1000) a 
read (1,*) alpha1
read (1,*) rho_f1
read (1,*) mi_f1
read (1,*) rho_f2
read (1,*) mi_f2
read (1,*) sigma
read (1,*) tipo 
read (1,*) ampl
read (1,*) lambdax
read (1,*) lambday
read (1,*) prof
read (1,*) m
read (1,1000) a 
read (1,*) numb_threads

close(1) 

dt = dt0
lx = dx*nx
ly = dy*ny
lz = dz*nz
nx1=nx+1
ny1=ny+1
nz1=nz+1
ts = ceiling(t_s/dt0)	
gx = 9.80665 * sin(atan(decliv)); gz = 9.80665 * cos(atan(decliv))

CALL init_variables1()
CALL init_variables2()
CALL init_variables3()
CALL init_variables4()
CALL init_variables5()
CALL init_variables6()
CALL init_variables7()

END SUBROUTINE parametros


SUBROUTINE iniciais()


	USE velpre
	USE obst
	USE tempo
	USE cond
	USE mms_m
	
	
	IMPLICIT NONE

	integer :: i, j, k
	
	write(*,*) "tamanho do domínio", nx*ny*nz
		
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
	u = uinicial
	v = 0.
	w = 0.
	prd0 = 0.
	prd1 = 0.
	ub = 0.
	vb = 0.
	wb = 0.
	d_max = 0.
	d_min = 0.

	!Utilizado para definir obstáculos de fundo
	CALL obstaculo()

	!Utilizado para definir onda
	if (wave_t > 0) then
		write(*,*) "Com entrada de onda."
		CALL waves_coef()
	else
		write(*,*) "Sem entrada de onda."
		bxx1(:,:) = uinicial
		bxx0(:,:) = uinicial
		
		if (obst_t .ne. 0) then
			do j = 0, ny+1
			bxx1(:,0:ku(1,j)) = 0.
			bxx0(:,0:ku(0,j)) = 0.
			enddo
		endif	
		bxy0 = 0.
		bxz0 = 0.
	endif

	CALL level_set_ini()

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE coef_tempo()

	USE disc
	!Método da integração temporal
	write(*,*) "Esquema temporal: "
	if (t_tempo == 0) then
		write(*,*) "Euler."
		ntt = 1
		a_dt = 1.0*dt
		
	elseif (t_tempo == 1) then
		write(*,*) "RK2."
		ntt = 2
		a_dt(1) = 0.5*dt
		a_dt(2) = 1.0*dt
		a_dt(3) = 1.0*dt
		
	elseif (t_tempo == 2) then
		write(*,*) "RK3."
		ntt = 3
		a_dt(1) = 0.5*dt
		a_dt(2) = 1.0*dt
		a_dt(3) = 1.0*dt
	elseif (t_tempo == 3) then
		write(*,*) "AB2."
		ntt = 1
		a_dt = 1.0*dt
    elseif (t_tempo == 4) then
		write(*,*) "AB3."
		ntt = 1
		a_dt = 1.0*dt
	endif

END SUBROUTINE coef_tempo

