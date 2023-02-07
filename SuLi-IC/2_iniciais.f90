!Subrotina: definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 13/01/2022

SUBROUTINE iniciais()

	USE velpre
	USE obst
	USE tempo
	USE cond
	USE mms_m
	USE ls_param, only: ls
	
	IMPLICIT NONE

	integer :: i, j, k, ik
	real(8),save :: umed
	real(8),save,dimension(nx1,ny,nz) :: ls_x
	
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
	
	CALL level_set_ini()

	call interpx_cf(ls,nx,ny,nz,ls_x) !(nx1,ny1,nz)
	
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		if (ls_x(i,j,k)>=0.) then
		
			call random_number(r)
				r=2.*(r-0.5)
			!#################
			!u = uinicial + t_iturb(r/100)
			u(i,j,k) = uinicial*(1.+r*iturb)
			!##############
			
			umed=u(i,j,k)+ umed
			ik = 1 + ik
			
		endif
	enddo
	enddo
	enddo
	
	umed=umed/ik
	
	u = u + (uinicial - umed)
	
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

