!Subrotina definir as condições de contorno das velocidades e suas influências
! Referência: Gotoh, 2013

!!! Implementação 15/04/2014
! Leonardo Romero Monteiro

!!! Modificações
! Leonardo Romero Monteiro - 13/01/2015

SUBROUTINE iniciais()

	USE dzs
	USE velpre
	USE obst
	USE tempo
        USE cond
	USE mms_m

	IMPLICIT NONE
	!===================================================================================================================
	real(8) :: x, y
	integer :: i, j, k

if (nx*ny*nz > 30000000) then
write(*,*) "tente outro computador melhor!"
STOP
endif

!###! Criar pastas se não existirem !###!

call system('mkdir -p arquivos')
call system('mkdir -p dados')

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
	
	dz(:,:,0:nz-1) = dz1
	dz(:,:,nz:nz+1) = dztopo


	if (mms_t > 0) then
        	CALL mms_i()
	else
	
	do j = 1, ny
	y = (real(j)-0.5)*dy
	do i = 1, nx
	x = (real(i)-0.5)*dx
		eta0(i,j) = 0.1 *cos(x*pi/(nx*dx)) *cos(y*pi/(ny*dy))
		!eta0(i,j) = 0. !-0.1 *cos(0.314159265*(real(i)-1.)*dx)
	enddo
	enddo
	
	   tf_press = 0.
	   tf_u     = 0.
	   tf_v     = 0.
	   tf_w     = 0.
	endif
	
	
	

	eta0(0,:) = eta0(1,:)
	eta0(:,0) = eta0(:,1)
	eta0(nx+1,:) = eta0(nx,:)
	eta0(:,ny+1) = eta0(:,ny)
	eta1 = eta0

	do j = 0, ny+1
	do i = 0, nx+1
		dz(i,j,nz) = dz(i,j,nz) + eta0(i,j) !o topo do domínio segue a superfície livre	
	enddo
	enddo

	!! Utilizado para definir obstáculos de fundo
	CALL obstaculo()
	CALL dz_i()


	!! Utilizado para definir onda
	if (wave_t .ne. 0) then
		CALL waves_coef()
	endif



	!inicialização do código
	it = 0
	cont = 10

	END SUBROUTINE iniciais

