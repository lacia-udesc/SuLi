!Subrotina: definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 13/01/2015

SUBROUTINE iniciais()

USE velpre
USE obst
USE tempo
	USE cond
USE mms_m

IMPLICIT NONE

integer :: i, j, k

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

!Método da integração temporal
write(*,*) "Esquema temporal: "
if ((t_tempo == 0)  .or. (t_tempo == 3)) then
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
	
endif

	END SUBROUTINE iniciais

