!Subrotina definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!!!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!!!Modificações
!Felipe Augusto R. Silva em 09/01/2022
!Bruna Fernanda Soares em 15/02/2023
!Karol Rocha Araújo em 20/02/2024


!Para sondas, nomear como: SNOME & unit=10000x...
!Para conservação de massa, nomear como: CONMASSA & unit=20000x...
!Para número de Courant, nomear como: NCOURANT & unit=30000x...
!Para contaegem temporal, nomear como: CONTEMP & unit=40000x...
!Para plotagens de planos (paraview), nomear como: PRTSCR & unit=50000x...
!Para o plot dos atributos, nomear como: PATRIB & unit=60000x...


SUBROUTINE PLOT()

USE velpre
USE cond
USE tempo
!$ USE omp_lib
	
IMPLICIT NONE

    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
    integer,dimension(8) :: values

    	
if (it == 0) then
	!Criar pastas se não existirem

	call system('mkdir -p arquivos')
	call system('mkdir -p dados')

	    ! using keyword arguments
	call date_and_time(date,time,zone,values)
	call date_and_time(DATE=date,ZONE=zone)
	call date_and_time(TIME=time)
	call date_and_time(VALUES=values)
	write(*,'(a,2x,a,2x,a)') " DATA-ANO/   MÊS/  DIA", "        TEMPO-HORA/MINUTO/SEGUNDO"
	write(*,'(A,3i5,A,3i5)') " DATA-",values(1:3), "		TEMPO-", values(5:7)
    
	!Modelo de turbulência
	if (m_turb == 0) then
		write(*,*) "Sem modelo de turbulência."		
	elseif (m_turb == 1) then		
		write(*,*) "Modelo de turbulência: LES Smagorinsky-Lilly Clássico."
	elseif (m_turb == 2) then		
		write(*,*) "Modelo de turbulência: RANS com energia cinemática turbulenta."		
	else	
		write(*,*) "Modelo de turbulência: DES."	
	endif

	!Termo Advectivo
	if (adv_type == 1) then
		write(*,*) "Termo advectivo: clássico."		
	elseif (adv_type == 2) then		
		write(*,*) "Termo advectivo: rotacional."
	elseif (adv_type == 3) then		
		write(*,*) "Termo advectivo: antissimétrico."			
	endif

	!Derivada Espacial para velocidade
	if (der == 1) then
		write(*,*) "Upwind."		
	elseif (der == 2) then		
		write(*,*) "Centrado."
	elseif (der == 3) then		
		write(*,*) "Upwind 2nd order."			
	endif


	!Método de pressão
	if (t_press == 0) then
		write(*,*) "Aproximação Hidrostática - não indicado para escoamentos complexos"		
	elseif (t_press == 1) then		
		write(*,*) "Aproximação Não-Hidrostática com Método da Projeção"	
	elseif (t_press == 2) then		
		write(*,*) "Aproximação Não-Hidrostática com IPA"		
	endif
	    
	!Terminal inicial
	write(*,*) " "
	write(*,*) "Parâmetros:"
	write(*,'(A,F10.4,A, A,F10.4,A, A,I10,A)') "  dt= ", dt, ";", "  dt_frame= ", dt_frame, ";", "  ts= ", ts, ";"
	write(*,'(A,I5,A, A,I5,A, A,I5,A)') "  nx= ", nx, ";", "  ny= ", ny, ";", "  nz= ", nz, ";"
	write(*,'(A,F10.4,A, A,F10.4,A, A,F10.4,A)') "  dx= ", dx, ";", "  dy= ", dy, ";", "  dz_repr= ", dza, ";"
	write(*,'(A,F10.4,A, A,F10.4,A, A,F10.4,A)') "  Lx= ", lx, ";", "  Ly= ", ly, ";", "  Lz= ", lz, ";"
	write(*,'(A,I2,A, A,I2,A, A,I2,A, A,I2,A, A,I2,A, A,I2,A)') "  ccx0= ", ccx0 , ";", "  ccxf= ", ccxf, ";", &
	"  ccy0= ", ccy0, ";", "  ccyf= ", ccyf, ";", "  ccz0= ", ccz0, ";", "  cczf= ", cczf, ";"
	!write(*,'(A,I2,A, A,I2,A)') "  der= ", der, ";", "  advectivo= ", adv_type, ";" !, "modelo de turb. = ", m_turb

	else

	!Cálculo para o a estimativa do tempo restante
	ciclo = (agora(5)-agora1(5)) * 60 * 60 + (agora(6)-agora1(6)) * 60 + (agora(7)-agora1(7)) + real(agora(8)-agora1(8))/1000
	prev = (prev*6 + (ts-it)*ciclo*1./(60.*60.))/(7.)
	agora1 = agora
	call date_and_time(values = agora)
	call cpu_time(t_a)
	!$ t_a = omp_get_wtime()

	endif

! a pressão é normalizada por um ponto

if (t_sonda == 1) then
	
	!Sondas
	CALL SPLONGDESH1() !Plot perfil longitudinal do desnível level-set-WENO-tests Wave
	CALL SPLONGDESH2() !Plot perfil longitudinal do desnível
	CALL SPLONGDESH3() !Plot perfil longitudinal do desnível level-set-WENO-tests DamBreak
	CALL SUMED()       !Velocidade média (umed) para Reynolds
	!*******************************************************************************
	!Simulação canal com rugosidade uniforme experimento GS_H128 (Zampiron et al., 2022)
	CALL SXH20() !Sonda x/H=20
	CALL SXH31() !Sonda x/H=31
	CALL SXH51() !Sonda x/H=51
	!*******************************************************************************

endif

!Conservação de massa
CALL CONMASSA()


!CALL NCOURANT() !! não está criando dados para o arquivo. Arrumar se necessário


!Contagem temporal
CALL CONTEMP()


!Plotagens de planos (paraview)
CALL PRTSCR()


if (it == ts) CALL PATRIB()

ENDSUBROUTINE PLOT

!#######################################################################################


SUBROUTINE SXH20() !Simulação canal com rugosidade uniforme experimento GS_H128 (Zampiron et al., 2022)
!Sonda x/H=20

USE ls_param
USE velpre

IMPLICIT NONE

!Declarado também no programa
real(8),dimension(nz) :: celula
integer :: i, j, k, nza
	
if (it == 0) open (unit=100001, action= 'write', file= 'dados//sondaxH20', status= 'unknown')

i = int(2.54/dx)
j = int(0.5/dy)   
do k = 1, nz
	if (ls(i,j,k)>=0.) then
		celula(k)=k
		nza=k
	endif
enddo
		
if (it == 0) write(100001,*) 't', celula(1:nza) 
write(100001,*) t, u(i,j,1:nza)

if (it == ts) 	close (unit=100001)

ENDSUBROUTINE SXH20

!#######################################################################################


SUBROUTINE SXH31() !Simulação canal com rugosidade uniforme experimento GS_H128 (Zampiron et al., 2022)
!Sonda x/H=31

USE ls_param
USE velpre

IMPLICIT NONE

!Declarado também no programa
real(8),dimension(nz) :: celula
integer :: i, j, k, nza
	
if (it == 0) open (unit=100002, action= 'write', file= 'dados//sondaxH31', status= 'unknown')

i = int(3.937/dx)
j = int(0.5/dy)   
do k = 1, nz 
	if (ls(i,j,k)>=0.) then
 		celula(k)=k
 		nza=k
 	endif
enddo
		
if (it == 0) write(100002,*) 't', celula(1:nza)
write(100002,*) t, u(i,j,1:nza)

if (it == ts) 	close (unit=100002)
	 
ENDSUBROUTINE SXH31

!#######################################################################################


SUBROUTINE SXH51() !Simulação canal com rugosidade uniforme experimento GS_H128 (Zampiron et al., 2022)
!Sonda x/H=51

USE ls_param
USE velpre

IMPLICIT NONE

!Declarado também no programa
real(8),dimension(nz) :: celula
integer :: i, j, k, nza
	
if (it == 0) open (unit=100003, action= 'write', file= 'dados//sondaxH51', status= 'unknown')
		
i = int(6.477/dx)
j = int(0.5/dy)  

do k = 1, nz
	if (ls(i,j,k)>=0.) then
 		celula(k)=k
 		nza=k
	endif
enddo
		
if (it == 0) write(100003,*) 't', celula(1:nza)
write(100003,*) t, u(i,j,1:nza)

if (it == ts) 	close (unit=100003)
	 
ENDSUBROUTINE SXH51

!#######################################################################################


SUBROUTINE SPLONGDESH1() !Plot level-set-WENO-tests Wave
!Perfil longitudinal do desnível 

USE ls_param

IMPLICIT NONE

!Declarado também no programa
integer :: i, j, k, ii
	
if (it == 0) open (unit=100004, action= 'write', file= 'dados//h1.txt', status= 'unknown') 

ii = 0
i = 1 !int(1.0/dx) 
j = 1 ! int(0.5/dy) 
do k = nz, 2, -1
	if (ii == 0) then
		if (ls(i,j,k-1)*ls(i,j,k) <= 0.) then
			write(100004,*) dt*it, z*(k) + ls(i,j,k)
			ii = 1 !Fará com que saida de todos os loops
		endif
	endif
enddo

if (it == ts) 	close (unit=100004)

ENDSUBROUTINE SPLONGDESH1

!#######################################################################################


SUBROUTINE SPLONGDESH2() 
!Perfil longitudinal do desnível 

USE ls_param

IMPLICIT NONE

!Declarado também no programa
integer :: i, j, k, ii
	
if (it == 0) open (unit=100005, action= 'write', file= 'dados//h2.txt', status= 'unknown') 

ii = 0
i = int(1.0/dx) 
j = int(0.5/dy) 
do k = nz, 2, -1
	if (ii == 0) then
		if (ls(i,j,k-1)*ls(i,j,k) <= 0.) then
			write(100005,*) dt*it, z*(k) + ls(i,j,k)
			ii = 1 !Fará com que saida de todos os loops
		endif
	endif
enddo

if (it == ts) 	close (unit=100005)

ENDSUBROUTINE SPLONGDESH2

!#######################################################################################


SUBROUTINE SPLONGDESH3() !Plot do DamBreak
!Perfil longitudinal do desnível 

USE ls_param

IMPLICIT NONE

!Declarado também no programa
integer :: i, j, k, ii
	
if (it == 0) open (unit=100006, action= 'write', file= 'dados//h3.txt', status= 'unknown') 

ii = 0
j =  int(real(ny)/2.) 
do i = nx, 2, -1 !Tem que ser o último a variar
do k = 1, nz
	if (ii == 0) then
		if (ls(i-1,j,k)*ls(i,j,k) <= 0.) then
			write(100006,*) dt*it, dx*(i-0.5) + ls(i,j,k)
			ii = 1 !Fará com que saida de todos os loops
		endif
	endif
enddo
enddo

if (it == ts) 	close (unit=100006)

ENDSUBROUTINE SPLONGDESH3

!#######################################################################################


SUBROUTINE SUMED() 
!Velocidade média (umed) para Reynolds

USE ls_param
USE velpre

IMPLICIT NONE

!Declarado também no programa
real(8) :: umed
integer :: i, j, k, n
	
if (it == 0) open (unit=100007, action= 'write', file= 'dados//umed.txt', status= 'unknown') 

umed = 0
n = 0
do k = 1, nz
do j = 1, ny
do i = 30, nx1-15		
	if (ls(i,j,k) >= 0.) then
		umed = umed + u(i,j,k) 
		n = n+1		
	endif
enddo
enddo
enddo

umed = umed / n 
	if (it == 0) write(100007,*) 't', "u_med"
	write(100007,*) t, umed

if (it == ts) 	close (unit=100007)

ENDSUBROUTINE SUMED

!#######################################################################################


SUBROUTINE CONMASSA()
!Conservação de massa [20000x...]

USE velpre
USE ls_param
USE tempo

IMPLICIT NONE
integer :: i,j,k
real(8) a, b, c


real(8),dimension(nx,ny,nz1) :: wdz
	
if (it == 0) then
	open (unit=200001, action= 'write', file= 'dados//conservacao_massa.txt', status= 'unknown')
	write(200001,*) "t", " ", "vol_ini", " ", "vol_ins", " ","vol_ini-vol_ins", " ","divergencia"
endif

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		div(i,j,k) = (u(i+1,j,k) - u(i,j,k))/dx +   (v(i,j+1,k) - v(i,j,k))/dy +  (w(i,j,k+1) - w(i,j,k))/dzz(k)
	enddo
	enddo
	enddo

	write(200001,*) t, vol_ini, vol_ins, vol_ini-vol_ins, maxval(abs(div))
	
	! parte da rotina que reavalia o dt caso o usuário escolha "tempo variável"
	if (t_tempo_var .eq. 1) then
	
	
	
		a = maxval(abs(u(1:nx1,1:ny,1:nz)))*dt/dx 
		b = maxval(abs(v(1:nx,1:ny1,1:nz)))*dt/dy 
		
		
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			wdz(i,j,k) = w(i,j,k)/dzz(k)
		enddo
		enddo
		enddo	
		c = maxval(abs(w(1:nx,1:ny,1:nz1)))*dt





		if (a > 0.3 .or. b > 0.3 .or. c > 0.3) then
			dt = dt - 0.1*dt0
			write(*,*) "t", "   ", "maxval(abs(div))", "   ", "a", "   ","b", "   ","c" ,  "   ", "dt"
			write(*,*) t, maxval(abs(div)), a,b,c   , dt
			call coef_tempo()
		elseif (a < 0.1 .and. b < 0.1 .and. c < 0.1) then
			dt = dt + 0.01*dt0
			write(*,*) "t", "   ", "maxval(abs(div))", "   ", "a", "   ","b", "   ","c" ,  "   ", "dt"
			write(*,*) t, maxval(abs(div)), a,b,c   , dt
			call coef_tempo()
		endif
	endif

	if (t > (cont-9)*dt_frame) then
		write(*,*) " "
	 	write(*,'(A,F10.4,A, A,F10.4,A, A,F10.4,A, A,F10.4,A)') "  *vol. inicial:", vol_ini, ";", "  vol. instant.:", vol_ins, ";", &
	 	"  erro (%):", (vol_ini-vol_ins)/vol_ini, ";", "  div.:", maxval(abs(div)), ";"
		write(*,*) 
	endif
		
if (it == ts) 	close (unit=200001)

ENDSUBROUTINE CONMASSA

!#######################################################################################


SUBROUTINE NCOURANT()
!Número de Courant [30000x...]

USE cond
USE disc

IMPLICIT NONE

open (unit=300001, action= 'write', file= 'dados//courant.txt', status= 'unknown')
	write(300001,*) "t", "  ", "ntal", "  ", "tal", "  ", "maxval(a)", " ", "maxval(d)", "  ", "loca", "  ", "locd"
	
if (it == ts) 	close (unit=300001)

ENDSUBROUTINE NCOURANT

!#######################################################################################


SUBROUTINE CONTEMP()
!Contaegem temporal [40000x...]

USE tempo
USE disc
!$ USE omp_lib

IMPLICIT NONE

if (it == 0) then
	open (unit=400001, file='dados//contagem')
	write(400001,*) "it", " ", "t", " ", "hoje", " ", "agora", " ", "duração da simulação (min)"

	call date_and_time(values = agora)
	agora1 = agora
	call cpu_time(t_i)
	!$ t_i = omp_get_wtime()
	
	write(400001,*) it, t, agora, t_i - t_i
else

        ciclo = (agora(5)-agora1(5)) * 60 * 60 + (agora(6)-agora1(6)) * 60 + (agora(7)-agora1(7)) + real(agora(8)-agora1(8))/1000
	prev = (prev*6 + (ts-it)*ciclo*1./(60.*60.))/(7.)
	agora1 = agora
	call date_and_time(values = agora)
	call cpu_time(t_a)
	!$ t_a = omp_get_wtime()

	
	write(400001,*) it, t, ciclo, prev, (t_a-t_i)/60.
	
	if (t > (cont-9)*dt_frame) then
	
		!Contagem temporal
		!write(*,*) "it,", " ", "t,", " ", "ciclo,", " ", "tempo restante aproximado (horas),", " ", "duração da simulação (min)"
		!write(*,'(I6, F10.4, F10.4, F10.4, F10.4)') it, t, ciclo, prev, (t_a-t_i)/60.
		
		write(*,'(A,I6,A, A,F10.4,A, A,F10.4,A, A,F10.4,A, A,F10.4,A)') "  it:", it, ";", "  t:", t, ";", &
		"  ciclo:", ciclo, ";", "  tempo restante aproximado (horas):", prev, ";", "  duração da simulação (min):", (t_a-t_i)/60., ";"
		
	
	
	endif
	
endif


if (it == ts) 	close (unit=400001)

ENDSUBROUTINE CONTEMP

!#######################################################################################


SUBROUTINE PRTSCR()
!Plotagens de planos (escolher um por simulação, fazer alterações dependendo do plano a ser extraído) [50000x...]

	USE ls_param
	USE velpre
	USE tempo
	USE obst
	USE mms_m
	USE cond
	USE les

	IMPLICIT NONE

	!Declarado também no programa
	real(8),dimension(nx1,ny1,nz1) :: uaux, vaux, waux, x1, y1, z1
	real(8),dimension(nx,ny,nz)    :: dudy, dudz, dvdx, dvdz, dwdx, dwdy, kaaux
	real(8),dimension(nx,ny,nz)    :: prdaux, kaux, vorti, vortj, vortk, lasux
	real(8),dimension(nx1,ny,nz1)  :: auxy
	real(8),dimension(nx1,ny1,nz)  :: auxz

	
	integer :: ifile, nfil, i, j, k, ii,nza
	
	!Número do arquivo de saída
	integer :: dig1, dig2, dig3, dig4, dig5

	!Nome do arquivo de saída
	character(5) chits

	open (unit=5468, file='dados//ku',status='unknown')
	do j = 1, ny
		write(5468,*) ku(:,j)
	enddo
	close(5468)
	
	open (unit=5468, file='dados//kv',status='unknown')
	do j = 1, ny1
		write(5468,*) kv(:,j)
	enddo
	close(5468)

	open (unit=5468, file='dados//kw',status='unknown')
	do j = 1, ny
		write(5468,*) kw(:,j)
	enddo
	close(5468)	
	
	
	
	if (t > (cont-9)*dt_frame .or. it == 0) then
	

		if (it == 0) then
			cont = 10
		else
			cont = cont + 1
			
			do ifile = 1, cont
				dig1 =    ifile/10000 + 48
				dig2 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
				dig3 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
				dig4 = ( ifile - 100*( ifile/100 ) )/10 + 48
				dig5 = ( ifile - 10*( ifile/10 ) )/1 + 48
				chits(1:5) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//char(dig5)
			enddo	
		endif
		

		do k = 1, nz1
		do j = 1, ny1
		do i = 1, nx1
		uaux(i,j,k) = (u(i,j,k) + u(i,j-1,k) + u(i,j,k-1) + u(i,j-1,k-1)) * 0.25
		vaux(i,j,k) = (v(i,j,k) + v(i-1,j,k) + v(i,j,k-1) + v(i-1,j,k-1)) * 0.25
		waux(i,j,k) = (w(i,j,k) + w(i,j-1,k) + w(i-1,j,k) + w(i-1,j-1,k)) * 0.25

		enddo
		enddo
		enddo

		!Posições nos nós
		do k = 1, nz1
		do j = 1, ny1
		do i = 1, nx1
			x1(i,j,k) = xx(i*2-1)
			y1(i,j,k) = yy(j*2-1)
			z1(i,j,k) = zz(k*2-1)
		enddo
		enddo
		enddo


		if (obst_t == 0) then
			obs_ls = 0
		endif

		prdaux(1:nx,1:ny,1:nz) = prd0(1:nx,1:ny,1:nz) + prd1(1:nx,1:ny,1:nz) - (prd0(nx,1,nz) + prd1(nx,1,nz))
		kaaux(1:nx,1:ny,1:nz) = ka(1:nx,1:ny,1:nz)


		if (mms_t > 0) then
			CALL interpy_cf(sqrt(erro_u),nx1,ny,nz,auxz)
			CALL interpz_cf(auxz,nx1,ny1,nz,dz,uaux)

			CALL interpx_cf(sqrt(erro_v),nx,ny1,nz,auxz)
			CALL interpz_cf(auxz,nx1,ny1,nz,dz,vaux)

			CALL interpx_cf(sqrt(erro_w),nx,ny,nz1,auxy)
			CALL interpy_cf(auxy,nx1,ny,nz1,waux)

			prdaux = sqrt(erro_p)
		endif	

		if (t_plot == 0) then
		
			if (it == 0) then
			 open (unit=cont, file='arquivos//campos_00010',form='unformatted',status='unknown')
			else
			 open (unit=cont, action= 'write', file= 'arquivos//campos_'//chits,form='unformatted',status='unknown')
			endif
						
			write(cont) real(x1,4),real(y1,4),real(z1,4),real(uaux,4), & 
			real(vaux,4),real(waux,4),real(ls,4),real(obs_ls,4)
			
		elseif (t_plot == 1) then
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
			dudy(i,j,k) = ((u(i,j+1,k)-u(i,j,k))/dy)
			dudz(i,j,k) = ((u(i,j,k+1)-u(i,j,k))/dz(k)) 
			dvdx(i,j,k) = ((v(i+1,j,k)-v(i,j,k))/dx)
			dvdz(i,j,k) = ((v(i,j,k+1)-v(i,j,k))/dz(k))
			dwdx(i,j,k) = ((w(i+1,j,k)-w(i,j,k))/dx)
			dwdy(i,j,k) = ((w(i,j+1,k)-w(i,j,k))/dy) 

			vorti(i,j,k) = dwdy(i,j,k) - dvdz(i,j,k)
			vortj(i,j,k) = dudz(i,j,k) - dwdx(i,j,k)
			vortk(i,j,k) = dvdx(i,j,k) - dudy(i,j,k)
			enddo
			enddo
			enddo

			if (it == 0) then
			 open (unit=cont, file='arquivos//campos_00010',form='unformatted',status='unknown')
			else
			 open (unit=cont, action= 'write', file= 'arquivos//campos_'//chits,form='unformatted',status='unknown')
			endif
			
			write(cont) real(x1,4),real(y1,4),real(z1,4),real(uaux,4), & 
			real(vaux,4),real(waux,4),real(ls,4),real(obs_ls,4),real(prdaux,4), &
			real(vorti,4),real(vortj,4),real(vortk,4),real(nut,4),real(kaaux,4)
		endif
		close (unit=cont)

		!Chama paraview
		if (it == 0) CALL visu ()
	endif

ENDSUBROUTINE PRTSCR

!#######################################################################################


SUBROUTINE PATRIB()

	USE velpre
	USE obst

	IMPLICIT NONE

	open(unit=600001, action= 'write', file= 'dados//atributos.txt', status= 'unknown')
	write(600001,*) "u inicial = ", uinicial
	write(600001,*) "dx = ", dx, ", dy = ", dy, ", dz_repr = ", dza
	write(600001,*) "nx = ", nx, "ny = ", ny, "nz = ", nz
	write(600001,*) "dt = ", dt, "ts = ", ts, "duração da sim =", (t_a-t_i)/60.
	close (unit=600001)

ENDSUBROUTINE PATRIB


