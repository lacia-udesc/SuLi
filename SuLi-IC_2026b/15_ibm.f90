!Subrotina definir IBM
!Leonardo Romero Monteiro


SUBROUTINE obstaculo()
!Subrotina para definir o formato do obstáculo

!Duna de forma senoidal (não varia em y, apenas repete)
!Daniel Rodrigues Acosta
!06/05/2016
!Tudo é resolvido em 2D e apenas repetido/transladado ao longo de Y, que para essa investigação terá espessura dy bem pequena (=3), para ignorar seus efeitos.
!Assim, tudo é definido sobre a projeção 2D e não é necessário utilizar índices j e k.

	USE velpre	
	USE obst
	USE disc
	USE cond
	
	!$ USE omp_lib

	IMPLICIT NONE

	real(8) :: x0, y0, a, b, c, d, sigx, sigy, sigz, tgaux, aux1, aux2, erro, raio, er 	!x mostra a localização atual
	real(8),dimension(-1:nx1*2+2,-1:ny1*2+2) :: auxx !!esse é o antigo
	real(8),dimension(-1:nxx+2,-1:nyy+2,-1:nzz+2) :: auxx_ls !esse é o novo
	real(8),dimension(nx,ny,nz) :: mod_ls !!esse é o antigo
	integer :: i,j,k,mm,n
	real(8),dimension(12,2) :: ponto

	!nxx = nx1*2
	!nyy = ny1*2
	!nzz = nz1*2

	!ku, kv, kw indicam até que altura as velocidades tem que ser zeradas (até qual índice k)

	!Desenhar um único objeto e fazer a divisão na hora de transformar em inteiro
	if (obst_t == 0) then
		write(*,*) "Sem obstáculo."
		! Inicializa os campos do IBM como "sem sólido" para que nenhuma
		! rotina posterior use valores não inicializados.
		ku = 0
		kv = 0
		kw = 0
		obs_ls  = -1.d0
		obs_lsx = -1.d0
		obs_lsy = -1.d0
		obs_lsz = -1.d0
		obs_lss = -1.d0
	return

	elseif (obst_t == 1) then !Rugosidade uniforme experimento GS_H128 (Zampiron et al., 2022)

        !$omp parallel do
		do k = -1,nzz+2
		do j = -1,nyy+2
		do i = -1,nxx+2
			auxx_ls(i,j,k)  = + 0.008*(1.+(dabs(dsin(pi*xx(i)/0.016))*dabs(dsin(pi*yy(j)/0.016)))**0.5)-zz(k)
			!auxx_ls(i,j,k)  = zz(5) -zz(k)
			
		enddo
		enddo
		enddo
        !$omp end parallel do

	elseif (obst_t == 2) then !(bd_koshizuka, 1995 e Kleefsman, 2005)

		!Condição inicial de barragem
		tgaux = 0.0!5715 ! dam depth
		sigx = 0.75!5715 ! dam x-length
		sigy = 0.5 ! dam y-length
		sigz = 0. ! dam y-length
		er = 0.0001 !para arrumar quando esta muito próximo do exato
		a = 0.08+er
		b = 0.2 +er
		c = 0.16+er
		mm = 50


        !$omp parallel do private(tgaux, aux1, aux2)
		do k = -1, nzz+2
		do j = -1, nyy+2
		do i = -1, nxx+2

			if (xx(i) <= 0.75) then
				tgaux = +xx(i) -0.75 + 0.08
			else
				tgaux = -xx(i) +0.75 + 0.08
			endif

			if (yy(j) <= 0.5) then
				aux1 = +yy(j) -0.5 + 0.20
			else
				aux1 = -yy(j) +0.5 + 0.20	
			endif

			aux2 = -zz(k) +0. + 0.16

			auxx_ls(i,j,k) = min(tgaux,aux1,aux2)

		enddo
		enddo
		enddo
        !$omp end parallel do

		! tentei criar o obstáculo com uma elipse, mas ela não cria uma função distância.
		!do k = -1, nzz+2
		!do j = -1, nyy+2
		!do i = -1, nxx+2
		!auxx_ls(i,j,k) = -((((xx(i)-sigx)/a)**mm + ((yy(j)-sigy)/b)**mm  + ((zz(k)-sigz)/c)**mm))**(1./mm)+1
		!enddo
		!enddo
		!enddo
		

				
	elseif (obst_t == 3) then !Degrau experimento T1 (Delft, 1980)


        !$omp parallel do
		do k = -1, nzz+2
		do j = -1, nyy+2
		do i = -1, nxx+2

			if (xx(i) <= (1.0+0.2)) then         ! +0.2m camada esponja
				auxx_ls(i,j,k) = 0.2 -zz(k)  !Raso plano

			else 
				auxx_ls(i,j,k) = -1. -zz(k)  !Fundo plano
			endif

		enddo
		enddo
		enddo
        !$omp end parallel do

				
	elseif (obst_t == 4) then !Corrego Grande
		er = 0.0
		sigy = 0. !correcao considerando o que foi anotado na folha
		sigz = 0. !correcao considerando o que foi anotado na folha
		n = 6
		ponto(1,1)  = 0.-sigy;     ponto(1,2)  = 0.5 + er -sigz
		ponto(2,1)  = 0.5-sigy;    ponto(2,2)  = 0.5 + er -sigz
		ponto(3,1)  = 0.7-sigy;   ponto(3,2)  = 0.0 + er -sigz
		ponto(4,1)  = 1.3-sigy;    ponto(4,2)  = 0.0 + er -sigz
		ponto(5,1)  = 1.5-sigy;   ponto(5,2)  = 0.50 + er -sigz
		ponto(6,1)  = 2.1-sigy;    ponto(6,2)  = 0.50 + er -sigz


	!corrego grande
	!ponto(1,1)  = 0.;      ponto(1,2)  = 0. + er
	!ponto(2,1)  = 0.001;   ponto(2,2)  = 0. + er
	!ponto(3,1)  = 0.1658;  ponto(3,2)  = 0. + er
	!ponto(4,1)  = 0.5;     ponto(4,2)  = 0.03 + er
	!ponto(5,1)  = 0.84;    ponto(5,2)  = 0.15 + er
	!ponto(6,1)  = 1.3;     ponto(6,2)  = 0.20 + er
	!ponto(7,1)  = 1.85;    ponto(7,2)  = 0.22 + er
	!ponto(8,1)  = 1.9;     ponto(8,2)  = 0.25 + er
	!ponto(9,1)  = 2.8;     ponto(9,2)  = 0.37 + er
	!ponto(10,1) = 4.8;     ponto(10,2) = 0.534 + er
	

	
		!$omp parallel do private(mm)
		do k = -1, nzz+2
		do j = -1, nyy+2
		do i = -1, nxx+2
			do mm = 2,n 
				if (yy(j) < ponto(mm,1) .and. yy(j) >= ponto(mm-1,1)) then    
					auxx_ls(i,j,k) = ponto(mm,2) + (ponto(mm,2)-ponto(mm-1,2))/(ponto(mm,1)-ponto(mm-1,1))*(yy(j)-ponto(mm,1))-zz(k) 
					EXIT ! intervalo encontrado, evita varrer os demais mm
				endif
			enddo
		enddo
		enddo
		enddo
		!$omp end parallel do
		
	endif

	! independente do IBM utilizado, este é o utilizado para plotar do paraview
	! Aqui definimos o IBM submalha e depois marcamos as duas células que estão dentro do obstáculo perto da interface.

    !$omp parallel do
    do k = 1, nz1
	do j = 1, ny1
	do i = 1, nx1
		!z = (k-1.0)*dz
		obs_ls(i,j,k) = auxx_ls(i*2-1,j*2-1,k*2-1) !- z
	enddo
	enddo
   	enddo
    !$omp end parallel do

	! Os índices ku/kv/kw são calculados mais abaixo, a partir da função
	! distância na malha deslocada (obs_lsx/obs_lsy/obs_lsz).

	if (ibm_t == 1 .or. ibm_t == 2 .or. ccz0 == 4) then !obs_ls* na malha deslocada (necessário ao IBM forçado e ao semi-slip)

		! em x   
        !$omp parallel
        !$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			obs_lsx(i,j,k) = auxx_ls(i*2-1,j*2,k*2)
		enddo
		enddo
	    enddo
	    !$omp end do nowait

		! em y
		!$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
			obs_lsy(i,j,k) = auxx_ls(i*2,j*2-1,k*2)
		enddo
		enddo
	    	enddo
		!$omp end do nowait

	  	! em z
		!$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			obs_lsz(i,j,k) = auxx_ls(i*2,j*2,k*2-1)
		enddo
		enddo
	    enddo
	    !$omp end do nowait

		! centro da célula
        !$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			obs_lss(i,j,k) = auxx_ls(i*2,j*2,k*2)
		enddo
		enddo
	    enddo
	    !$omp end do
	    !$omp end parallel
	endif

   	if (ibm_t == 1 .or. ibm_t == 2 .or. ccz0.eq.4) then

   	! ku/kv/kw: índice do nó MAIS ALTO de cada coluna que está dentro do
   	! sólido (obs_ls* > 0). Obtido diretamente da função distância na
   	! malha deslocada, em acordo com a forçagem direta. A varredura por
   	! coluna (seriada em k) evita as ambiguidades da antiga busca por
   	! mudança de sinal (que falhava para sólidos "flutuantes" e para
   	! colunas totalmente bloqueadas).
   	ku = 0
   	kv = 0
   	kw = 0

   	!$omp parallel
   	!$omp do
   	do j = 1, ny
	do i = 1, nx1
		do k = 1, nz
			if (obs_lsx(i,j,k) > 0.) ku(i,j) = k
		enddo
	enddo
	enddo
	!$omp end do nowait

	!$omp do
   	do j = 1, ny1
	do i = 1, nx
		do k = 1, nz
			if (obs_lsy(i,j,k) > 0.) kv(i,j) = k
		enddo
	enddo
	enddo
	!$omp end do nowait

	!$omp do
   	do j = 1, ny
	do i = 1, nx
		do k = 1, nz1
			if (obs_lsz(i,j,k) > 0.) kw(i,j) = k
		enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel

	! Extensão para as colunas fantasma (usadas nas cópias periódicas do
	! semi_slip). Em x segue o mapeamento periódico de contorno():
	!   u(0)=u(nx1-1), u(nx1+1)=u(2); v(0)=v(nx), v(nx+1)=v(1); w idem a v.
	!$omp parallel do
	do j = 1, ny
		if (ccx0 == 0 .and. ccxf == 0) then
			ku(0,j)     = ku(nx1-1,j)
			ku(nx1+1,j) = ku(2,j)
			kw(0,j)     = kw(nx,j)
			kw(nx+1,j)  = kw(1,j)
		else
			ku(0,j)     = ku(1,j)
			ku(nx1+1,j) = ku(nx1,j)
			kw(0,j)     = kw(1,j)
			kw(nx+1,j)  = kw(nx,j)
		endif
	enddo
	!$omp end parallel do

	!$omp parallel do
	do j = 1, ny1
		if (ccx0 == 0 .and. ccxf == 0) then
			kv(0,j)    = kv(nx,j)
			kv(nx+1,j) = kv(1,j)
		else
			kv(0,j)    = kv(1,j)
			kv(nx+1,j) = kv(nx,j)
		endif
	enddo
	!$omp end parallel do

	!$omp parallel do
	do i = 0, nx1+1
		ku(i,0)     = ku(i,1)
		ku(i,ny+1)  = ku(i,ny)
	enddo
	!$omp end parallel do

	!$omp parallel do
	do i = 0, nx+1
		kv(i,0)     = kv(i,1)
		kv(i,ny1+1) = kv(i,ny1)
		kw(i,0)     = kw(i,1)
		kw(i,ny+1)  = kw(i,ny)
	enddo
	!$omp end parallel do

	endif

	if (ibm_t == 2) then
	    	! para manter a função distância (muito provavelmente é desnecessário)
	    	!CALL reinic_weno(obs_ls,nx1,ny1,nz1)
	    	!CALL reinic_weno(obs_lsx,nx1,ny,nz)
	    	!CALL reinic_weno(obs_lsy,nx,ny1,nz)
	    	!CALL reinic_weno(obs_lsz,nx,ny,nz1)
	    	!CALL reinic_weno(obs_lss,nx,ny,nz)

	    	! cria marcadores que identificam as duas primeiras células para dentro do obstáculo
		CALL ibm_posicoes(id_ibmx,obs_lsx,nx1,ny,nz)
		CALL ibm_posicoes(id_ibmy,obs_lsy,nx,ny1,nz)
		CALL ibm_posicoes(id_ibmz,obs_lsz,nx,ny,nz1)
		CALL ibm_posicoes(id_ibm, obs_lss, nx,ny,nz)

		!CALL mod_ls1(obs_lss,mod_ls,norm_xx,norm_yy,norm_zz,nx,ny,nz)   !deletar depois
	endif

END SUBROUTINE obstaculo

!##################################################################################################################

!Subrotina de forçagem direta para o IBM do tipo 1
!Referência: Mohd-Yusof (1997) e Fadlun et al. (2000)
!
!A forçagem é aplicada como reconstrução da velocidade nos nós da malha
!deslocada que interagem com o sólido:
! 1) todos os nós internos ao sólido (obs_ls* > 0) são zerados;
! 2) o nó sólido adjacente à interface (nó fantasma) é reconstruído por
!    interpolação linear de modo que a velocidade interpolada na posição
!    exata da parede (obs_ls* = 0) seja nula (no-slip/impermeabilidade).
!    Quando ccz0 == 4 (semi-slip), o nó fantasma é posteriormente
!    substituído pela reconstrução com lei da parede em semi_slip().
! 3) A fração f = d_s/(d_s+d_f) é limitada a fmax = 0.9 (e d_f a eps).
!    Sem esse limite, quando a interface fica muito perto de um nó fluido
!    (d_f -> 0), o fantasma -f/(1-f)*u_f cresce sem limite (ou divide por
!    zero), contaminando o termo viscoso do 1º nó fluido e a reconstrução
!    de w em posdin(), o que desestabiliza a superfície livre após muitos
!    passos de tempo.
!
!Os loops são por coluna (i,j) e seriais em k, de modo que a zeragem dos
!nós sólidos e a reconstrução do fantasma não apresentam condição de
!corrida em OpenMP.
!
!Implementação em 08/2026

SUBROUTINE ibm_forcado()

	USE disc, only : nx, ny, nz, nx1, ny1, nz1
	USE obst, only : obs_lsx, obs_lsy, obs_lsz
	USE velpre, only : u, v, w
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k
	real(8) :: d_s, d_f, f, eps
	real(8), parameter :: fmax = 0.9d0 ! limita a amplificação do nó fantasma a f/(1-f) <= 9

	eps = 1.0d-12

	! componente u (nós da malha u: nx1 x ny x nz)
	!$omp parallel do private(k,d_s,d_f,f)
	do j = 1, ny
	do i = 1, nx1
		do k = 1, nz
			if (obs_lsx(i,j,k) > eps) then
				! nó dentro do sólido
				u(i,j,k) = 0.d0
				if (k > 1 .and. obs_lsx(i,j,k-1) <= -eps) then
					! sólido com fluido abaixo (parede inferior de objeto "flutuante")
					d_s = obs_lsx(i,j,k)
					d_f = max(-obs_lsx(i,j,k-1), eps)
					f = min(d_s/(d_s+d_f), fmax)
					u(i,j,k) = -f/(1.d0-f)*u(i,j,k-1)
				endif
			elseif (k > 1 .and. obs_lsx(i,j,k-1) > eps) then
				! nó fluido com sólido abaixo: reconstrução do nó fantasma
				! (no-slip na posição exata da interface)
				d_s = obs_lsx(i,j,k-1)
				d_f = max(-obs_lsx(i,j,k), eps)
				f = min(d_s/(d_s+d_f), fmax)
				u(i,j,k-1) = -f/(1.d0-f)*u(i,j,k)
			endif
		enddo
	enddo
	enddo
	!$omp end parallel do

	! componente v (nós da malha v: nx x ny1 x nz)
	!$omp parallel do private(k,d_s,d_f,f)
	do j = 1, ny1
	do i = 1, nx
		do k = 1, nz
			if (obs_lsy(i,j,k) > eps) then
				v(i,j,k) = 0.d0
				if (k > 1 .and. obs_lsy(i,j,k-1) <= -eps) then
					d_s = obs_lsy(i,j,k)
					d_f = max(-obs_lsy(i,j,k-1), eps)
					f = min(d_s/(d_s+d_f), fmax)
					v(i,j,k) = -f/(1.d0-f)*v(i,j,k-1)
				endif
			elseif (k > 1 .and. obs_lsy(i,j,k-1) > eps) then
				d_s = obs_lsy(i,j,k-1)
				d_f = max(-obs_lsy(i,j,k), eps)
				f = min(d_s/(d_s+d_f), fmax)
				v(i,j,k-1) = -f/(1.d0-f)*v(i,j,k)
			endif
		enddo
	enddo
	enddo
	!$omp end parallel do

	! componente w (nós da malha w: nx x ny x nz1)
	!$omp parallel do private(k,d_s,d_f,f)
	do j = 1, ny
	do i = 1, nx
		do k = 1, nz1
			if (obs_lsz(i,j,k) > eps) then
				w(i,j,k) = 0.d0
				if (k > 1 .and. obs_lsz(i,j,k-1) <= -eps) then
					d_s = obs_lsz(i,j,k)
					d_f = max(-obs_lsz(i,j,k-1), eps)
					f = min(d_s/(d_s+d_f), fmax)
					w(i,j,k) = -f/(1.d0-f)*w(i,j,k-1)
				endif
			elseif (k > 1 .and. obs_lsz(i,j,k-1) > eps) then
				d_s = obs_lsz(i,j,k-1)
				d_f = max(-obs_lsz(i,j,k), eps)
				f = min(d_s/(d_s+d_f), fmax)
				w(i,j,k-1) = -f/(1.d0-f)*w(i,j,k)
			endif
		enddo
	enddo
	enddo
	!$omp end parallel do

END SUBROUTINE ibm_forcado

!##################################################################################################################

!Subrotina definir IBM
!Referência: Auguste et al., 2019

!Implementação em 06/2023
!Leonardo Romero Monteiro

!############################################################

SUBROUTINE ibm_posicoes(id_ibm,obs_lsg,dimx,dimy,dimz)

!	USE disc
!	USE obst
	!$ USE omp_lib

	IMPLICIT NONE

	integer,intent(IN) :: dimx,dimy,dimz
	real(8),intent(IN),dimension(dimx,dimy,dimz) :: obs_lsg
	integer,intent(OUT),dimension(dimx,dimy,dimz) :: id_ibm
  	real(8) :: v_sol
	integer :: i,j,k
	
	! defini como zero aqui para quando fizermos objetos móveis, já ter o lugar para modificar.
	v_sol = 0.
	! zera todos os ids para iniciar sem obstáculo
	id_ibm = 0
	
	! O ideal era fazer só uma vez no começo da simulação ou quando o objeto se mexer.	
	! primeiro pensei em ser binário o id sim ou não, mas pode ser interessante que ele defina se ele usa uma célula +2 ou não
	
	!obs_lsg positivo é dentro do obstáculo. Só serão avaliadas as células de dentro do obstáculo, pois elas que serão modificadas
	! limite na direção x


	!id_ibm define a posição da célula que será modificada, ou seja, a primeira célula dentro do obstáculo. Assim, se id_ibm = 1, a célula i,j,k é a primeira dentro do obstáculo e a célula i+1,j,k é a segunda. Se id_ibm = 0, a célula i,j,k não é a primeira dentro do obstáculo e não é necessário modificar a célula i+1,j,k.

	!$omp parallel
	!$omp do
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx-2
		if(obs_lsg(i,j,k)>0) then ! garante que apenas as células de dentro do obstáculo serão marcadas
			if(obs_lsg(i,j,k)*obs_lsg(i+1,j,k)<0. .or. obs_lsg(i,j,k)*obs_lsg(i+2,j,k)<0.) then ! mudança de sinal representa presença da interface
				id_ibm(i,j,k) =  1
			endif
		endif	

	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do
	do k = 1, dimz
	do j = 1, dimy
	do i = 3, dimx
		
		if(obs_lsg(i,j,k)>0) then
			if(obs_lsg(i,j,k)*obs_lsg(i-1,j,k)<0. .or. obs_lsg(i,j,k)*obs_lsg(i-2,j,k)<0.) then
				id_ibm(i,j,k) =  1
			endif
		endif	
	enddo
	enddo
	enddo
	!$omp end do nowait

	! limite na direção y
	!$omp do
	do k = 1, dimz
	do j = 1, dimy-2
	do i = 1, dimx
		if(obs_lsg(i,j,k)>0) then
			if(obs_lsg(i,j,k)*obs_lsg(i,j+1,k)<0. .or. obs_lsg(i,j,k)*obs_lsg(i,j+2,k)<0.) then
				id_ibm(i,j,k) =  1
			endif
		endif	
	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do
	do k = 1, dimz
	do j = 3, dimy
	do i = 1, dimx
		if(obs_lsg(i,j,k)>0) then
			if(obs_lsg(i,j,k)*obs_lsg(i,j-1,k)<0. .or. obs_lsg(i,j,k)*obs_lsg(i,j-2,k)<0.) then
				id_ibm(i,j,k) =  1
			endif
		endif	
	enddo
	enddo
	enddo
	!$omp end do nowait


	! limite na direção z
	!$omp do
	do k = 1, dimz-2
	do j = 1, dimy
	do i = 1, dimx
		if(obs_lsg(i,j,k)>0) then
			if(obs_lsg(i,j,k)*obs_lsg(i,j,k+1)<0. .or. obs_lsg(i,j,k)*obs_lsg(i,j,k+2)<0.) then
				id_ibm(i,j,k) =  1
			endif
		endif	
	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do
	do k = 3, dimz
	do j = 1, dimy
	do i = 1, dimx
				
		if(obs_lsg(i,j,k)>0) then
			if(obs_lsg(i,j,k)*obs_lsg(i,j,k-1)<0. .or. obs_lsg(i,j,k)*obs_lsg(i,j,k-2)<0.) then
				id_ibm(i,j,k) =  1
			endif
		endif	

	enddo
	enddo
	enddo
	!$omp end do
    !$omp end parallel


END SUBROUTINE ibm_posicoes

!############################################################
SUBROUTINE ibm_images(obs_lsg,a1,v_sol,dimx,dimy,dimz,id_ibm,der_t)

	USE disc, only : nx,ny,nz,dx,dy,dz,dzz,deltai,ibm_t
	USE param, only : z0, cka
	USE cond, only : ccz0
	USE velpre, only : rho
	USE les, only : nut
	USE ls_param, only : mu_f2, rho_f2

	IMPLICIT NONE
	
		integer :: i,j,k,l,ll,ii,jj,kk,i1,der_t,n
		integer :: dimx,dimy,dimz
		integer,dimension(dimx,dimy,dimz) :: id_ibm
		real(8),dimension(dimx,dimy,dimz) :: obs_lsg
	
	real(8),dimension(dimx,dimy) :: v_sol
	real(8),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: a1
	
	real(8),dimension(dimx,dimy,dimz) :: a1_image, norm_x, norm_y, norm_z, mod_ls
	
	real(8),dimension(dimx,dimy,dimz) :: a1_v, norm_x1, norm_y1, norm_z1, mod_lsa
	real(8),dimension(2) :: a1_imagep
	real(8),dimension(dimx) :: x
	real(8),dimension(dimy) :: y
	real(8),dimension(dimz) :: z
		integer,dimension(8) :: ic, jc, kc
		real(8),dimension(8) :: lint
			real(8) :: lg, li1, li2, sumlint,er, xl, yl, zl, dist, image_dist
		logical :: exact_hit

	

	
	ll = 2 ! número de camadas imaginárias
	er = minval(deltai)*0.001 ! erro de aproximação da célula
	
	!v_sol = 0. !velocidade alvo
	
	!calcula a normal em todo o domínio ! coloquei o menos para que a normal saia da superfície.
	CALL mod_ls1(-obs_lsg,mod_ls,norm_x,norm_y,norm_z,dimx,dimy,dimz)
	
	!if (der_t == 2) then
	!	a1_v(1:dimx,1:dimy,1:dimz) = a1(1:dimx,1:dimy,1:dimz) ! ativado para o level_set, retira os fantasmas
	!	CALL mod_ls1(a1_v,mod_lsa,norm_x1,norm_y1,norm_z1,dimx,dimy,dimz) ! ativado para level_set, calcula normal
	!endif
	
	!CALL mod_ls1(obs_lsg,mod_lsa,norm_x1,norm_y1,norm_z1,dimx,dimy,dimz)
	
	! calcula a distância de acordo com a posição que o ponto está na malha

	if (dimx == nx) then
		do i = 1, dimx
			x(i) = (i-0.5) * dx
		enddo
	else
		do i = 1, dimx
			x(i) = (i-1.) * dx
		enddo	
	endif

	if (dimy == ny) then
		do j = 1, dimy
			y(j) = (j-0.5) * dy
		enddo
	else
		do j = 1, dimy
			y(j) = (j-1.) * dy
		enddo	
	endif

	if (dimz == nz) then
		z(1) = 0.5*dz(1)
		do k = 2, dimz
			z(k) = z(k-1)+dzz(k)
		enddo
	else
		z(1) = 0.
		do k = 2, dimz
			z(k) = z(k-1)+dz(k)
		enddo
	endif
	
	
	! inicia o cálculo das propriedades imaginárias

	!$omp parallel do private(image_dist,xl,yl,zl,ii,jj,kk,ic,jc,kc,lint,sumlint,&
	!$omp n,dist,exact_hit,li1,li2,lg,l,a1_imagep) schedule(dynamic,32)
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx

    ! inicia zerando tudo dentro do obstáculo !! não está parecendo ser necesssário
    if (obs_lsg(i,j,k) > 0. .and. id_ibm(i,j,k) .ne. 1) a1(i,j,k) = 0.
    !if (obs_lsg(i,j,k) > 0. .and. id_ibm(i,j,k) .ne. 1 .and. der_t ==2) a1(i,j,k) = -deltai(k)*2.    
       
        
        ! identifica que está nas primeiras células dentro da interface do obstáculo
	if (id_ibm(i,j,k) == 1) then 

		! se a posição da célular for muito próxima da interface do obstáculo, existe esta condição especial, sen ão a interpoação pode dar problema
		if (obs_lsg(i,j,k) < er) then
			a1_image(i,j,k) = a1(i,j,k)

        else
         	
         	! está querendo ir para baixo
			if (k == 1) norm_z(i,j,k) = norm_z(i,j,k+1)  	
		
			! cria as duas imagens que servirão para a extrapolação
			do l = 1, ll
				lint = 0.
				sumlint = 0.
				
				! cálculo da variável interpolada para as posições imaginárias 1 e 2
				!distâncias em que os valores imaginários estão localizados com relação ao ponto dentro do IBM (Figura 2)
					image_dist = l*deltai(k) + obs_lsg(i,j,k)
					xl = x(i) + image_dist*norm_x(i,j,k)
					yl = y(j) + image_dist*norm_y(i,j,k)
					zl = z(k) + image_dist*norm_z(i,j,k)
				
					!arredondou para baixo. Assim sabemos que a média deve ser feita com o +1.
					!Os índices são limitados para manter os oito vértices dentro do domínio,
					!inclusive quando a imagem fica perto de fronteiras ou em malha z variável.
					ii = lower_index_uniform(xl,x,dx,dimx)
					jj = lower_index_uniform(yl,y,dy,dimy)
					kk = lower_index_z(zl,z,dimz)

					ic = (/ ii, ii+1, ii,   ii,   ii+1, ii+1, ii,   ii+1 /)
					jc = (/ jj, jj,   jj+1, jj,   jj+1, jj,   jj+1, jj+1 /)
					kc = (/ kk, kk,   kk,   kk+1, kk,   kk+1, kk+1, kk+1 /)

					! Interpolação por distância inversa (Auguste et al., 2019) usando
					! apenas nós no fluido. id_ibm marca só a camada fantasma; obs_lsg
					! identifica todo o sólido.
					lint = 0.d0
					sumlint = 0.d0
					a1_imagep(l) = 0.d0
					exact_hit = .false.
					do n = 1, 8
						if (obs_lsg(ic(n),jc(n),kc(n)) <= 0.d0) then
							dist = dsqrt((xl-x(ic(n)))**2 + (yl-y(jc(n)))**2 + (zl-z(kc(n)))**2)
							if (dist <= er) then
								a1_imagep(l) = a1(ic(n),jc(n),kc(n))
								sumlint = 1.d0
								exact_hit = .true.
								exit
							else
								lint(n) = 1.d0/dist
								sumlint = sumlint + lint(n)
								a1_imagep(l) = a1_imagep(l) + lint(n)*a1(ic(n),jc(n),kc(n))
							endif
						endif
					enddo

					if (.not. exact_hit) then
						if (sumlint > 0.d0) then
							a1_imagep(l) = a1_imagep(l)/sumlint
						else
							a1_imagep(l) = 0.d0
						endif
					endif

				enddo
			
			!!@@## interpolador a (eq. 14 do artigo)
			li1 = (2.*deltai(k) - obs_lsg(i,j,k))/deltai(k) * (2.*obs_lsg(i,j,k))/(deltai(k) + obs_lsg(i,j,k))
			li2 = (obs_lsg(i,j,k) - deltai(k))/deltai(k) * (2.*obs_lsg(i,j,k))/(2.*deltai(k) + obs_lsg(i,j,k))
			lg  = (obs_lsg(i,j,k) - deltai(k))/(obs_lsg(i,j,k) + deltai(k))*(obs_lsg(i,j,k)-2.*deltai(k))/(obs_lsg(i,j,k)+2.*deltai(k))
		
			! cálculo das imagens que juntam as imagens 1 e 2 (eq. 12 do artigo)
			if (der_t == 1) a1_image(i,j,k) = (2.*lg*v_sol(i,j)   + li1*a1_imagep(1) + li2*a1_imagep(2))/(1.+lg)
			if (der_t == 2) a1_image(i,j,k) = (2.*lg*a1_imagep(1) + li1*a1_imagep(1) + li2*a1_imagep(2))/(1.+lg)
			
			!!@@## interpolador b
			!li1 = (2.*deltai(k) - obs_lsg(i,j,k))/deltai(k)*obs_lsg(i,j,k)/deltai(k)
			!li2 = (obs_lsg(i,j,k) - deltai(k))/deltai(k) * (obs_lsg(i,j,k))/(2.*deltai(k))
			!lg  = (obs_lsg(i,j,k))/deltai(k) *(obs_lsg(i,j,k)-2.*deltai(k))/(2.*deltai(k))
		
			! cálculo das imagens que juntam as imagens 1 e 2
			!if (der_t == 1) a1_image(i,j,k) = lg*v_sol(i,j) + li1*a1_imagep(1) + li2*a1_imagep(2)
			!if (der_t == 2) a1_image(i,j,k) = lg*a1_imagep(1) + li1*a1_imagep(1) + li2*a1_imagep(2)
			
			!!@@## interpolador c
			!li1 = 1./(deltai(k)    -obs_lsg(i,j,k))
			!li2 = 1./(2.*deltai(k) -obs_lsg(i,j,k))
			!lg  = 1./obs_lsg(i,j,k) 
		
			! cálculo das imagens que juntam as imagens 1 e 2
			!if (der_t == 1)	a1_image(i,j,k) = (lg*v_sol(i,j) + li1*a1_imagep(1) + li2*a1_imagep(2))/(li1+li2+lg)
			!if (der_t == 2) a1_image(i,j,k) = (lg*a1_imagep(1) + li1*a1_imagep(1) + li2*a1_imagep(2))/(li1+li2+lg)
    	endif


		if (der_t == 1) then
			a1(i,j,k) =  2.*v_sol(i,j) - a1_image(i,j,k)
			
		elseif (der_t == 2) then ! ativado apenas para semislip com obstáculo
			a1(i,j,k) =  a1_image(i,j,k) + 2.*v_sol(i,j)*obs_lsg(i,j,k) !!! a distância está ao quadrado pois trouxe a outra do v_sol
			
		endif
!		if (der_t == 2) then
!			norm_x1(i,j,k) = (norm_x1(i,j,k)*mod_lsa(i,j,k)*obs_lsg(i,j,k))*norm_x(i,j,k)
!			norm_y1(i,j,k) = (norm_y1(i,j,k)*mod_lsa(i,j,k)*obs_lsg(i,j,k))*norm_y(i,j,k)
!			norm_z1(i,j,k) = (norm_z1(i,j,k)*mod_lsa(i,j,k)*obs_lsg(i,j,k))*norm_z(i,j,k)	

!			a1(i,j,k) =  2.*(norm_x1(i,j,k) + norm_y1(i,j,k) + norm_z1(i,j,k)) + a1_image(i,j,k)
!		endif

	endif
		
	enddo
	enddo
	enddo
	!$omp end parallel do

	CONTAINS

		INTEGER FUNCTION lower_index_uniform(xp,coord,ds,dim)
			IMPLICIT NONE
			integer,intent(IN) :: dim
			real(8),intent(IN) :: xp, ds
			real(8),intent(IN),dimension(dim) :: coord

			lower_index_uniform = floor((xp - coord(1))/ds) + 1
			lower_index_uniform = max(1, min(dim-1, lower_index_uniform))
		END FUNCTION lower_index_uniform

		INTEGER FUNCTION lower_index_z(zp,coord,dim)
			IMPLICIT NONE
			integer,intent(IN) :: dim
			real(8),intent(IN) :: zp
			real(8),intent(IN),dimension(dim) :: coord
			integer :: kk_local

			if (zp <= coord(1)) then
				lower_index_z = 1
			elseif (zp >= coord(dim)) then
				lower_index_z = dim-1
			else
				lower_index_z = 1
				do kk_local = 1, dim-1
					if (zp >= coord(kk_local) .and. zp < coord(kk_local+1)) then
						lower_index_z = kk_local
						exit
					endif
				enddo
			endif
		END FUNCTION lower_index_z

END SUBROUTINE ibm_images
