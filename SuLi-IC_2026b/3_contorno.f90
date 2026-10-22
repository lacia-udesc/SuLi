!Subrotina definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!!!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!!!Modificações
!Leonardo Romero Monteiro em 18/08/2016
!Bruna Fernanda Soares em 23/02/2024

SUBROUTINE contorno(nlock)

	USE velpre
	USE obst
	USE cond
	USE ls_param
	USE param
	!$ USE omp_lib

	IMPLICIT NONE

	!Declarado também no programa
	integer :: i, j, k, niv, nlock
	real(8) :: zi, zj, zk

	real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: dpdx
	real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: dpdy
	real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: dpdz
	real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: ls_v

	real(8),dimension(nx1,ny) :: utau2
	real(8),dimension(nx,ny1) :: vtau2
	real(8),dimension(nx,ny)  :: wtau2


	!RESOLUÇÃO DO PROBLEMA
	if ((nlock == 1).or.(nlock == 2)) then

		
		if (ibm_t == 1 .and. obst_t /= 0) then
			! Forçagem direta (Fadlun et al., 2000):
			! zera os nós dentro do sólido e reconstrói o nó fantasma da
			! interface por interpolação linear, impondo no-slip na
			! posição exata da parede (obs_ls* = 0).
			! Para ccz0 == 4 o semi_slip() complementa a reconstrução
			! substituindo o fantasma pela lei da parede.
			CALL ibm_forcado()
		elseif (ibm_t == 2 .and. ccz0 .ne. 4) then
			
			utau2 = 0.
			vtau2 = 0.
			wtau2 = 0.
			
			!cria as imagens interpoladas
			CALL ibm_images(obs_lsx,u,utau2,nx1,ny,nz,id_ibmx,1)
			CALL ibm_images(obs_lsy,v,vtau2,nx,ny1,nz,id_ibmy,1)	
			CALL ibm_images(obs_lsz,w,wtau2,nx,ny,nz1,id_ibmz,1)
			
		endif

		!if (mms_t > 0) call mms_bc() !desatualizada
		
		CALL cc_geral(prd1,nx,ny,nz,1)



		!Semi-slip condition (Ferziger et al., 2020; Chow et al., 2005)
		!A camada limite também deve ser aplicada quando há IBM ativo, porque a interface
		!imersa passa a atuar como uma parede com efeito de rugosidade.
		! Paredes semi-slip (fundo e laterais): aplica a lei da parede nas
		! componentes tangenciais de cada parede ativa (ccz0/ccy0/ccyf == 4).
		! A rotina semi_slip() trata cada direção de forma genérica.
		if (ccz0 == 4 .or. ccy0 == 4 .or. ccyf == 4) CALL semi_slip()




		!!##!! Parede do fundo (k = 1)  !!##!!
		!Free-slip
		if (ccz0.eq.1) then
			u(:,:,0) = u(:,:,1)
			v(:,:,0) = v(:,:,1)
			w(:,:,0) = -w(:,:,2)

		!No-slip condition
		elseif (ccz0.eq.2) then
			u(:,:,0) = -u(:,:,1)
			v(:,:,0) = -v(:,:,1)
			w(:,:,0) = w(:,:,2)

		!Velocidade prescrita para manufaturada
		elseif (ccz0.eq.3) then
			u(:,:,0) = bzx0(:,:)
			v(:,:,0) = bzy0(:,:)
			w(:,:,0) = bzz0(:,:)
		!! retirei o ccz0 .eq. 4 daqui. Ele vai chamar o semislip uma vez e avaliará lá dentro quem é semislip
		endif


		!!##!! Superfície Livre (k = nz ou nz1) !!##!!
		!Sempre na superfície livre será free-slip(em teste)
		if (cczf.eq.1) then
			u(:,:,nz+1)  = u(:,:,nz)
			v(:,:,nz+1)  = v(:,:,nz)
			w(:,:,nz1+1) = -w(:,:,nz1-1)

		!Sempre na superfície livre será no-slip(em teste)
		elseif (cczf.eq.2) then
			u(:,:,nz+1)  = -u(:,:,nz)
			v(:,:,nz+1)  = -v(:,:,nz)
			w(:,:,nz1+1) = w(:,:,nz1-1)

		!Velocidade prescrita para manufaturada
		elseif (cczf.eq.3) then
			u(:,:,nz+1)  = bzxf(:,:)
			v(:,:,nz+1)  = bzyf(:,:)
			w(:,:,nz1+1) = bzzf1(:,:)
		endif
		
		!!##!! Parede esquerda (j = 1) !!##!!
		!Periodica
		if (ccy0.eq.0.and.ccyf.eq.0) then
			u(:,0,:) = u(:,ny,:)    
			v(:,0,:) = v(:,ny1-1,:)
			w(:,0,:) = w(:,ny,:)  
		  			
		!Free-slip condition
		elseif (ccy0.eq.1) then 
			u(:,0,:) = u(:,1,:)
			v(:,0,:) = -v(:,2,:)
			w(:,0,:) = w(:,1,:)

		!No-slip condition
		elseif (ccy0.eq.2) then 
			u(:,0,:) = -u(:,1,:)
			v(:,0,:) = v(:,2,:)
			w(:,0,:) = -w(:,1,:)

		!Semi-slip condition (tangenciais u e w impostas em semi_slip())
		elseif (ccy0.eq.4) then 
			v(:,0,:) = v(:,2,:)
			
		!Prescrita
		elseif (ccy0.eq.3) then
			u(:,0,:) = byx0(:,:)
			v(:,0,:) = byy0(:,:)
			w(:,0,:) = byz0(:,:)
		endif


		!!##!! Parede direita (j = ny ou ny1)  !!##!!
		!Periodica
		if (ccy0.eq.0.and.ccyf.eq.0) then
			u(:,ny+1,:)  = u(:,1,:)
			v(:,ny1+1,:) = v(:,2,:)
			w(:,ny+1,:)  = w(:,1,:)

		!Free-slip condition
		elseif (ccyf.eq.1) then
			u(:,ny+1,:)  = u(:,ny,:)
			v(:,ny1+1,:) = -v(:,ny1-1,:)
			w(:,ny+1,:)  = w(:,ny,:)
			
		!No-slip condition
		elseif (ccyf.eq.2) then 
			u(:,ny+1,:) = -u(:,ny,:)
			v(:,ny1+1,:) = v(:,ny1-1,:)
			w(:,ny+1,:) = -w(:,ny,:)

		!Semi-slip condition (tangenciais u e w impostas em semi_slip())
		elseif (ccyf.eq.4) then 
			v(:,ny1+1,:) = v(:,ny1-1,:)

		!Prescrita
		elseif (ccyf.eq.3) then 
			u(:,ny+1,:)  = byxf(:,:)
			v(:,ny1+1,:) = byyf1(:,:)
			w(:,ny+1,:)  = byzf(:,:)
		endif
		
		!!##!! Parede frente (i = 1)  !!##!!
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(0,:,:) = u(nx1-1,:,:)
			v(0,:,:) = v(nx,:,:)
			w(0,:,:) = w(nx,:,:)

		!Free-slip condition
		elseif (ccx0.eq.1) then
			u(0,:,:) = -u(2,:,:)
			v(0,:,:) = v(1,:,:)
			w(0,:,:) = w(1,:,:)

		!No-slip condition
		elseif (ccx0.eq.2) then
			u(0,:,:) = u(2,:,:)
			v(0,:,:) = -v(1,:,:)
			w(0,:,:) = -w(1,:,:)

		!Prescrita
		elseif (ccx0.eq.3) then  !! rever organização!!!
			CALL preescrita(0,dpdx)
		endif

		 !!##!! Parede de trás (i = nx ou nx1) !!##!!
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(nx1+1,:,:) = u(2,:,:)
			v(nx+1,:,:)  = v(1,:,:)
			w(nx+1,:,:)  = w(1,:,:)

		!Free-slip condition
		elseif (ccxf.eq.1) then
			u(nx1+1,:,:) = -u(nx1-1,:,:)
			v(nx+1,:,:)  = v(nx,:,:)
			w(nx+1,:,:)  = w(nx,:,:)

		!No-slip condition
		elseif (ccxf.eq.2) then
			u(nx1+1,:,:) = u(nx1-1,:,:)
			v(nx+1,:,:)  = -v(nx,:,:)
			w(nx+1,:,:)  = -w(nx,:,:)
		!Prescrita
		elseif (ccxf.eq.3) then
			u(nx1+1,:,:) = bxxf1(:,:)
			v(nx+1,:,:)  = bxyf(:,:)
			w(nx+1,:,:)  = bxzf(:,:)

		!Saida livre
		elseif (ccxf.eq.4) then
			u(nx1+1,:,:) = u(nx1,:,:)
			v(nx+1,:,:)  = v(nx,:,:)
			w(nx+1,:,:)  = w(nx,:,:)
		endif

	endif

	if (nlock == 2) then

		if (ccx0 .eq. 3 .or. ccxf .eq. 3 .or. ccy0 .eq. 3 .or. ccyf .eq. 3 .or. ccz0 .eq. 3 .or. cczf .eq. 3 .or. &
		    ccz0 .eq. 4 .or. ccy0 .eq. 4 .or. ccyf .eq. 4) then
			CALL prd_corr(dpdx,dpdy,dpdz)
		endif

		!!##!! Parede do fundo (k = 1) !!##!!
		!Free-slip
		if (ccz0.eq.1) then
			w(:,:,1) = 0.

		!No-slip condition
		elseif (ccz0.eq.2) then
			w(:,:,1) = 0.

		!Velocidade prescrita para manufaturada
		elseif (ccz0.eq.3) then
			w(:,:,1) = bzz1(:,:) + dpdz(:,:,1)
			
		!Semi-slip condition
		elseif (ccz0.eq.4) then
			w(:,:,1) = 0.
		endif			
					
		!!##!! Superfície Livre (k = nz ou nz1) !!##!!
		!Sempre na superfície livre será free-slip(em teste)
		if (cczf.eq.1) then
			w(:,:,nz1)   = 0.
			
		!Sempre na superfície livre será no-slip(em teste)
		elseif (cczf.eq.2) then
			w(:,:,nz1)   = 0.

		!Velocidade prescrita para manufaturada
		elseif (cczf.eq.3) then
			w(:,:,nz1)   = bzzf(:,:)    + dpdz(:,:,nz1)
		endif



		!!##!! Parede direita (j = 1 ) !!##!!
		!Periodica
		if (ccy0.eq.0.and.ccyf.eq.0) then
			v(:,1,:) = v(:,ny1,:)     
			
		!Free-slip condition
		elseif (ccy0.eq.1) then 
			v(:,1,:) = 0.

		!No-slip condition
		elseif (ccy0.eq.2) then 
			v(:,1,:) = 0.
			
		!Semi-slip condition (componente normal nula na parede)
		elseif (ccy0.eq.4) then 
			v(:,1,:) = 0.

		!Prescrita
		elseif (ccy0.eq.3) then
			v(:,1,:) = byy1(:,:) + dpdy(:,1,:)
		endif

		!!##!! Parede direita (j = ny ou ny1) !!##!!
		!Periodica
		if (ccy0.eq.0.and.ccyf.eq.0) then
			v(:,ny1,:)   = v(:,1,:)

		!Free-slip condition
		elseif (ccyf.eq.1) then
			v(:,ny1,:)   = 0.
			
		!No-slip condition
		elseif (ccyf.eq.2) then 
			v(:,ny1,:) = 0.

		!Semi-slip condition (componente normal nula na parede)
		elseif (ccyf.eq.4) then 
			v(:,ny1,:) = 0.

		!Prescrita
		elseif (ccyf.eq.3) then 
			v(:,ny1,:)   = byyf(:,:) + dpdy(:,ny1,:)
		
		endif

		!!##!! Parede frente (i = 1) !!##!!
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(1,:,:) = u(nx1,:,:)

		!Free-slip condition
		elseif (ccx0.eq.1) then
			u(1,:,:) = 0.

		!No-slip condition
		elseif (ccx0.eq.2) then
			u(1,:,:) = 0.

		!Prescrita
		elseif (ccx0.eq.3) then !! rever organização!!!
					
			CALL preescrita(1,dpdx)
		endif


		!!##!! Parede de trás (i = nx ou nx1) !!##!!
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(nx1,:,:)   = u(1,:,:)

		!Free-slip condition
		elseif (ccxf.eq.1) then
			u(nx1,:,:)   = 0.

		!No-slip condition
		elseif (ccxf.eq.2) then
			u(nx1,:,:)   = 0.

		!Prescrita
		elseif (ccxf.eq.3) then
			u(nx1,:,:)   = bxxf(:,:)   + dpdx(nx1,:,:)

		!Saida livre
		elseif (ccxf.eq.4) then
			u(nx1,:,:)   = u(nx1-1,:,:)
		endif
			
			
	elseif ((nlock == 3) .and. (ibm_t .ne. 0)) then

		if (ibm_t == 1) then

		!Contorno para Level Set, o ideal é que o objeto seja representado por pelo menos dois grid por direção			

			! OBS: há dependência espacial via ls(i-1,*,k) e ls(*,j-1,k); portanto não é seguro paralelizar.
			! Reordena-se apenas para melhorar acesso à memória (i como índice interno) preservando a varredura em j/i.
			do k = 1, nz-1
			do j = 2,ny-1
			do i = 2,nx-1
				if (k < kw(i,j)) then

				if ((kw(i-1,j) < k) .and. ((kw(i,j-1) < k)) ) then ! 5ª aqui faz interpolação 
					ls(i,j,k) = (ls(i-1,j,k) + ls(i,j-1,k))*0.5
				
				elseif ((kw(i+1,j) < k) .and. ((kw(i,j-1) < k)) ) then ! 6ª aqui faz interpolação 
					ls(i,j,k) = (ls(i+1,j,k) + ls(i,j-1,k))*0.5

				elseif ((kw(i+1,j) < k) .and. ((kw(i,j+1) < k)) ) then ! 7ª aqui faz interpolação 
					ls(i,j,k) = (ls(i+1,j,k) + ls(i,j+1,k))*0.5

				elseif ((kw(i-1,j) < k) .and. ((kw(i,j+1) < k)) ) then ! 8ª aqui faz interpolação 
					ls(i,j,k) = (ls(i-1,j,k) + ls(i,j+1,k))*0.5

				elseif ((kw(i-1,j) < k) .and. ((kw(i,j-1) >= k) .and. (kw(i,j+1) >= k)) ) then ! 1ª aqui faz interpolação 
					ls(i,j,k) = ls(i-1,j,k) !2.*ls(i-1,j,k) - ls(i-2,j,k)

				elseif ((kw(i,j-1) < k) .and. ((kw(i-1,j) >= k) .and. (kw(i+1,j) >= k)) ) then ! 2ª aqui faz interpolação 
					ls(i,j,k)= ls(i,j-1,k) !2.*ls(i,j-1,k) - ls(i,j-2,k)

				elseif ((kw(i+1,j) < k) .and. ((kw(i,j-1) >= k) .and. (kw(i,j+1) >= k)) ) then ! 3ª aqui faz interpolação 
					ls(i,j,k)= ls(i+1,j,k) !2.*ls(i+1,j,k) - ls(i+2,j,k)

				elseif ((kw(i,j+1) < k) .and. ((kw(i-1,j) >= k) .and. (kw(i+1,j) >= k)) ) then ! 4ª aqui faz interpolação 
					ls(i,j,k)= ls(i,j+1,k) !2.*ls(i,j+1,k) - ls(i,j+2,k)

				elseif ( (k == kw(i,j)-1) .and. (k < nz) ) then !Condição de topo
					ls(i,j,k) = ls(i,j,k+1)

				!Caso dentro do obstáculo tenha água, é para interpolar.
				!else	
				! 	ls(i,j,k) = (ls(i+1,j,k)+ls(i-1,j,k)+ls(i,j+1,k)+ls(i,j-1,k)) *0.25
				endif
				endif
			enddo
			enddo
			enddo

			! CC do Level Set nas colunas de fronteira. Nas direcoes periodicas
			! nao se impoe gradiente nulo: as celulas de contorno sao interiores
			! ao wrap periodico usado pelo WENO (ihs = 0).
			if (ccx0 .ne. 0) then
				ls(1,:,:)  = ls(2,:,:)
				ls(nx,:,:) = ls(nx-1,:,:)
			endif

			if (ccy0 .ne. 0) then
				ls(:,1,:)  = ls(:,2,:)
				ls(:,ny,:) = ls(:,ny-1,:)
			endif
			
			
		elseif (ibm_t == 22222) then !!!!! função desativada do ibm_t == 2. Quando foi avaliada a planície de inundação percebeu-se que ela maximiza um efeito possívelmente falso. O ideal é achar algum caso de validação para ter certeza.
		
			!! IBM para Level-Set
!!			CALL fantasmas_cf(ls,nx,ny,nz,ls_v)
!!			CALL ibm_images(obs_lss,ls_v,nx,ny,nz,id_ibm,2)

!!			!$omp parallel do
!!			do k = 1, nz
!!			do j = 1, ny
!!			do i = 1, nx
!!				ls(i,j,k) = ls_v(i,j,k)
!!			enddo
!!			enddo
!!			enddo
!!			!$omp end parallel do

			!CALL reinic_weno(ls,nx,ny,nz)
		endif

		CALL heaviside()
	endif

END SUBROUTINE contorno

!##################################################################################################################

SUBROUTINE contorno_les()

	USE velpre
	USE cond
	USE les
	USE obst
	
	IMPLICIT NONE
	
	integer :: i, j, k	
	
	!EIXO Z
	
	!Free-slip condition
	if (ccz0.eq.1) then
		ka(:,:,0) = ka(:,:,1)
	endif

	!No-slip condition
	if (ccz0.eq.2) then
		ka(:,:,0) = -ka(:,:,1)
	endif
	
	!Semi-slip condition
	if (ccz0.eq.4) then
		ka(:,:,0) = -ka(:,:,1)
		!ka(1:nx,1:ny,1) = 0.02*0.02*uinicial*uinicial*1.5 
	endif	
	
		
	!Free-slip
	if (cczf.eq.1) then
		ka(:,:,nz+1) = ka(:,:,nz)
	endif
	
	!Preescrita
	!if (ccy0.eq.3) then
	!	ka(1:nx,1:ny,0) = (u(1:nx,1:ny,0)*u(1:nx,1:ny,0)+v(1:nx,1:ny,0)*v(1:nx,1:ny,0)+w(1:nx,1:ny,0)*w(1:nx,1:ny,0))/6.
	!endif
	
	!if (ccyf.eq.3) then
	!	ka(1:nx,1:ny,nz+1) = (u(1:nx,1:ny,nz1+1)*u(1:nx,1:ny,nz1+1)+v(1:nx,1:ny,nz1+1)* &
	!			     v(1:nx,1:ny,nz1+1)+w(1:nx,1:ny,nz1+1)*w(1:nx,1:ny,nz1+1))/6.
	!endif


	!if (ibm_t == 2) then
		!! IBM para o ka
	!	CALL ibm_images(obs_lss,ka,0,nx,ny,nz,id_ibm,3)
	!endif

	! IBM forçado (ibm_t == 1): ka = 0 dentro do sólido (condição de parede
	! para a energia cinética turbulenta), evitando geração de nut/mut no
	! interior do obstáculo.
	if (ibm_t == 1 .and. obst_t /= 0) then
		!$omp parallel do
		do j = 1, ny
		do i = 1, nx
			do k = 1, nz
				if (obs_lss(i,j,k) > 0.) ka(i,j,k) = 0.
			enddo
		enddo
		enddo
		!$omp end parallel do
	endif

	!EIXO Y	
	
	!Periodica
	if (ccy0.eq.0.and.ccyf.eq.0) then
		ka(:,ny+1,:) = ka(:,1,:)
		ka(:,0,:) = ka(:,ny,:)
	endif
	
	!Free-slip condition
	if (ccy0.eq.1) then
		ka(:,0,:) = ka(:,1,:)
	endif

	!No-slip condition
	if (ccy0.eq.2) then
		ka(:,0,:) = -ka(:,1,:)
	endif
	
	!Semi-slip condition
	if (ccy0.eq.4) then
		ka(:,0,:) = -ka(:,1,:)
	endif
		
	!Free-slip
	if (ccyf.eq.1) then
		ka(:,ny+1,:) = ka(:,ny,:)
	endif
	
	if (ccyf.eq.2) then
		ka(:,ny+1,:) = -ka(:,ny,:)
	endif
	
	!Semi-slip condition
	if (ccyf.eq.4) then
		ka(:,ny+1,:) = -ka(:,ny,:)
	endif
	
	!Preescrita
	!if (ccy0.eq.3) then
	!	ka(1:nx,0,1:nz) = (u(1:nx,0,1:nz)*u(1:nx,0,1:nz)+v(1:nx,0,1:nz)*v(1:nx,0,1:nz)+w(1:nx,0,1:nz)*w(1:nx,0,1:nz))/6.
	!endif
	
	!if (ccyf.eq.3) then
	!	ka(1:nx,ny+1,1:nz) = (u(1:nx,ny1+1,1:nz)*u(1:nx,ny1+1,1:nz)+v(1:nx,ny1+1,1:nz) &
	!			     *v(1:nx,ny1+1,1:nz)+w(1:nx,ny1+1,1:nz)*w(1:nx,ny1+1,1:nz))/6.
	!endif	


	!EIXO X
	
	!Periodica
	if (ccx0.eq.0.and.ccxf.eq.0) then
		ka(nx+1,:,:) = ka(1,:,:)
		ka(0,:,:) = ka(nx,:,:)
	endif
	
	!Free-slip condition
	if (ccx0.eq.1) then
		ka(0,:,:) = ka(1,:,:)
	endif
	
	!Free-slip ou saída livre
	if (ccxf.eq.1 .or. ccxf.eq.4) then
		ka(nx+1,:,:) = ka(nx,:,:)
	endif

	!No-slip condition
	if (ccx0.eq.2) then
		ka(0,:,:) = -ka(1,:,:)
	endif
	
	if (ccxf.eq.2) then
		ka(nx+1,:,:) = -ka(nx,:,:)
	endif
	
	!Preescrita
	if (ccx0.eq.3) then
		ka(0,1:ny,1:nz) = iturb*iturb*uinicial*uinicial*1.5 
		
	endif
	
	if (ccxf.eq.3) then
		ka(nx+1,:,:) = ka(nx,:,:)
	endif
	
END SUBROUTINE contorno_les

!##################################################################################################################

SUBROUTINE contorno_sst()

	USE velpre
	USE cond
	USE les
	USE obst, only: obs_lss, id_ibm
	USE disc, only: nx, ny, nz, ibm_t, obst_t
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k
	logical :: near_obs
	real(8) :: nu_loc, ydist_loc, omega_w_loc
	real(8), parameter :: omega_min = 1.0d-6
	real(8), parameter :: beta1_sst = 0.075d0
	real(8), save, allocatable :: vsol0(:,:)
	

	!====================
	! Periodicidade / X
	!====================
	if (ccx0.eq.0.and.ccxf.eq.0) then
		ka(nx+1,:,:)  = ka(1,:,:)
		ka(0,:,:)     = ka(nx,:,:)
		omg(nx+1,:,:) = omg(1,:,:)
		omg(0,:,:)    = omg(nx,:,:)
	endif

	if (ccx0.eq.1) then
		ka(0,:,:)  = ka(1,:,:)
		omg(0,:,:) = omg(1,:,:)
	endif

	if (ccxf.eq.1 .or. ccxf.eq.4) then
		ka(nx+1,:,:)  = ka(nx,:,:)
		omg(nx+1,:,:) = omg(nx,:,:)
	endif

	if (ccx0.eq.2) then
		ka(0,:,:)  = -ka(1,:,:)
		omg(0,:,:) = omg(1,:,:)
		! omega na parede x=0 (no-slip)
		!$omp parallel do private(ydist_loc,nu_loc,omega_w_loc)
		do k = 1, nz
		do j = 1, ny
			ydist_loc = distw(1,j,k)
			if (ydist_loc /= ydist_loc) ydist_loc = 1.0d-12
			ydist_loc = max(ydist_loc, 1.0d-12)
			nu_loc = max(ls_mu(1,j,k)/rho(1,j,k), 1.0d-12)
			omega_w_loc = omega_wall(nu_loc, ydist_loc)
			omg(1,j,k) = max(omg(1,j,k), omega_w_loc)
		enddo
		enddo
		!$omp end parallel do
	endif

	if (ccxf.eq.2) then
		ka(nx+1,:,:)  = -ka(nx,:,:)
		omg(nx+1,:,:) = omg(nx,:,:)
		! omega na parede x=xf (no-slip)
		!$omp parallel do private(ydist_loc,nu_loc,omega_w_loc)
		do k = 1, nz
		do j = 1, ny
			ydist_loc = distw(nx,j,k)
			if (ydist_loc /= ydist_loc) ydist_loc = 1.0d-12
			ydist_loc = max(ydist_loc, 1.0d-12)
			nu_loc = max(ls_mu(nx,j,k)/rho(nx,j,k), 1.0d-12)
			omega_w_loc = omega_wall(nu_loc, ydist_loc)
			omg(nx,j,k) = max(omg(nx,j,k), omega_w_loc)
		enddo
		enddo
		!$omp end parallel do
	endif

	if (ccx0.eq.3) then
		ka(0,1:ny,1:nz) = iturb*iturb*uinicial*uinicial*1.5d0
		omg(0,1:ny,1:nz) = omg(1,1:ny,1:nz)
	endif

	if (ccxf.eq.3) then
		ka(nx+1,:,:)  = ka(nx,:,:)
		omg(nx+1,:,:) = omg(nx,:,:)
	endif

	!====================
	! Periodicidade / Y
	!====================
	if (ccy0.eq.0.and.ccyf.eq.0) then
		ka(:,ny+1,:)  = ka(:,1,:)
		ka(:,0,:)     = ka(:,ny,:)
		omg(:,ny+1,:) = omg(:,1,:)
		omg(:,0,:)    = omg(:,ny,:)
	endif

	if (ccy0.eq.1) then
		ka(:,0,:)  = ka(:,1,:)
		omg(:,0,:) = omg(:,1,:)
	endif

	if (ccyf.eq.1) then
		ka(:,ny+1,:)  = ka(:,ny,:)
		omg(:,ny+1,:) = omg(:,ny,:)
	endif

	if (ccy0.eq.2 .or. ccy0.eq.4) then
		ka(:,0,:)  = -ka(:,1,:)
		omg(:,0,:) = omg(:,1,:)
		! omega na parede y=0 (no-slip/semi-slip)
		!$omp parallel do private(ydist_loc,nu_loc,omega_w_loc)
		do k = 1, nz
		do i = 1, nx
			ydist_loc = distw(i,1,k)
			if (ydist_loc /= ydist_loc) ydist_loc = 1.0d-12
			ydist_loc = max(ydist_loc, 1.0d-12)
			nu_loc = max(ls_mu(i,1,k)/rho(i,1,k), 1.0d-12)
			omega_w_loc = omega_wall(nu_loc, ydist_loc)
			omg(i,1,k) = max(omg(i,1,k), omega_w_loc)
		enddo
		enddo
		!$omp end parallel do
	endif

	if (ccyf.eq.2 .or. ccyf.eq.4) then
		ka(:,ny+1,:)  = -ka(:,ny,:)
		omg(:,ny+1,:) = omg(:,ny,:)
		! omega na parede y=yf (no-slip/semi-slip)
		!$omp parallel do private(ydist_loc,nu_loc,omega_w_loc)
		do k = 1, nz
		do i = 1, nx
			ydist_loc = distw(i,ny,k)
			if (ydist_loc /= ydist_loc) ydist_loc = 1.0d-12
			ydist_loc = max(ydist_loc, 1.0d-12)
			nu_loc = max(ls_mu(i,ny,k)/rho(i,ny,k), 1.0d-12)
			omega_w_loc = omega_wall(nu_loc, ydist_loc)
			omg(i,ny,k) = max(omg(i,ny,k), omega_w_loc)
		enddo
		enddo
		!$omp end parallel do
	endif

	!====================
	! Z
	!====================
	if (ccz0.eq.1) then
		ka(:,:,0)  = ka(:,:,1)
		omg(:,:,0) = omg(:,:,1)
	endif

	if (cczf.eq.1) then
		ka(:,:,nz+1)  = ka(:,:,nz)
		omg(:,:,nz+1) = omg(:,:,nz)
	endif

	if (ccz0.eq.2 .or. ccz0.eq.4) then
		ka(:,:,0) = -ka(:,:,1)
		
		! omega na parede (aprox. Menter SST) usando distância ao fundo y=z(1)
		!$omp parallel do private(ydist_loc,nu_loc,omega_w_loc)
		do j = 1, ny
		do i = 1, nx
			ydist_loc = distw(i,j,1)
			if (ydist_loc /= ydist_loc) ydist_loc = 1.0d-12
			ydist_loc = max(ydist_loc, 1.0d-12)
			nu_loc = max(ls_mu(i,j,1)/rho(i,j,1), 1.0d-12)
			omega_w_loc = omega_wall(nu_loc, ydist_loc)
			omg(i,j,1) = max(omg(i,j,1), omega_w_loc)
		enddo
		enddo
		!$omp end parallel do
		omg(:,:,0) = omg(:,:,1)
	endif

	!====================
	! IBM (obstáculos internos)
	!====================
	if (ibm_t == 2 .and. obst_t /= 0) then
		! k na parede: condição Dirichlet k=0 aplicada por método das imagens
		if (.not. allocated(vsol0)) allocate(vsol0(nx,ny))
		vsol0 = 0.d0
		CALL ibm_images(obs_lss, ka, vsol0, nx, ny, nz, id_ibm, 1)

		! omega na parede: impõe limite inferior (Menter SST) nas células do fluido adjacentes ao obstáculo
		!$omp parallel do private(near_obs,ydist_loc,nu_loc,omega_w_loc)
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			near_obs = .false.
			if (obs_lss(i,j,k) <= 0.d0) then
				if (i > 1)  near_obs = near_obs .or. (obs_lss(i-1,j,k) > 0.d0)
				if (i < nx) near_obs = near_obs .or. (obs_lss(i+1,j,k) > 0.d0)
				if (j > 1)  near_obs = near_obs .or. (obs_lss(i,j-1,k) > 0.d0)
				if (j < ny) near_obs = near_obs .or. (obs_lss(i,j+1,k) > 0.d0)
				if (k > 1)  near_obs = near_obs .or. (obs_lss(i,j,k-1) > 0.d0)
				if (k < nz) near_obs = near_obs .or. (obs_lss(i,j,k+1) > 0.d0)
			elseif (id_ibm(i,j,k) == 1) then
				! primeira camada dentro do obstáculo (valor "fantasma" influencia os estênceis no fluido)
				near_obs = .true.
			endif

			if (near_obs) then
				ydist_loc = distw(i,j,k)
				if (ydist_loc /= ydist_loc) ydist_loc = 1.0d-12
				ydist_loc = max(ydist_loc, 1.0d-12)
				nu_loc = max(ls_mu(i,j,k)/rho(i,j,k), 1.0d-12)
				omega_w_loc = omega_wall(nu_loc, ydist_loc)
				omg(i,j,k) = max(omg(i,j,k), omega_w_loc)
			endif
		enddo
		enddo
		enddo
		!$omp end parallel do
	endif

	!====================
	! IBM forçado (ibm_t == 1): classificação das células via obs_lss
	!====================
	if (ibm_t == 1 .and. obst_t /= 0) then
		!$omp parallel do private(k,ydist_loc,nu_loc,omega_w_loc)
		do j = 1, ny
		do i = 1, nx
			! zera ka dentro do sólido (evita geração de nut/mut no sólido)
			do k = 1, nz
				if (obs_lss(i,j,k) > 0.d0) ka(i,j,k) = 0.d0
			enddo
			! célula fantasma da interface (1ª célula sólida abaixo da parede)
			do k = 1, nz-1
				if (obs_lss(i,j,k) > 0.d0 .and. obs_lss(i,j,k+1) <= 0.d0) then
					ka(i,j,k)  = -ka(i,j,k+1)
					omg(i,j,k) = omg(i,j,k+1)

					! condição de omega na parede (Menter SST) no 1º nó fluido
					ydist_loc = distw(i,j,k+1)
					if (ydist_loc /= ydist_loc) ydist_loc = 1.0d-12
					ydist_loc = max(ydist_loc, 1.0d-12)
					nu_loc = max(ls_mu(i,j,k+1)/rho(i,j,k+1), 1.0d-12)
					omega_w_loc = omega_wall(nu_loc, ydist_loc)
					omg(i,j,k+1) = max(omg(i,j,k+1), omega_w_loc)
				endif
			enddo
		enddo
		enddo
		!$omp end parallel do
	endif

	!====================
	! Robustez numérica
	!====================
	! Garante omega finito/positivo também em células fantasmas
	!$omp parallel do
	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, nx+1
		if (omg(i,j,k) /= omg(i,j,k)) then
			omg(i,j,k) = omega_min
		else
			omg(i,j,k) = max(omg(i,j,k), omega_min)
		endif
	enddo
	enddo
	enddo
	!$omp end parallel do
	! k só é limitado no interior; no sólido mantém-se o fantasma
	! anti-simétrico definido acima (obs_lss sempre inicializado).
	!$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		if (ka(i,j,k) /= ka(i,j,k)) ka(i,j,k) = 0.d0
		if (obs_lss(i,j,k) <= 0.d0) ka(i,j,k) = max(ka(i,j,k), 0.d0)
	enddo
	enddo
	enddo
	!$omp end parallel do

CONTAINS

	REAL(8) FUNCTION omega_wall(nu_in, ydist_in)

		IMPLICIT NONE
		REAL(8), INTENT(IN) :: nu_in, ydist_in
		REAL(8) :: ysafe
		ysafe = max(ydist_in, 1.0d-12)
		omega_wall = 60.d0*max(nu_in, 1.0d-12)/(beta1_sst*ysafe*ysafe)
	END FUNCTION omega_wall

END SUBROUTINE contorno_sst

!##################################################################################################################

SUBROUTINE sponge_layer(epis_z)

!Subrotina para representar a camada esponja, quando se quer um domínio com saída livre.
!Não funciona com escoamentos, porque a velocidade é forçada a zero.

	USE velpre
	USE ls_param
	!$ USE omp_lib

	IMPLICIT NONE

	!Declarado também no programa

	real(8),intent(out),dimension(nx,ny,nz1) :: epis_z

	!Declarado apenas na rotina

	!Profundidade do domínio
	real(8),dimension(nx) :: zpfs
	real(8),dimension(nx) :: xp
        real(8),dimension(nx) :: zpfs_local

	!Parâmetro da camada esponja (permeabilidade)
	real(8),save :: alfa, alfa0, alfa1, aux1, alt1, alt2

	!Comprimento da camada esponja e posição do seu começo
	real(8),save :: l_sponge, x_inicial
	integer :: i, j, k

	!Inicialização do código
	epis_z = 0.
	l_sponge  = 1. !Comprimento da camada esponja
	x_inicial = nx*dx - l_sponge

	alfa  = 60.          !A ser calibrado
	alfa0 = 0.           !A ser calibrado
	alfa1 = 30. !2.*alfa !A ser calibrado

	!Resolução do problema

	if (esp_type == 1) then

		! A redução de max por coluna i é exata (max é associativo/independe
		! da ordem); paraleliza por i mantendo a varredura (j,k) por coluna.
		!$omp parallel do private(j,k,alt1,alt2)
		do i = nint(x_inicial/dx)+1, nx
		do j = 1, ny
		do k = 1, nz-1
		
		if (ls(i,j,k)*ls(i,j,k+1) <= 0) then
			alt1 = z(k) + ls(i,j,k)
			alt2 = z(k+1) + ls(i,j,k+1)
			zpfs(i) = max(alt1,alt2) ! altura da interface
		endif
		enddo
		enddo
		enddo
		!$omp end parallel do

	!Camada esponja para a direção x

		!$omp parallel do private(j,k,aux1)
		do i = nint(x_inicial/dx)+1, nx
		
		xp(i) = real(i-0.5)*dx
		aux1 = (xp(i)-x_inicial)/l_sponge

		do k = 1, nz
		do j = 1, ny

                if (ls(i,j,k) >= 0.) then !Na água 
                        epis_z(i,j,k) = alfa * aux1 * aux1 * (0.-z(k)) / (0.-zpfs(i))
                else
                        epis_z(i,j,k) = alfa * aux1 * aux1 * (z(k)-lz) / (zpfs(i)-lz)
                endif

		enddo
		enddo
		enddo
		!$omp end parallel do

	elseif (esp_type == 2) then


        !$omp parallel do private(aux1)
		do k = 1, nz
		do j = 1, ny
		do i = nint(x_inicial/dx)+1, nx
		xp(i) = real(i-0.5)*dx
		aux1 = (xp(i)-x_inicial)/l_sponge

		epis_z(i,j,k) = alfa0 + aux1 * (alfa1 - alfa0)

		enddo
		enddo
		enddo
        !$omp end parallel do

	endif

END SUBROUTINE sponge_layer

!##################################################################################################################

SUBROUTINE prd_corr(dpdx,dpdy,dpdz) !! arrumar rotina para eficiência!!
!Derivadas das pressões para adicionar nas condições de contorno (aproximar o valor em u^n+1 ...)

	USE velpre
	USE param
	!$ USE omp_lib

	IMPLICIT NONE
	!Declarado também no programa
	real(8),intent(out),dimension(0:nx1+1,0:ny+1,0:nz+1) :: dpdx
	real(8),intent(out),dimension(0:nx+1,0:ny1+1,0:nz+1) :: dpdy
	real(8),intent(out),dimension(0:nx+1,0:ny+1,0:nz1+1) :: dpdz

	real(8),dimension(nx1,ny,nz) :: rhox
	real(8),dimension(nx,ny1,nz) :: rhoy
	real(8),dimension(nx,ny,nz1) :: rhoz

	integer :: i, j, k

	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)


    !$omp parallel
    !$omp do
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
	        dpdz(i,j,k) = (prd(i,j,k) - prd(i,j,k-1))/ rhoz(i,j,k) *dt /dzz(k)
	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
	        dpdy(i,j,k) = (prd(i,j,k) - prd(i,j-1,k))/ rhoy(i,j,k) *dt /dy
	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
	        dpdx(i,j,k) = (prd(i,j,k) - prd(i-1,j,k))/ rhox(i,j,k) *dt /dx
	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel


	CALL cc_geral(dpdx,nx1,ny,nz,2)
	CALL cc_geral(dpdy,nx,ny1,nz,3)
	CALL cc_geral(dpdz,nx,ny,nz1,1)


END SUBROUTINE prd_corr

!##################################################################################################################
!Reference: Coubilla, 2015 (Thesis)
SUBROUTINE boundary_waves()

	USE wave_c
	USE param
	USE velpre
	USE ls_param
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k, iik
	real(8),save :: aux1, aux2, aux3, aux4, aux5, h_fa, l_wa

	real(8),dimension(0:1) :: h_f

	real(8),dimension(0:1,0:ny+1,0:nz+1) :: u1
	real(8),dimension(0:1,0:ny1+1,0:nz+1) :: v1
	real(8),dimension(0:1,0:ny+1,0:nz1+1) :: w1
!===================================================================
	real(8),dimension(0:1,0:ny+1,0:nz+1) :: ls1

	u1 = 0.
	v1 = 0.
	w1 = 0.
	ls1 = 0.
	i = 1
	

        !$omp parallel do
        do k = 1, nz
        do j = 1, ny
                ls1(i,j,k) = ls(i,j,k)
        enddo
        enddo
        !$omp end parallel do

	ls1(0,:,:)    = ls1(1,:,:)

	ls1(:,0,:)    = ls1(:,1,:)
	ls1(:,ny+1,:) = ls1(:,ny,:)

	ls1(:,:,0)    = ls1(:,:,1)
	ls1(:,:,nz+1) = ls1(:,:,nz)

	!Stokes I
	if (wave_t == 1) then


!       !$omp parallel do
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, 1
		h_f(i)     = h0_f   + a_w*dcos(n_w*dx*(i-0.5)-f_w*t)
		ls1(i,j,k) = h_f(i) - kp(k)

		if (ls1(i,j,k) >= 0.) then
			u1(i,j,k) = a_w * f_w * dcosh(n_w*kp(k))/dsinh(n_w*h0_f) * dcos(n_w*dx*(i-1.)-f_w*t)
			w1(i,j,k) = a_w * f_w * dsinh(n_w*(kp(k)-0.5*dzz(k)))/dsinh(n_w*h0_f) * dsin(n_w*dx*(i-0.5)-f_w*t)
			
			iik  = k
		else
			u1(i,j,k) = u1(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
			w1(i,j,k) = w1(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
		endif

		enddo
		enddo
		enddo
!       !$omp end parallel do


	!Stokes II
	elseif (wave_t == 2 ) then

		aux1 = dtanh(n_w*h0_f) !Sigma

        !!$omp parallel do
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, 1

		aux2 = n_w*dx*(i-0.5)-f_w*t !Ômega para w e eta
		aux3 = n_w*dx*(i-1.)-f_w*t !Ômega para u

		h_f(i)     = h0_f + a_w*dcos(aux2) +n_w*a_w*a_w*(3.-aux1*aux1)/(4.*aux1*aux1*aux1)*dcos(2.*aux2)
		ls1(i,j,k) = h_f(i) - kp(k)

		if (ls1(i,j,k) >= 0.) then
			u1(i,j,k) = a_w*f_w*dcosh(n_w*kp(k))/dsinh(n_w*h0_f)*dcos(aux3) &
			+3./4.*a_w*a_w*f_w*n_w*dcosh(2.*n_w*kp(k))/((dsinh(n_w*h0_f))**4.)*dcos(2*aux3)

			w1(i,j,k) = a_w*f_w*dsinh(n_w*(kp(k)-0.5*dzz(k)))/dsinh(n_w*h0_f)*dsin(aux2) &
			+3./4.*a_w*a_w*f_w*n_w*dsinh(2.*n_w*(kp(k)-0.5*dzz(k)))/((dsinh(n_w*h0_f))**4.)*dsin(2.*aux3)

			iik  = k
		else
			u1(i,j,k) = u1(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
			w1(i,j,k) = w1(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
		endif
		
		enddo
		enddo
		enddo
        !!$omp end parallel do



			
	!Stokes V
	elseif (wave_t == 5 ) then

        !!$omp parallel do
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, 1

		aux2 = n_w*dx*(i-0.5)-f_w*t !Ômega para w e eta
		aux3 = n_w*dx*(i-1.)-f_w*t !Ômega para u

		aux4 = n_w*kp(k) !Para u
		aux5 = n_w*(kp(k)-0.5*dzz(k)) !Para w

		h_f(i) = h0_f +aeta1*dcos(aux2) +aeta2*dcos(2.*aux2) +aeta3*dcos(3.*aux2) +aeta4*dcos(4.*aux2) +aeta5*dcos(5.*aux2)

		ls1(i,j,k) = h_f(i) - kp(k)

		if (ls1(i,j,k) >= 0.) then
			u1(i,j,k) = avel1*dcosh(aux4)*dcos(aux3) +avel2*dcosh(2*aux4)*dcos(2*aux3) +avel3*dcosh(3*aux4)*dcos(3*aux3) &
			+ avel4*dcosh(4*aux4)*dcos(4*aux3) +avel5*dcosh(5*aux4)*dcos(5*aux3)

			w1(i,j,k) = avel1*dsinh(aux5)*dsin(aux2) +avel2*dsinh(2*aux5)*dsin(2*aux2) +avel3*dsinh(3*aux5)*dsin(3*aux2) &
			+ avel4*dsinh(4*aux5)*dsin(4*aux2) +avel5*dsinh(5*aux5)*dsin(5*aux2)

			iik  = k
		else
			u1(i,j,k) = u1(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
			w1(i,j,k) = w1(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
		endif

		enddo
		enddo
		enddo
        !!$omp end parallel do

	endif

! ==========================================
	blx1(:,:) = ls(1,:,:)
	bxx0(:,:) = u1(0,:,:)
	bxy0(:,:) = v1(0,:,:)
	bxz0(:,:) = w1(0,:,:)
	bxx1(:,:) = u1(1,:,:)
! ==========================================


!        !$omp parallel
!        !$omp do
!        do k = 1, nz
!        do j = 1, ny
!        do i = 1, nx
!                ls(i,j,k) = ls1(i,j,k)
!        end do
!        end do
!        end do
!        !$omp end do nowait

!        !$omp do
!        do k = 1, nz
!        do j = 1, ny
!                bxx0(j,k) = u1(0,j,k)
!                bxy0(j,k) = v1(0,j,k)
!                bxz0(j,k) = w1(0,j,k)
!        end do
!        end do
!        !$omp end do nowait

!        !$omp do
!        do k = 1, nz
!        do j = 1, ny
!                bxx1(j,k) = u1(1,j,k)
!        end do
!        end do
!        !$omp end do
!        !$omp end parallel

	!Streamfunction

	!Solitary wave

END SUBROUTINE boundary_waves

!##################################################################################################################

SUBROUTINE waves_coef()
	USE wave_c
	USE param
	IMPLICIT NONE

	integer :: i, k, ii, nxii
	real(8),save :: aux1, aux2, aux3, aux4, aux5, h_fa, l_wa, lamb, lamb1, lamb2, erro, erro0, l_w1
	real(8),save :: s, c, a11, a13, a15, a22, a24, a33, a35, a44, a55, b22, b24, b33, b35, b44, b55, c1, c2

	real(8),dimension(nx) :: h_f

	!Reference: Coubilla, 2015 (Thesis)
	do k = 1, nz
		kp(k) = z(k)
	enddo
	kp(0)    = -dz(1)*0.5 
	kp(nz+1) = z(nz) + dz(nz)
	

	a_w  = ampl !Amplitude da onda
	p_w  = 1.  !Período da onda 
	h0_f = prof   !Profundidade do escoamento sem onda
	h_f  = h0_f

	l0_w = gz * p_w * p_w / (2.*pi)

	!Stokes I e Stokes II
	if (wave_t <= 2) then

		l_w = gz * p_w * p_w / (2.*pi)
		do i = 1, 1000
		l_w = gz * p_w * p_w / (2.*pi) * dtanh(2.*pi*h0_f /l_w)
		enddo
		c_w = l_w   / p_w !Celeridade
		f_w = 2.*pi / p_w !Frequencia angular
		n_w = 2.*pi / l_w !Número de onda
		
	elseif (wave_t == 5 ) then

		erro0 = 999.
		nxii = 10000. ! número de intervalos
		do ii = 1, nxii
		lamb1 = 0.
		lamb2 = 0.
		l_w1 = l0_w* 0.5 +l0_w*(real(ii)/nxii) *2.5
		s = dsinh(2.*pi*h0_f/l_w1)
		c = dcosh(2.*pi*h0_f/l_w1)

		b33 = 3.*(8.*c**6. +1.)/(64.*s**6.)
		b35 = (88128.*c**14. -208224.*c**12 +70848.*c**10. +54000.*c**8. -21816.*c**6. +6264.*c**4. -54.*c*c -81.)/&
		(12288.*s**12. * (6.*c*c -1.))
		b55 = (192000.*c**16. -262720.*c**14. +83680.*c**12 +20160.*c**10. -7280.*c**8. + 7160.*c**6. -1800.*c**4. -1050.*c*c +225.)/&
		(12288*s**10. * (6.*c*c -1.) * (8.*c**4. -11.*c*c +3.)) 

		c1 = (8.*c**4. -8.*c*c +9.)/(8.*s**4.)
		c2 = (3840.*c**12. -4096.*c**10. +2592.*c**8. -1008.*c**6. +5944.*c**4. -1830.*c*c +147.)/(512.*s**10. * (6.*c*c -1.))

		!Isolado1
		do i = 1, 10000
		lamb1 = pi*a_w*2. / l_w1  - lamb1**3.*b33 -lamb1**5.*(b35+b55) 
		enddo
		!Isolado2
		do i = 1, 10000
		lamb2 = dsqrt((l_w1 /(l0_w*dtanh(2.*pi*h0_f/l_w1)) -1.) / (c1 +lamb2**2.*c2))
		enddo

		erro = dabs(lamb1 - lamb2)

		if (erro < erro0) then
			erro0 = erro
			lamb = lamb1
			l_w = l_w1
		endif
		enddo

		!Aqui já temos lamb e l_w
		c_w = l_w   / p_w !celeridade
		f_w = 2.*pi / p_w !frequencia angular
		n_w = 2.*pi / l_w !número de onda

		s = dsinh(2.*pi*h0_f/l_w)
		c = dcosh(2.*pi*h0_f/l_w)

		b33 = 3.*(8.*c**6. +1.)/(64.*s**6.)
		b35 = (88128.*c**14. -208224.*c**12 +70848.*c**10. +54000.*c**8. -21816.*c**6. +6264.*c**4. -54.*c*c -81.)/&
			(12288.*s**12. * (6.*c*c -1.))
		b55 = (192000.*c**16. -262720.*c**14. +83680.*c**12 +20160.*c**10. -7280.*c**8. + 7160.*c**6. -1800.*c**4. -1050.*c*c +225.)/&
			(12288*s**10. * (6.*c*c -1.) * (8.*c**4. -11.*c*c +3.)) 

		c1 = (8.*c**4. -8.*c*c +9.)/(8.*s**4.)
		c2 = (3840.*c**12. -4096.*c**10. +2592.*c**8. -1008.*c**6. +5944.*c**4. -1830.*c*c +147.)/(512.*s**10. * (6.*c*c -1.))

		a11 = 1./s
		a13 = -c*c*(5.*c*c +1.) / (8.*s**5.)
		a15 = -(1184.*c**10. -1440.*c**8. -1992.*c**6. +2641.*c**4. -249.*c**2. +18.)/(1536.*s**11.)
		a22 = 3./(8.*s**4.)
		a24 = (192.*c**8. -424.*c**6. -312.*c**4. +480.*c*c -17.)/(768.*s**10.)
		a33 = (13.-4.*c*c)/(64.*s**7.)
		a35 = (512.*c**12. +4224.*c**10. -6800.*c**8. -12808.*c**6. +16704.*c**4. -3154.*c*c +107.)/(4096.*s**13. * (6.*c*c -1.))
		a44 = (80.*c**6. -816.*c**4. +1338.*c*c -197.)/(1536.*s**10. * (6.*c*c -1.))
		a55 = -(2880.*c**10.-72480.*c**8.+324000.*c**6.-432000.*c**4.+163470.*c*c-16245.)&
			/(61440.*s**11.*(6.*c*c-1.)*(8.*c**4.-11.*c*c+3.))
		b22 = c*(2.*c*c +1.)/(4.*s**3.)
		b24 = c*(272.*c**8. -504.*c**6. -192.*c**4. +322*c*c +21.)/(384.*s**9.)
		b44 = c*(768.*c**10. -448.*c**8. -48.*c**6. +48.*c**4. +106.*c*c -21.)/(384.*s**9. *(6.*c*c -1.)) 	

		aeta1 = lamb/n_w
		aeta2 = (lamb*lamb*b22 + lamb**4.*b24)/n_w
		aeta3 = (lamb**3.*b33  + lamb**5.*b35)/n_w
		aeta4 = lamb**4.*b44/n_w
		aeta5 = lamb**5.*b55/n_w

		avel1 = 1.*2.*pi/(p_w*n_w) * (lamb*a11     + lamb**3.*a13 + lamb**5.*a15)
		avel2 = 2.*2.*pi/(p_w*n_w) * (lamb**2.*a22 + lamb**4.*a24)
		avel3 = 3.*2.*pi/(p_w*n_w) * (lamb**3.*a33 + lamb**5.*a35)
		avel4 = 4.*2.*pi/(p_w*n_w) * (lamb**4.*a44)
		avel5 = 5.*2.*pi/(p_w*n_w) * (lamb**5.*a55)

	endif

	!write(*,*) erro0, lamb, h0_f/l_w
	!write(*,*) b22, b24
	!Streamfunction
	!Solitary wave

END SUBROUTINE waves_coef



!##################################################################################################################

SUBROUTINE preescrita(i,dpdx)

	USE velpre
	USE obst
	USE param
	USE ls_param
	
	!$ USE omp_lib

	IMPLICIT NONE

	!Declarado também no programa
	integer :: i, j, k, ik, iik
	real(8), save :: coef, ta, umed, uaux, rr
	real(8),dimension(nx1,ny,nz) :: ls_x
	real(8),dimension(0:ny+1,0:nz+1) :: bxx
	real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: dpdx
	
	if (i==0) then
		bxx(:,:) = bxx0(:,:)		
	elseif (i == 1) then
		bxx(:,:) = bxx1(:,:) + dpdx(1,:,:)
	endif
		
		
	if (wave_t > 0) then
	
		call boundary_waves() !For wave propagation
		u(i,:,:) = bxx(:,:) 
		
	else
			
			
		! Para o caso da enxurrada que entra com um tempo de ascensão
		ta = 1. !(tempo de ascensão)
		if (t .le. ta) then
			coef = 1.*it*dt/ta 
		else
			coef = 1.
		endif
	
		!Velocidade com incremento de turbulência *********************		
		call interpx_cf(ls,nx,ny,nz,ls_x) !(nx1,ny1,nz)

		ik = 0
		iik = 0
		umed = 0.

		uaux = sqrt(9.806*(ampl+prof))*1.5

		ta = 1. !(tempo de ascensão)
		if (t .le. ta) then
			coef = 1.*it*dt/ta 
		else
			coef = 1.
		endif

		do k = 1, nz
		do j = 1, ny
			if (ls_x(i,j,k)>=-2*dz(k) .and. k > ku(1,j)) then
			
				u(i,j,k) = bxx(j,k)
				
				call random_number(rr)

				rr = 2.*(rr-0.5)
				u(i,j,k) = u(i,j,k)*(1.+rr*iturb)+ (uaux-u(i,j,k))*coef
				iik  = k
				
			elseif (ls_x(i,j,k)<-2*dz(k)) then
				u(i,j,k) = u(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
			endif
		enddo
		enddo


		if (i == 0) then
			!$omp parallel do private(rr)
			do k = 1, nz
			do j = 1, ny
			
				if (ls_x(i,j,k)>=-2*dz(k) .and. k > ku(1,j)) then
					
					v(i,j,k) = bxy0(j,k)
					call random_number(rr)
					rr = 2.*(rr-0.5)
					v(i,j,k) = v(i,j,k)+rr*iturb*u(i,j,k)
									
					w(i,j,k) = bxz0(j,k)
							
					call random_number(rr)
					rr = 2.*(rr-0.5)
					w(i,j,k) = w(i,j,k)+rr*iturb*u(i,j,k)
				
				endif		
				
			enddo
			enddo
			!$omp end parallel do
		endif
		
		v(i,:,:) = bxy0(:,:)
		w(i,:,:) = bxz0(:,:)

	endif

	
	if (i == 1) then
	do k = 1, nz
		ls(i,:,k) = blx1(:,k)
	enddo 

	ls(i,:,:) = blx1(:,:)
	
	endif
	
END SUBROUTINE preescrita


	
!##################################################################################################################

SUBROUTINE cc_geral(a,dimx,dimy,dimz,tipo)

	USE cond
	!$ USE omp_lib

	IMPLICIT NONE

	!Declarado também no programa
	integer :: i, j, k, dimx, dimy, dimz,tipo
	real(8),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: a
			

	!$omp parallel
	if (ccx0.eq.0) then  ! Condição periódica

		if (tipo .ne. 2) then
			!$omp do
			do k = 0, dimz+1
			do j = 0, dimy+1		
				a(0,j,k)    = a(dimx,j,k)
				a(dimx+1,j,k) = a(1,j,k)
			enddo
			enddo
			!$omp end do
		else
			!$omp do
			do k = 0, dimz+1
			do j = 0, dimy+1	
				a(0,j,k) = a(dimx-1,j,k)
				a(1,j,k) = a(dimx,j,k)
				a(dimx+1,j,k) = a(2,j,k)
			enddo
			enddo
			!$omp end do
		endif
	
	else
        !$omp do
		do k = 0, dimz+1
		do j = 0, dimy+1	
			a(0,j,k)    = a(1,j,k)
			a(dimx+1,j,k) = a(dimx,j,k)
		enddo
		enddo
		!$omp end do
	endif

	if (ccy0.eq.0) then  ! Condição periódica
		if (tipo .ne. 3) then
			!$omp do
			do k = 0, dimz+1	
			do i = 0, dimx+1
				a(i,0,k)    = a(i,dimy,k)
				a(i,dimy+1,k) = a(i,1,k)
			enddo
			enddo
			!$omp end do
		else	
			!$omp do
			do k = 0, dimz+1	
			do i = 0, dimx+1
				a(i,0,k)    = a(i,dimy-1,k)
				a(i,1,k)    = a(i,dimy,k)
				a(i,dimy+1,k) = a(i,2,k)
			enddo
			enddo
			!$omp end do
		endif
	else
        !$omp do
		do k = 0, dimz+1
		do i = 0, dimx+1
			a(i,0,k)    = a(i,1,k)
			a(i,dimy+1,k) = a(i,dimy,k)
		enddo
		enddo
		!$omp end do
	endif

    !$omp do
	do j = 0, dimy+1
	do i = 0, dimx+1
		a(i,j,0)    = a(i,j,1)
		a(i,j,dimz+1) = a(i,j,dimz)
	enddo
	enddo
	!$omp end do
	!$omp end parallel



END SUBROUTINE cc_geral

!##################################################################################################################

SUBROUTINE semi_slip()

	USE velpre
	USE obst
	USE param
	USE les
	USE ls_param
	USE cond
	!$ USE omp_lib

	IMPLICIT NONE

	!Declarado também no programa
	integer :: i, j, k, k0, t_tau, ki
	real(8), save :: cd
	
	real(8),dimension(nx1,ny,nz) :: vx
	real(8),dimension(nx,ny1,nz) :: uy
		
	real(8),dimension(nx1,ny) :: xtau, utau, utau_ibm, lyx
	real(8),dimension(nx,ny1) :: ytau, vtau, vtau_ibm, lyy
	real(8),dimension(nx,ny)  :: ztau, wtau

	! auxiliares da reconstrução do nó fantasma no IBM (forçagem direta + lei da parede)
	real(8) :: d_s, d_f, f, u_w, v_w

	real(8),dimension(nx1,ny,nz) :: mux, rhox
	real(8),dimension(nx,ny1,nz) :: muy, rhoy
	real(8),dimension(nx,ny,nz1) :: muz, rhoz
	real(8),dimension(nx1,ny,nz) :: mux_mol
	real(8),dimension(nx,ny1,nz) :: muy_mol
	real(8),dimension(nx,ny,nz1) :: muz_mol

	real(8), parameter :: eps_min = 1.0d-10

	utau = 1.0d-6  ! Inicializar com valor pequeno positivo para evitar dlog(0)
	vtau = 1.0d-6
	wtau = 0.0d0

	utau_ibm = 1.0d-6
	vtau_ibm = 1.0d-6

	! levando os rhos e mus para as faces para trabalhar com as velocidades
	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)
	call interpx_cf(ls_mu,nx,ny,nz,mux_mol)
	call interpy_cf(ls_mu,nx,ny,nz,muy_mol)
	call interpz_cf(ls_mu,nx,ny,nz,dz,muz_mol)

	call interpx_cf(mu,nx,ny,nz,mux) !(nx1,ny,nz)
	call interpy_cf(mu,nx,ny,nz,muy) !(nx,ny1,nz)
	call interpz_cf(mu,nx,ny,nz,dz,muz) !(nx,ny,nz1)



	! Aplica semi-slip por direção normal da parede. O seletor mapeia
	! automaticamente as duas componentes tangenciais, sem espalhar a
	! associação u/v/w pela rotina principal.
	if (ccz0 == 4) CALL semi_slip_wall(3, 0, dz(1), u, v, w, rhox, rhoy, rhoz, mux_mol, muy_mol, muz_mol, &
	                                    mux, muy, muz, nx, ny, nz, nx1, ny1, nz1)
	if (ccy0 == 4) CALL semi_slip_wall(2, 0, dy, u, v, w, rhox, rhoy, rhoz, mux_mol, muy_mol, muz_mol, &
	                                    mux, muy, muz, nx, ny, nz, nx1, ny1, nz1)
	if (ccyf == 4) CALL semi_slip_wall(2, 1, dy, u, v, w, rhox, rhoy, rhoz, mux_mol, muy_mol, muz_mol, &
	                                    mux, muy, muz, nx, ny, nz, nx1, ny1, nz1)



	! Distância da rugosidade de fato até o primeiro ponto de velocidade
	! acima do IBM (necessária à reconstrução com lei da parede da interface
	! imersa, aplicada apenas quando o fundo é semi-slip).
	if (ccz0 == 4 .and. (ibm_t == 1 .or. ibm_t == 2)) then

		!$omp parallel
		!$omp do private(ki)
		do j=1,ny
		do i=1,nx1
			ki = ku(i,j)
			if (ki > 0 .and. ki < nz) then !! distância do 1º nó fluido até a parede
				lyx(i,j) = dabs(obs_lsx(i,j,ki+1)) + z0 ! distância efetiva (a partir da rugosidade)
			elseif (ki == nz) then !! caso bata no "teto"
				lyx(i,j) = dabs(obs_lsx(i,j,ki)) + z0
			else !! caso não tenha ibm nesta posição
				lyx(i,j) = dz(1)*0.5d0 + z0
			endif
		enddo
		enddo
		!$omp end do nowait
					
		!$omp do private(ki)
		do j=1,ny1
		do i=1,nx
			ki = kv(i,j)
			if (ki > 0 .and. ki < nz) then
				lyy(i,j) = dabs(obs_lsy(i,j,ki+1)) + z0 !! distância efetiva do 1º nó fluido até a parede
			elseif (ki == nz) then !! caso bata no "teto"
				lyy(i,j) = dabs(obs_lsy(i,j,ki)) + z0
			else !! caso não tenha ibm nesta posição
				lyy(i,j) = dz(1)*0.5d0 + z0
			endif
		enddo
		enddo
		!$omp end do
		!$omp end parallel

		! calculo de utau_ibm/vtau_ibm (em cima do IBM)
	CALL lei_parede_ibm(lyx, lyy, utau, vtau, utau_ibm, vtau_ibm,rhox,rhoy,mux_mol,muy_mol,mux,muy)

	endif

	if (ccz0 == 4 .and. ibm_t == 1) then

		! Fantasma com lei da parede (semi-slip, Chow et al., 2005; Ferziger et al., 2020):
		! o nó fantasma (k = ku) é reconstruído por interpolação linear na
		! interface de modo que a velocidade interpolada na parede seja
				!   tau_w = rho*u_tau*|u_tau|; extrapola-se u_w do primeiro
				!   ponto fluido até a interface usando a viscosidade efetiva.
		! (escorregamento da lei logarítmica avaliada no 1º nó fluido).
		! O fantasma de ibm_forcado() (no-slip puro) é substituído aqui.
		!$omp parallel
		!$omp do private(ki,d_s,d_f,f,u_w)
		do j=1,ny
		do i=1,nx1
			ki = ku(i,j)
			if (ki > 0 .and. ki < nz) then
				d_s = max(obs_lsx(i,j,ki),    eps_min) ! distância do nó sólido até a parede
				d_f = max(-obs_lsx(i,j,ki+1), eps_min) ! distância do 1º nó fluido até a parede
				f   = d_s/(d_s+d_f)                       ! fração da parede a partir do nó sólido
				f   = min(f, 0.9d0)                       ! limita a amplificação 1/(1-f) (estabilidade do IBM1)
				! Converte a tensão de parede em gradiente usando a mesma
				! viscosidade dinâmica efetiva usada no operador de momento.
				u_w = u(i,j,ki+1) - abs(utau_ibm(i,j))*utau_ibm(i,j)*dz(ki)*rhox(i,j,ki) / &
				      max(mux(i,j,ki),eps_min)  ! velocidade de escorregamento na parede
				u(i,j,ki) = (u_w - f*u(i,j,ki+1))/(1.d0-f)
				!u(i,j,ki) = -u(i,j,ki+1) +2.*0.043*0.043*dz(ki)/mux(i,j,ki)*rhox(i,j,ki) !-u(i,j,ki+1) + 2.*0.043*(1.0d0/cka * log(lyx(i,j)/z0) + 8.0d0)
				!u(i,j,ki) = u(i,j,ki+1) - utau_ibm(i,j)*dabs(utau_ibm(i,j))*dz(ki)*rhox(i,j,ki+1)/mux(i,j,ki+1)
			endif
		enddo
		enddo
		!$omp end do nowait

		!$omp do private(ki,d_s,d_f,f,v_w)
		do j=1,ny1
		do i=1,nx
			ki = kv(i,j)
			if (ki > 0 .and. ki < nz) then
				d_s = max(obs_lsy(i,j,ki),    eps_min)
				d_f = max(-obs_lsy(i,j,ki+1), eps_min)
				f   = d_s/(d_s+d_f)
				f   = min(f, 0.9d0)                       ! limita a amplificação 1/(1-f) (estabilidade do IBM1)
				v_w = v(i,j,ki+1) - abs(vtau_ibm(i,j))*vtau_ibm(i,j)*dz(ki)*rhoy(i,j,ki) / &
				      max(muy(i,j,ki),eps_min)
				v(i,j,ki) = (v_w - f*v(i,j,ki+1))/(1.d0-f)
				!v(i,j,ki) = -v(i,j,ki+1)  +2.*0.043*0.043*dz(ki)/muy(i,j,ki)*rhoy(i,j,ki) !-v(i,j,ki+1) + 2.*0.043*(1.0d0/cka * log(lyy(i,j)/z0) + 8.0d0)
				!v(i,j,ki) = v(i,j,ki+1) - vtau_ibm(i,j)*dabs(vtau_ibm(i,j))*dz(ki)*rhoy(i,j,ki+1)/muy(i,j,ki+1)
			endif
		enddo
		enddo
		!$omp end do nowait

		! Copias das colunas de contorno (periodicidade em x/y).
		! ANTES: secoes de tamanhos diferentes -> atribuicao ilegal quando ku muda
		!$omp do private(k0)
		do j=0,ny+1
			k0 = ku(nx1,j)
			do k = 0, k0
				u(nx1+1,j,k) = u(nx1,j,k)
			enddo
			k0 = ku(1,j)
			do k = 0, k0
				u(0,j,k) = u(1,j,k)
			enddo
		enddo
		!$omp end do nowait

		!$omp do private(k0)
		do i=0,nx+1
			k0 = kv(i,ny1)
			do k = 0, k0
				v(i,ny1+1,k) = v(i,ny1,k)
			enddo
			k0 = kv(i,1)
			do k = 0, k0
				v(i,0,k) = v(i,1,k)
			enddo
		enddo
		!$omp end do nowait

		!$omp do
		do j=0,ny+1
		do i=0,nx+1
			w(i,j,0:kw(i,j))=0.
		enddo
		enddo
		!$omp end do
		!$omp end parallel

	elseif (ccz0 == 4 .and. ibm_t == 2) then

		! Lei da parede no IBM induzido: calcula-se tau_w pela lei de parede
		! usando nu molecular em y+, e reconstrói-se a velocidade alvo a partir
		! do gradiente correspondente usando mu efetiva, como no operador de
		! momento. Isso evita misturar a viscosidade molecular nesta conversão.
		!
		! OBS: a versão anterior aplicava v_sol = u_tau*(1/ka*ln(phi/z0)+8)
		! com u_tau fixo (0.043) e o MESMO valor positivo para u e v em todo
		! o domínio. Isso impõe um escorregamento espúrio +y em todas as
		! faces transversais do obstáculo, empurrando o escoamento para a
		! parede ccyf (y=yf), onde surge a oscilação indevida. A velocidade
		! de parede deve ser calculada componente a componente, com a
		! velocidade de atrito local (com sinal) obtida em lei_parede():
		!   u_w = u_1 - tau_w*dy/mu_eff
		!   v_w = v_1 - tau_w*dy/mu_eff
		!$omp parallel
		!$omp do private(ki)
		do j=1,ny
		do i=1,nx1
			ki = ku(i,j)
			if (ki > 0 .and. ki < nz) then !! REVER
				!xtau(i,j) = u(i,j,ki+1) - utau_ibm(i,j)*(1.0d0/cka *dlog(lyx(i,j)/z0) + 8.5d0)
				xtau(i,j) = u(i,j,ki+1) - abs(utau_ibm(i,j))*utau_ibm(i,j)*dz(ki)/mux(i,j,ki)*rhox(i,j,ki)
			else
				xtau(i,j) = 0.d0
			endif
		enddo
		enddo
		!$omp end do nowait

		!$omp do private(ki)
		do j=1,ny1
		do i=1,nx
			ki = kv(i,j)
			if (ki > 0 .and. ki < nz) then !! REVER
				!ytau(i,j) = v(i,j,ki+1) - vtau_ibm(i,j)*(1.0d0/cka*dlog(lyy(i,j)/z0) + 8.5d0
				ytau(i,j) = v(i,j,ki+1) - abs(vtau_ibm(i,j))*vtau_ibm(i,j)*dz(ki)/muy(i,j,ki)*rhoy(i,j,ki)
			else
				ytau(i,j) = 0.d0
			endif
		enddo
		enddo
		!$omp end do
		!$omp end parallel

		CALL ibm_images(obs_lsx,u,xtau,nx1,ny,nz,id_ibmx,1)
		CALL ibm_images(obs_lsy,v,ytau,nx,ny1,nz,id_ibmy,1)
		CALL ibm_images(obs_lsz,w,wtau,nx,ny,nz1,id_ibmz,1)

	endif


END SUBROUTINE semi_slip





!##################################################################################################################

SUBROUTINE lei_parede_ibm(lyx, lyy, utau, vtau, utau_ibm, vtau_ibm,rhox,rhoy,mux_mol,muy_mol,mux,muy)

USE velpre
USE obst
USE param
USE les
USE ls_param
!$ USE omp_lib

IMPLICIT NONE

integer :: i, j, k, iter, ki,kk
real(8) :: arg_log

real(8),dimension(nx1,ny) :: utau, utau_ibm, lyx
real(8),dimension(nx,ny1) :: vtau, vtau_ibm, lyy
real(8),dimension(nx1,ny,nz) :: mux, rhox
real(8),dimension(nx,ny1,nz) :: muy, rhoy
real(8),dimension(nx1,ny,nz) :: mux_mol
real(8),dimension(nx,ny1,nz) :: muy_mol

real(8), parameter :: eps_ln = 1.0d-10
real(8), parameter :: utau_ini = 1.0d-6

utau = utau_ini
vtau = utau_ini
utau_ibm = utau_ini
vtau_ibm = utau_ini
iter = 10

!agora fazemos os cálculos considerando a presença do IBM
if (ibm_t == 1 .or. ibm_t == 2) then
!aqui temos uma primeira estimativa para u* de 0.041. Caso isso dê menor que 5 com uma parede hidraulicamente lisa.
! aqui temos k = iter apenas para o procedimento ser repetido 20 vezes pois é um processo iterativo

	if(z0*0.041d0*rho_f2/mu_f2 < 5.d0) then

		!$omp parallel
		!$omp do private(ki,arg_log)
		do j=1,ny
		do i=1,nx1
			ki = ku(i,j)
			if (ki > 0 .and. ki < nz) then
				utau_ibm(i,j) = utau_ini
				do kk = 1, iter
					arg_log = lyx(i,j)*utau_ibm(i,j)*rhox(i,j,ki+1)/max(mux_mol(i,j,ki+1),eps_ln)
					utau_ibm(i,j) = dabs(u(i,j,ki+1)) / (1.0d0/cka * dlog(max(arg_log,eps_ln)) + 5.0d0)
				enddo
				utau_ibm(i,j) = dsign(1.0d0,u(i,j,ki+1)) * utau_ibm(i,j)
			endif
		enddo
		enddo
		!$omp end do nowait

		!$omp do private(ki,arg_log)
		do j=1,ny1
		do i=1,nx
			ki = kv(i,j)
			if (ki > 0 .and. ki < nz) then
				vtau_ibm(i,j) = utau_ini
				do kk = 1, iter
					arg_log = lyy(i,j)*vtau_ibm(i,j)*rhoy(i,j,ki+1)/max(muy_mol(i,j,ki+1),eps_ln)
					vtau_ibm(i,j) = dabs(v(i,j,ki+1)) / (1.0d0/cka * dlog(max(arg_log,eps_ln)) + 5.0d0)
				enddo
				vtau_ibm(i,j) = dsign(1.0d0,v(i,j,ki+1)) * vtau_ibm(i,j)
			endif
		enddo
		enddo
		!$omp end do
		!$omp end parallel
	else
		!$omp parallel
		!$omp do private(ki)
		do j=1,ny
		do i=1,nx1
			ki = ku(i,j)
			if (ki > 0 .and. ki < nz) then
				utau_ibm(i,j) = u(i,j,ki+1) / (1.0d0/cka * dlog(max(lyx(i,j)/max(z0,eps_ln),eps_ln))+8.5d0)
			endif
		enddo
		enddo
		!$omp end do nowait

		!$omp do private(ki)
		do j=1,ny1
		do i=1,nx
		ki = kv(i,j)
			if (ki > 0 .and. ki < nz) then
				vtau_ibm(i,j)= v(i,j,ki+1) / (1.0d0/cka * dlog(max(lyy(i,j)/max(z0,eps_ln),eps_ln))+8.5d0)
			endif
		enddo
		enddo
		!$omp end do
		!$omp end parallel
	endif
endif

END SUBROUTINE lei_parede_ibm


!##################################################################################################################

! Aplica semi-slip a todas as componentes tangenciais de uma parede.
! idir/lado identificam exclusivamente a orientação e o lado da parede;
! as componentes e as dimensões das respectivas malhas deslocadas são
! escolhidas aqui.
SUBROUTINE semi_slip_wall(idir, lado, hnorm, u, v, w, rhox, rhoy, rhoz, mux_mol, muy_mol, muz_mol, &
                          mux, muy, muz, &
                          nx, ny, nz, nx1, ny1, nz1)
	IMPLICIT NONE
	integer,intent(IN) :: idir,lado,nx,ny,nz,nx1,ny1,nz1
	real(8),intent(IN) :: hnorm
	real(8),intent(INOUT) :: u(0:nx1+1,0:ny+1,0:nz+1)
	real(8),intent(INOUT) :: v(0:nx+1,0:ny1+1,0:nz+1)
	real(8),intent(INOUT) :: w(0:nx+1,0:ny+1,0:nz1+1)
	real(8),intent(IN) :: rhox(nx1,ny,nz),mux(nx1,ny,nz)
	real(8),intent(IN) :: rhoy(nx,ny1,nz),muy(nx,ny1,nz)
	real(8),intent(IN) :: rhoz(nx,ny,nz1),muz(nx,ny,nz1)
	real(8),intent(IN) :: mux_mol(nx1,ny,nz),muy_mol(nx,ny1,nz),muz_mol(nx,ny,nz1)

	select case(idir)
	case(1) ! parede normal a x: v e w são tangenciais !!ainda não temos casos prontos para isso
		CALL semi_slip_face(v,rhoy,muy_mol,muy,nx,ny1,nz,1,lado,hnorm)
		CALL semi_slip_face(w,rhoz,muz_mol,muz,nx,ny,nz1,1,lado,hnorm)
	case(2) ! parede normal a y: u e w são tangenciais
		CALL semi_slip_face(u,rhox,mux_mol,mux,nx1,ny,nz,2,lado,hnorm)
		CALL semi_slip_face(w,rhoz,muz_mol,muz,nx,ny,nz1,2,lado,hnorm)
	case(3) ! parede normal a z: u e v são tangenciais
		CALL semi_slip_face(u,rhox,mux_mol,mux,nx1,ny,nz,3,lado,hnorm)
		CALL semi_slip_face(v,rhoy,muy_mol,muy,nx,ny1,nz,3,lado,hnorm)
	end select
END SUBROUTINE semi_slip_wall

!##################################################################################################################


! Aplica a condição semi-slip (lei da parede) a UMA componente de velocidade
! em uma direção genérica. A célula fantasma da parede é reconstruída de modo
! a impor a tensão de cisalhamento tau_w = rho*utau*|utau|, com utau obtido
! da lei logarítmica avaliada no primeiro nó fluido.
!
! idir : direção normal à parede (1=x, 2=y, 3=z)
! lado : 0 = parede inferior (fantasma no índice 0)
!        1 = parede superior (fantasma no índice dim+1)
! hnorm: espaçamento da malha na direção normal à parede (dx, dy ou dz(1))

SUBROUTINE semi_slip_face(a1, rho_a, mu_mol_a, mu_eff_a, dimx, dimy, dimz, idir, lado, hnorm)

	USE param, only : z0, cka
	USE ls_param, only : rho_f2, mu_f2
	!$ USE omp_lib

	IMPLICIT NONE

	integer,intent(IN) :: dimx, dimy, dimz, idir, lado
	real(8),intent(IN) :: hnorm
	real(8),intent(INOUT),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: a1
	real(8),intent(IN),dimension(dimx,dimy,dimz) :: rho_a, mu_mol_a, mu_eff_a

	integer :: i, j, k, n
	real(8) :: utau, tau_w, arg_log, lyw
	real(8), parameter :: utau_ini = 1.0d-6
	real(8), parameter :: eps_ln  = 1.0d-10
	integer, parameter :: iter = 10
	logical :: parede_lisa

	parede_lisa =  (1.d0 < 5.d0)!!(z0*0.041d0*rho_f2/mu_f2 < 5.d0) ! forçar a parede lisa já que o fundo vai ser por IBM e as laterais são lisas.


	! A lei logarítmica só é válida para y > z0, então a distância efetiva usada
	! nela é medida a partir da rugosidade (y + z0). Sem esse offset o argumento
	! do log fica abaixo de z0 e o denominador (1/ka*ln(y/z0) + 8.5) tende a zero
	! => u_tau e os nós fantasma explodem já na 1ª iteração (ccz0 = 4).
	lyw = 0.5d0*hnorm + z0  ! distância efetiva da parede até o 1º nó fluido

	!$omp parallel

	!====================
	! parede normal a X
	!====================
	if (idir == 1) then

		if (lado == 0) then
			!$omp do private(utau,tau_w,arg_log,n)
			do k = 1, dimz
			do j = 1, dimy
				utau = utau_ini
				if (parede_lisa) then
					do n = 1, iter
						arg_log = lyw*utau*rho_a(1,j,k)/max(mu_mol_a(1,j,k),eps_ln)
						utau = dabs(a1(1,j,k)) / (1.0d0/cka*dlog(max(arg_log,eps_ln)) + 5.0d0)
					enddo
					utau = dsign(1.0d0,a1(1,j,k))*utau
				else
					utau = a1(1,j,k) / (1.0d0/cka*dlog(max(lyw/max(z0,eps_ln),eps_ln)) + 8.5d0)
				endif
				tau_w = rho_a(1,j,k)*utau*dabs(utau)
				a1(0,j,k) = a1(1,j,k) - tau_w/max(mu_eff_a(1,j,k),eps_ln)*hnorm
			enddo
			enddo
			!$omp end do
		else
			!$omp do private(utau,tau_w,arg_log,n)
			do k = 1, dimz
			do j = 1, dimy
				utau = utau_ini
				if (parede_lisa) then
					do n = 1, iter
						arg_log = lyw*utau*rho_a(dimx,j,k)/max(mu_mol_a(dimx,j,k),eps_ln)
						utau = dabs(a1(dimx,j,k)) / (1.0d0/cka*dlog(max(arg_log,eps_ln)) + 5.0d0)
					enddo
					utau = dsign(1.0d0,a1(dimx,j,k))*utau
				else
					utau = a1(dimx,j,k) / (1.0d0/cka*dlog(max(lyw/max(z0,eps_ln),eps_ln)) + 8.5d0)
				endif
				tau_w = rho_a(dimx,j,k)*utau*dabs(utau)
				a1(dimx+1,j,k) = a1(dimx,j,k) + tau_w/max(mu_eff_a(dimx,j,k),eps_ln)*hnorm
			enddo
			enddo
			!$omp end do
		endif

	!====================
	! parede normal a Y
	!====================
	elseif (idir == 2) then

		if (lado == 0) then
			!$omp do private(utau,tau_w,arg_log,n)
			do k = 1, dimz
			do i = 1, dimx
				utau = utau_ini
				if (parede_lisa) then
					do n = 1, iter
						arg_log = lyw*utau*rho_a(i,1,k)/max(mu_mol_a(i,1,k),eps_ln)
						utau = dabs(a1(i,1,k)) / (1.0d0/cka*dlog(max(arg_log,eps_ln)) + 5.0d0)
					enddo
					utau = dsign(1.0d0,a1(i,1,k))*utau
				else
					utau = a1(i,1,k) / (1.0d0/cka*dlog(max(lyw/max(z0,eps_ln),eps_ln)) + 8.5d0)
				endif
				tau_w = rho_a(i,1,k)*utau*dabs(utau)
				a1(i,0,k) = a1(i,1,k) - tau_w/max(mu_eff_a(i,1,k),eps_ln)*hnorm
			enddo
			enddo
			!$omp end do
		else
			!$omp do private(utau,tau_w,arg_log,n)
			do k = 1, dimz
			do i = 1, dimx
				utau = utau_ini
				if (parede_lisa) then
					do n = 1, iter
						arg_log = lyw*utau*rho_a(i,dimy,k)/max(mu_mol_a(i,dimy,k),eps_ln)
						utau = dabs(a1(i,dimy,k)) / (1.0d0/cka*dlog(max(arg_log,eps_ln)) + 5.0d0)
					enddo
					utau = dsign(1.0d0,a1(i,dimy,k))*utau
				else
					utau = a1(i,dimy,k) / (1.0d0/cka*dlog(max(lyw/max(z0,eps_ln),eps_ln)) + 8.5d0)
				endif
				tau_w = rho_a(i,dimy,k)*utau*dabs(utau)
				a1(i,dimy+1,k) = a1(i,dimy,k) + tau_w/max(mu_eff_a(i,dimy,k),eps_ln)*hnorm
			enddo
			enddo
			!$omp end do
		endif

	!====================
	! parede normal a Z
	!====================
	elseif (idir == 3) then

		if (lado == 0) then
			!$omp do private(utau,tau_w,arg_log,n)
			do j = 1, dimy
			do i = 1, dimx
				utau = utau_ini
				if (parede_lisa) then
					do n = 1, iter
						arg_log = lyw*utau*rho_a(i,j,1)/max(mu_mol_a(i,j,1),eps_ln)
						utau = dabs(a1(i,j,1)) / (1.0d0/cka*dlog(max(arg_log,eps_ln)) + 5.0d0)
					enddo
					utau = dsign(1.0d0,a1(i,j,1))*utau
				else
					utau = a1(i,j,1) / (1.0d0/cka*dlog(max(lyw/max(z0,eps_ln),eps_ln)) + 8.5d0)
				endif
				tau_w = rho_a(i,j,1)*utau*dabs(utau)
				a1(i,j,0) = a1(i,j,1) - tau_w/max(mu_eff_a(i,j,1),eps_ln)*hnorm
			enddo
			enddo
			!$omp end do
		else
			!$omp do private(utau,tau_w,arg_log,n)
			do j = 1, dimy
			do i = 1, dimx
				utau = utau_ini
				if (parede_lisa) then
					do n = 1, iter
						arg_log = lyw*utau*rho_a(i,j,dimz)/max(mu_mol_a(i,j,dimz),eps_ln)
						utau = dabs(a1(i,j,dimz)) / (1.0d0/cka*dlog(max(arg_log,eps_ln)) + 5.0d0)
					enddo
					utau = dsign(1.0d0,a1(i,j,dimz))*utau
				else
					utau = a1(i,j,dimz) / (1.0d0/cka*dlog(max(lyw/max(z0,eps_ln),eps_ln)) + 8.5d0)
				endif
				tau_w = rho_a(i,j,dimz)*utau*dabs(utau)
				a1(i,j,dimz+1) = a1(i,j,dimz) + tau_w/max(mu_eff_a(i,j,dimz),eps_ln)*hnorm
			enddo
			enddo
			!$omp end do
		endif

	endif
	!$omp end parallel

END SUBROUTINE semi_slip_face

