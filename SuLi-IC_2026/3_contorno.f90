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
	real(8), save :: zi, zj, zk

	real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: dpdx
	real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: dpdy
	real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: dpdz
	real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: ls_v

	real(8),dimension(nx1,ny) :: utau2
	real(8),dimension(nx,ny1) :: vtau2
	real(8),dimension(nx,ny)  :: wtau2


	!RESOLUÇÃO DO PROBLEMA
	if ((nlock == 1).or.(nlock == 2)) then

		
		if (ibm_t == 1 .and. ccz0 .ne. 4) then
			!$omp parallel do
			do j=0,ny+1
			do i=0,nx+1
				u(i,j,0:ku(i,j))=0. !+ dpdx(i,j,0:ku(i,j))
				v(i,j,0:kv(i,j))=0. !+ dpdy(i,j,0:kv(i,j))
				w(i,j,0:kw(i,j))=0. !+ dpdz(i,j,0:kw(i,j))
			enddo
			enddo
			!$omp end parallel do
					
			i = nx1+1 !j=todos
			!$omp parallel do
			do j=0,ny+1
				u(i,j,0:ku(i,j))=0. !+ dpdx(i,j,0:ku(i,j))
			enddo
			!$omp end parallel do

			j=ny1+1 !i=todos
			!$omp parallel do
			do i=0,nx+1
				v(i,j,0:kv(i,j))=0. !+ dpdy(i,j,0:kv(i,j))
			enddo
			!$omp end parallel do
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
		
		CALL cc_geral(prd1,nx,ny,nz)

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

		!Prescrita
		elseif (ccyf.eq.3) then 
			u(:,ny+1,:)  = byxf(:,:)
			v(:,ny1+1,:) = byyf1(:,:)
			w(:,ny+1,:)  = byzf(:,:)
		endif


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

		!Semi-slip condition (Ferziger et al., 2020; Chow et al., 2005)
		!Rugosidade sendo aplicada apenas no fundo do canal
		elseif (ccz0.eq.4) then									
			CALL semi_slip()
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

		if (ccx0 .eq. 3 .or. ccxf .eq. 3 .or. ccy0 .eq. 3 .or. ccyf .eq. 3 .or. ccz0 .eq. 3 .or. cczf .eq. 3 .or. ccz0 .eq. 4) then
			CALL prd_corr(dpdx,dpdy,dpdz)
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

		!Prescrita
		elseif (ccyf.eq.3) then 
			v(:,ny1,:)   = byyf(:,:) + dpdy(:,ny1,:)
		
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

			ls(1,:,:)  = ls(2,:,:)
			ls(nx,:,:) = ls(nx-1,:,:)

			ls(:,1,:)  = ls(:,2,:)
			ls(:,ny,:) = ls(:,ny-1,:)			
			
			
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
	
	!Free-slip
	if (ccyf.eq.1) then
		ka(:,ny+1,:) = ka(:,ny,:)
	endif

	!No-slip condition
	if (ccy0.eq.2) then
		ka(:,0,:) = -ka(:,1,:)
	endif
	
	if (ccyf.eq.2) then
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

	!EIXO Z
	
	!Free-slip condition
	if (ccz0.eq.1) then
		ka(:,:,0) = ka(:,:,1)
	endif
	
	!Free-slip
	if (cczf.eq.1) then
		ka(:,:,nz+1) = ka(:,:,nz)
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
	return

	elseif (obst_t == 1) then !Rugosidade uniforme experimento GS_H128 (Zampiron et al., 2022)

        !$omp parallel do
		do k = -1,nzz+2
		do j = -1,nyy+2
		do i = -1,nxx+2
			!auxx_ls(i,j,k)  = zz(5) + 0.008*(1.+(dabs(dsin(pi*xx(i)/0.016))*dabs(dsin(pi*yy(j)/0.016)))**0.5)-zz(k)
			auxx_ls(i,j,k)  = zz(5) -zz(k)
			
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


        !$omp parallel do
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
	

	
		!$omp parallel do
		do k = -1, nzz+2
		do j = -1, nyy+2
		do i = -1, nxx+2
			do mm = 2,n 
				if (yy(j) < ponto(mm,1) .and. yy(j) >= ponto(mm-1,1)) then    
					auxx_ls(i,j,k) = ponto(mm,2) + (ponto(mm,2)-ponto(mm-1,2))/(ponto(mm,1)-ponto(mm-1,1))*(yy(j)-ponto(mm,1))-zz(k) 
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

   	if (ibm_t == 1 .or. ccz0.eq.4) then
   	! EM TESTE ! IBM COM DZ VARIÁVEL

   	ku = 0
   	kv = 0
   	kw = 0


!!        !$omp parallel do

   	do k = 1, nz1
	do j = 0, ny1
	do i = 0, nx1+1
		if (	(auxx_ls(i*2-1,j*2,k*2)*auxx_ls(i*2-1,j*2,k*2-1)) < 0.) then
			ku(i,j) = k-1
		elseif ((auxx_ls(i*2-1,j*2,k*2)*auxx_ls(i*2-1,j*2,k*2+1)) <= 0.) then
			ku(i,j) = k
		endif
	enddo
	enddo
	
	
	do j = 0, ny1+1
	do i = 0, nx1
		if (	(auxx_ls(i*2,j*2-1,k*2)*auxx_ls(i*2,j*2-1,k*2-1)) < 0.) then
			kv(i,j) = k-1
		elseif ((auxx_ls(i*2,j*2-1,k*2)*auxx_ls(i*2,j*2-1,k*2+1)) <= 0.) then
			kv(i,j) = k
		endif
	enddo
	enddo
	
	do j = 0, ny1
	do i = 0, nx1
		if (	(auxx_ls(i*2,j*2,k*2-1)*auxx_ls(i*2,j*2,k*2-1-1)) < 0.) then
			kw(i,j) = k-1
		elseif ((auxx_ls(i*2,j*2,k*2-1)*auxx_ls(i*2,j*2,k*2-1+1)) <= 0.) then
			kw(i,j) = k
		endif
	enddo
	enddo
   	enddo
!!        !$omp end parallel do

	endif
	
	!if (ibm_t == 1 ) then
	!	auxx = -10.
		! converter a visão nova para a antiga (para manter o ibm antigo vivo)
	!	do k = -1, nzz+1 !reduziu um por ser avaliado de forma progressiva
	!	do j = -1, nyy+2
	!	do i = -1, nxx+2
	!		if ((auxx_ls(i,j,k)*auxx_ls(i,j,k+1)) <= 0.) then !condição para definir a posição da interface
	!			auxx(i,j) = zz(k)
	!		endif
	!	enddo
	!	enddo
	!	enddo


	!	do k = 0, ny+1
	!	do j = 0, ny+1
	!	do i = 0, nx1+1
	!		if ((auxx_ls(i*2,j*2,k*2)*auxx_ls(i*2,j*2,k*2+1)) <= 0.) then !condição para definir a posição da interface
	!			auxx(i,j) = zz(k)
	!		endif
	!	enddo
	!	enddo
	!	enddo

	!	do j = 0, ny+1
	!	do i = 0, nx1+1
	!	 	ku(i,j) = max(0,nint(auxx(i*2-1,j*2)/dz+0.5)) 
	!		if (ku(i,j) > nz) ku(i,j) = nz
	!	enddo
	!	enddo

	!	do j = 0, ny1+1
	!	do i = 0, nx+1
	!		kv(i,j) = max(0,nint(auxx(i*2,j*2-1)/dz+0.5))
	!		if (kv(i,j) > nz) kv(i,j) = nz
	!	enddo
	!	enddo

	!	do j = 0, ny+1
	!	do i = 0, nx+1
	!		kw(i,j) = max(0,nint(auxx(i*2,j*2)/dz+1.))
	!		if (kw(i,j) > nz1) kw(i,j) = nz1
	!	enddo
	!	enddo


	if (ibm_t == 2) then

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

		do j = 1, ny
		do i = nint(x_inicial/dx)+1, nx
		do k = 1, nz-1
		
		if (ls(i,j,k)*ls(i,j,k+1) <= 0) then
			alt1 = z(k) + ls(i,j,k)
			alt2 = z(k+1) + ls(i,j,k+1)
			zpfs(i) = max(alt1,alt2) ! altura da interface
		endif
		enddo
		enddo
		enddo

	!Camada esponja para a direção x

		do k = 1, nz
		do j = 1, ny
		do i = nint(x_inicial/dx)+1, nx
		
		xp(i) = real(i-0.5)*dx
		aux1 = (xp(i)-x_inicial)/l_sponge

                if (ls(i,j,k) >= 0.) then !Na água 
                        epis_z(i,j,k) = alfa * aux1 * aux1 * (0.-z(k)) / (0.-zpfs(i))
                else
                        epis_z(i,j,k) = alfa * aux1 * aux1 * (z(k)-lz) / (zpfs(i)-lz)
                endif

		enddo
		enddo
		enddo

	elseif (esp_type == 2) then


        !$omp parallel do
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


	CALL cc_geral(dpdx,nx1,ny,nz)
	CALL cc_geral(dpdy,nx,ny1,nz)
	CALL cc_geral(dpdz,nx,ny,nz1)


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

SUBROUTINE semi_slip()

	USE velpre
	USE obst
	USE param
	!$ USE omp_lib

	IMPLICIT NONE

	!Declarado também no programa
	integer :: i, j, k
	real(8), save :: cd
	
	real(8),dimension(nx1,ny,nz) :: vx
	real(8),dimension(nx,ny1,nz) :: uy
		
	real(8),dimension(nx1,ny) :: xtau, utau
	real(8),dimension(nx,ny1) :: ytau, vtau
	real(8),dimension(nx,ny)  :: ztau, wtau

	utau = 0.
	vtau = 0.
	wtau = 0.

	!Interpolação velocidades	
	CALL interpxy_ff(u,nx,ny,nz,nx1,ny1,uy)
	CALL interpyx_ff(v,nx,ny,nz,nx1,ny1,vx)

		
	!Coeficiente de arrasto
	cd = (1./cka * dlog((dzmin*0.5 + z0) / z0))**(-2.) 
						 							
	!$omp parallel
	!$omp do
	do j=1,ny
	do i=1,nx1
		if (ku(i,j) == 0) then
			!Cisalhamento de fundo (Ferziger et al., 2020)
			xtau(i,j) = cd * u(i,j,1) * dsqrt( u(i,j,1)*u(i,j,1) + vx(i,j,1)*vx(i,j,1) )
			
			!Velocidade de cisalhamento
			utau(i,j) = dsign(dsqrt(dabs(xtau(i,j))), xtau(i,j))
			u(i,j,0) = 2.*utau(i,j) - u(i,j,1)
		endif
	enddo
	enddo
	!$omp end do nowait
	
	!$omp do
	do j=1,ny1
	do i=1,nx
		if (kv(i,j) == 0) then
			!Cisalhamento de fundo (Ferziger et al., 2020)
			ytau(i,j) = cd * v(i,j,1) * dsqrt( uy(i,j,1)*uy(i,j,1) + v(i,j,1)*v(i,j,1) )

			!Velocidade de cisalhamento
			vtau(i,j) = dsign(dsqrt(dabs(ytau(i,j))), ytau(i,j))
			v(i,j,0) = 2.*vtau(i,j) - v(i,j,1)
		endif
	enddo
	enddo
	!$omp end do nowait

	!$omp do
	do j=0,ny+1
	do i=0,nx+1
		if (kw(i,j) == 0) then
			w(i,j,0) = w(i,j,2)
		endif
	enddo
	enddo
	!$omp end do
	!$omp end parallel


									 							
	!Obstáculo (ibm forçado)
	if (ibm_t == 1 .or. ibm_t == 2) then
		!!$omp parallel 
		!!$omp do 
		do j=1,ny
		do i=1,nx1
		
			if (ku(i,j) > 0) then
				!Cisalhamento de fundo (Ferziger et al., 2020)
				xtau(i,j) = cd * u(i,j,ku(i,j)) * dsqrt( u(i,j,ku(i,j))*u(i,j,ku(i,j)) + vx(i,j,ku(i,j))*vx(i,j,ku(i,j)) )

				!Velocidade de cisalhamento
				utau(i,j) = dsign(dsqrt(dabs(xtau(i,j))), xtau(i,j))
			endif
		enddo
		enddo
		!!$omp end do nowait

		!!$omp do
		do j=1,ny1
		do i=1,nx
			if (kv(i,j) > 0) then
				!Cisalhamento de fundo (Ferziger et al., 2020)
				ytau(i,j) = cd * v(i,j,kv(i,j)) * dsqrt( uy(i,j,kv(i,j))*uy(i,j,kv(i,j)) + v(i,j,kv(i,j))*v(i,j,kv(i,j)) )

				!Velocidade de cisalhamento
				vtau(i,j) = dsign(dsqrt(dabs(ytau(i,j))), ytau(i,j))
		
			endif
		enddo
		enddo
		!!$omp end do
		!!$omp end parallel

		wtau = 0.
	endif
				
	! será que não devemos impor esta velocidade em vez de calcular como ali em cima?
	!utau = 0.05516
	!vtau = 0.
	!wtau = 0.

	if (ibm_t == 1) then	
		do j=1,ny
		do i=1,nx1
			if (ku(i,j) > 0) then
				u(i,j,ku(i,j)-1) = 2.*utau(i,j) - u(i,j,ku(i,j))
				u(i,j,0:ku(i,j)-2)=0. !+ dpdx(i,j,0:ku(i,j))
			endif
		enddo
		enddo
				
		do j=1,ny1
		do i=1,nx
			if (kv(i,j) > 0) then
				v(i,j,kv(i,j)-1) = 2.*vtau(i,j) - v(i,j,kv(i,j))
				v(i,j,0:kv(i,j)-2)=0. !+ dpdy(i,j,0:kv(i,j))		
			endif
		enddo
		enddo				

		i = nx1+1 !j=todos
		do j=0,ny+1
			u(i,j,0:ku(i,j))=u(i-1,j,0:ku(i-1,j)) !+ dpdx(i,j,0:ku(i,j))
		enddo
		i = 0 !j=todos
		do j=0,ny+1
			u(i,j,0:ku(i,j))=u(i+1,j,0:ku(i+1,j)) !+ dpdx(i,j,0:ku(i,j))
		enddo

		j=ny1+1 !i=todos
		do i=0,nx+1
			v(i,j,0:kv(i,j))= v(i,j-1,0:kv(i,j-1)) !+ dpdy(i,j,0:kv(i,j))
		enddo
		j=0 !i=todos
		do i=0,nx+1
			v(i,j,0:kv(i,j))= v(i,j+1,0:kv(i,j+1)) !+ dpdy(i,j,0:kv(i,j))
		enddo
		
		!$omp parallel do
		do j=0,ny+1
		do i=0,nx+1
			w(i,j,0:kw(i,j))=0.   !+ dpdz(i,j,0:kw(i,j))
		enddo
		enddo
		!$omp end parallel do
				
	elseif (ibm_t == 2) then
					
		!cria as imagens interpoladas
		CALL ibm_images(obs_lsx,u,utau,nx1,ny,nz,id_ibmx,1)
		CALL ibm_images(obs_lsy,v,vtau,nx,ny1,nz,id_ibmy,1)	
		CALL ibm_images(obs_lsz,w,wtau,nx,ny,nz1,id_ibmz,1)
	endif



END SUBROUTINE semi_slip



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

				umed = u(i,j,k) + umed
				ik   = 1 + ik
				iik  = k
			elseif (ls_x(i,j,k)<0) then
				u(i,j,k) = u(i,j,iik)*(1.+tanh(pi*(iik-k-1)/2.))
			endif
		enddo
		enddo


		if (i == 0) then
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

SUBROUTINE cc_geral(a,dimx,dimy,dimz)

	USE cond
	!$ USE omp_lib

	IMPLICIT NONE

	!Declarado também no programa
	integer :: i, j, k, dimx, dimy, dimz
	real(8),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: a
			

	!$omp parallel
	if (ccx0.eq.0) then  ! Condição periódica
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
			a(0,j,k)    = a(1,j,k)
			a(dimx+1,j,k) = a(dimx,j,k)
		enddo
		enddo
		!$omp end do
	endif

	if (ccy0.eq.0) then  ! Condição periódica
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





