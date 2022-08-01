!Subrotina definir as condições de contorno das velocidades e suas influências
!Referência: Gotoh, 2013

!Implementação em 15/04/2014
!Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 18/08/2016

SUBROUTINE contorno(nlock)

USE velpre
USE obst
   USE cond
USE ls_param

IMPLICIT NONE

!Declarado também no programa
real(8) :: zi, zj, zk
integer :: i, j, k, niv, ii,nlock

real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: dpdx
real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: dpdy
real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: dpdz

!RESOLUÇÃO DO PROBLEMA

if ((nlock == 1).or.(nlock == 2)) then
	if (mms_t > 0) call mms_bc()
		!Velocidades de atrito
		!ub = 0.
		!vb = 0.
		!pressão
		!grad 0
				   
		if (ccx0.eq.0.and.ccxf.eq.0) then
		   prd1(0,:,:)    = prd1(nx,:,:)
		   prd1(nx+1,:,:) = prd1(1,:,:)
		else
		   prd1(0,:,:)    = prd1(1,:,:)
		   prd1(nx+1,:,:) = prd1(nx,:,:)
		endif

		if (ccy0.eq.0.and.ccyf.eq.0) then
		   prd1(:,0,:)    = prd1(:,ny,:)
		   prd1(:,ny+1,:) = prd1(:,1,:)
		else
		   prd1(:,0,:)    = prd1(:,1,:)
		   prd1(:,ny+1,:) = prd1(:,ny,:)
		endif

		prd1(:,:,0)    = prd1(:,:,1)
		prd1(:,:,nz+1) = prd1(:,:,nz)

		!Parede esquerda (j = 1)
		!Periodica
		if (ccy0.eq.0.and.ccyf.eq.0) then
			u(:,0,:) = u(:,ny,:)    
			v(:,0,:) = v(:,ny1-1,:)
			w(:,0,:) = w(:,ny,:)    
		endif
		
		!Free-slip condition
		if (ccy0.eq.1) then 
			u(:,0,:) = u(:,1,:)
			v(:,0,:) = -v(:,2,:)
			w(:,0,:) = w(:,1,:)
		endif
		!No-slip condition
		if (ccy0.eq.2) then 
			u(:,0,:) = -u(:,1,:)
			v(:,0,:) = v(:,2,:)
			w(:,0,:) = -w(:,1,:)
		endif
		
		!Prescrita
		if (ccy0.eq.3) then
			u(:,0,:) = byx0(:,:)
			v(:,0,:) = byy0(:,:)
			w(:,0,:) = byz0(:,:)
		endif

		!Parede direita (j = ny ou ny1)
		!Periodica
		if (ccy0.eq.0.and.ccyf.eq.0) then
			u(:,ny+1,:)  = u(:,1,:)
			v(:,ny1+1,:) = v(:,2,:)
			w(:,ny+1,:)  = w(:,1,:)
		endif

		!Free-slip condition
		if (ccyf.eq.1) then
			u(:,ny+1,:)  = u(:,ny,:)
			v(:,ny1+1,:) = -v(:,ny1-1,:)
			w(:,ny+1,:)  = w(:,ny,:)
		endif
		
		!No-slip condition
		if (ccyf.eq.2) then 
			u(:,ny+1,:) = -u(:,ny,:)
			v(:,ny1+1,:) = v(:,ny1-1,:)
			w(:,ny+1,:) = -w(:,ny,:)
		endif

		!Prescrita
		if (ccyf.eq.3) then 
			u(:,ny+1,:)  = byxf(:,:)
			v(:,ny1+1,:) = byyf1(:,:)
			w(:,ny+1,:)  = byzf(:,:)
		endif

		!Parede frente (i = 1)
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(0,:,:) = u(nx1-1,:,:)
			v(0,:,:) = v(nx,:,:)
			w(0,:,:) = w(nx,:,:)
		endif
		
		!Free-slip condition
		if (ccx0.eq.1) then
			u(0,:,:) = -u(2,:,:)
			v(0,:,:) = v(1,:,:)
			w(0,:,:) = w(1,:,:)
		endif

		!No-slip condition
		if (ccx0.eq.2) then
			u(0,:,:) = u(2,:,:)
			v(0,:,:) = -v(1,:,:)
			w(0,:,:) = -w(1,:,:)
		endif

		!Prescrita
		if (ccx0.eq.3) then
			u(0,:,:) = bxx0(:,:)
			v(0,:,:) = bxy0(:,:)
			w(0,:,:) = bxz0(:,:)
		endif
		
		!Para validação
		if (ccx0.eq.4) then 
			v(0,:,:) = bxy0(:,:)
			w(0,:,:) = w(1,:,:)
			!do i = 0, 1
				!u(i,:,11+1) =	0.15946 + dpdx(i,:,12)
				!u(i,:,12+1) =	0.2873  + dpdx(i,:,13)
				!u(i,:,13+1) =	0.328   + dpdx(i,:,14)
				!u(i,:,14+1) =	0.36376 + dpdx(i,:,15)
				!u(i,:,15+1) =	0.39226 + dpdx(i,:,16)
				!u(i,:,16+1) =	0.41742 + dpdx(i,:,17)
				!u(i,:,17+1) =	0.44166 + dpdx(i,:,18)
				!u(i,:,18+1) =	0.46318 + dpdx(i,:,19)
				!u(i,:,19+1) =	0.48141 + dpdx(i,:,20)
				!u(i,:,20+1) =	0.4867  + dpdx(i,:,21)
			!enddo
		endif

		!Parede de trás (i = nx ou nx1)	
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(nx1+1,:,:) = u(2,:,:)
			v(nx+1,:,:)  = v(1,:,:)
			w(nx+1,:,:)  = w(1,:,:)
		endif

		!Free-slip condition
		if (ccxf.eq.1) then
			u(nx1+1,:,:) = -u(nx1-1,:,:)
			v(nx+1,:,:)  = v(nx,:,:)
			w(nx+1,:,:)  = w(nx,:,:)
		endif

		!No-slip condition
		if (ccxf.eq.2) then
			u(nx1+1,:,:) = u(nx1-1,:,:)
			v(nx+1,:,:)  = -v(nx,:,:)
			w(nx+1,:,:)  = -w(nx,:,:)
		endif

		!Prescrita
		if (ccxf.eq.3) then
			u(nx1+1,:,:) = bxxf1(:,:)
			v(nx+1,:,:)  = bxyf(:,:)
			w(nx+1,:,:)  = bxzf(:,:)
		endif

		!Saida livre
		if (ccxf.eq.4) then
			u(nx1+1,:,:) = u(nx1,:,:)
			v(nx+1,:,:)  = v(nx,:,:)
			w(nx+1,:,:)  = w(nx,:,:)
		endif

		!Parede do fundo (k = 1)
		!Free-slip
		if (ccz0.eq.1) then
			u(:,:,0) = u(:,:,1)
			v(:,:,0) = v(:,:,1)
			w(:,:,0) = -w(:,:,2)
		endif

		!No-slip condition
		if (ccz0.eq.2) then
			u(:,:,0) = -u(:,:,1)
			v(:,:,0) = -v(:,:,1)
			w(:,:,0) = w(:,:,2)
		endif

		!Velocidade prescrita para manufaturada
		if (ccz0.eq.3) then
			u(:,:,0) = bzx0(:,:)
			v(:,:,0) = bzy0(:,:)
			w(:,:,0) = bzz0(:,:)
		endif

		!Superfície Livre (k = nz ou nz1)
		!Sempre na superfície livre será free-slip(em teste)
		if (cczf.eq.1) then
			u(:,:,nz+1)  = u(:,:,nz)
			v(:,:,nz+1)  = v(:,:,nz)
			w(:,:,nz1+1) = -w(:,:,nz1-1)
		endif
		
		!Sempre na superfície livre será no-slip(em teste)
		if (cczf.eq.2) then
			u(:,:,nz+1)  = -u(:,:,nz)
			v(:,:,nz+1)  = -v(:,:,nz)
			w(:,:,nz1+1) = w(:,:,nz1-1)
		endif

		!Velocidade prescrita para manufaturada
		if (cczf.eq.3) then
			u(:,:,nz+1)  = bzxf(:,:)
			v(:,:,nz+1)  = bzyf(:,:)
			w(:,:,nz1+1) = bzzf1(:,:)
		endif
	
		!Obstáculo de fundo
		!Ativar se tiver osbstáculo de fundo
		!ku, kv, kw indicam até que altura as velocidades tem que ser zeradas (até qual índice k)

		if (obst_t == 0) then
			if (ccz0.eq.2) then
				ub(:,:,2) = u(:,:,2)
				vb(:,:,2) = v(:,:,2)
				wb(:,:,2) = w(:,:,2)
			endif
				
			if (ccy0.eq.2) then
				ub(:,2,:) = u(:,2,:)
				vb(:,2,:) = v(:,2,:)
				wb(:,2,:) = w(:,2,:)
			endif
				
			if (ccyf.eq.2) then
				ub(:,ny1-1,:) = u(:,ny1-1,:)
				vb(:,ny1-1,:) = v(:,ny1-1,:)
				wb(:,ny1-1,:) = w(:,ny1-1,:)
			endif	
			
		else
			do j=0,ny+1
			do i=0,nx+1
				u(i,j,0:ku(i,j)-1)=0. !+ dpdx(i,j,0:ku(i,j))
				v(i,j,0:kv(i,j)-1)=0. !+ dpdy(i,j,0:kv(i,j))
				w(i,j,0:kw(i,j)-1)=0. !+ dpdz(i,j,0:kw(i,j))
				!Rugosidade Interna
				ub(i,j,ku(i,j)) = u(i,j,ku(i,j)+1)
				vb(i,j,kv(i,j)) = v(i,j,kv(i,j)+1)
				wb(i,j,kw(i,j)) = w(i,j,kw(i,j)+1)
				!ub(i,j,ku(i,j)+1) = 0.
				!vb(i,j,kv(i,j)+1) = 0.
				!wb(i,j,kw(i,j)+1) = 0.
				 
			enddo
			enddo
			
			i = nx1+1 !j=todos
				do j=0,ny+1
					u(i,j,0:ku(i,j))=0. !+ dpdx(i,j,0:ku(i,j))
				enddo

			j=ny1+1 !i=todos
				do i=0,nx+1
					v(i,j,0:kv(i,j))=0. !+ dpdy(i,j,0:kv(i,j))
				enddo
		endif


endif

if (nlock == 2) then

	if (ccx0 .eq. 3 .or. ccxf .eq. 3 .or. ccy0 .eq. 3 .or. ccyf .eq. 3 .or. ccz0 .eq. 3 .or. cczf .eq. 3) then
		CALL prd_corr(dpdx,dpdy,dpdz)
	endif

	if (ccy0.eq.0.and.ccyf.eq.0) then
		v(:,1,:) = v(:,ny1,:)     
	endif
		
		!Free-slip condition
		if (ccy0.eq.1) then 
			v(:,1,:) = 0.
		endif
		!No-slip condition
		if (ccy0.eq.2) then 
			v(:,1,:) = 0.
		endif
		
		!Prescrita
		if (ccy0.eq.3) then
			v(:,1,:) = byy1(:,:) + dpdy(:,1,:)
		endif

		!Parede direita (j = ny ou ny1)
		!Periodica
		if (ccy0.eq.0.and.ccyf.eq.0) then
			v(:,ny1,:)   = v(:,1,:)
		endif

		!Free-slip condition
		if (ccyf.eq.1) then
			v(:,ny1,:)   = 0.
		endif
		
		!No-slip condition
		if (ccyf.eq.2) then 
			v(:,ny1,:) = 0.
		endif

		!Prescrita
		if (ccyf.eq.3) then 
			v(:,ny1,:)   = byyf(:,:) + dpdy(:,ny1,:)
		endif

		!Parede frente (i = 1)
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(1,:,:) = u(nx1,:,:)
		endif
		
		!Free-slip condition
		if (ccx0.eq.1) then
			u(1,:,:) = 0.
		endif

		!No-slip condition
		if (ccx0.eq.2) then
			u(1,:,:) = 0.
		endif

		!Prescrita
		if (ccx0.eq.3) then
			u(1,:,:) = bxx1(:,:) + dpdx(1,:,:)
			do k = 1, nz
			ls(1,:,k) = blx1(1:ny,k)
			enddo 
		endif
		
		!Para validação
		if (ccx0.eq.4) then 
			!do i = 0, 1
				!u(i,:,11+1) =	0.15946 + dpdx(i,:,12)
				!u(i,:,12+1) =	0.2873  + dpdx(i,:,13)
				!u(i,:,13+1) =	0.328   + dpdx(i,:,14)
				!u(i,:,14+1) =	0.36376 + dpdx(i,:,15)
				!u(i,:,15+1) =	0.39226 + dpdx(i,:,16)
				!u(i,:,16+1) =	0.41742 + dpdx(i,:,17)
				!u(i,:,17+1) =	0.44166 + dpdx(i,:,18)
				!u(i,:,18+1) =	0.46318 + dpdx(i,:,19)
				!u(i,:,19+1) =	0.48141 + dpdx(i,:,20)
				!u(i,:,20+1) =	0.4867  + dpdx(i,:,21)
			!enddo
		endif

		!Parede de trás (i = nx ou nx1)	
		!Periodica
		if (ccx0.eq.0.and.ccxf.eq.0) then
			u(nx1,:,:)   = u(1,:,:)
		endif

		!Free-slip condition
		if (ccxf.eq.1) then
			u(nx1,:,:)   = 0.
		endif

		!No-slip condition
		if (ccxf.eq.2) then
			u(nx1,:,:)   = 0.
		endif

		!Prescrita
		if (ccxf.eq.3) then
			u(nx1,:,:)   = bxxf(:,:)   + dpdx(nx1,:,:)
		endif

		!Saida livre
		if (ccxf.eq.4) then
			u(nx1,:,:)   = u(nx1-1,:,:)
		endif

		!Parede do fundo (k = 1)
		!Free-slip
		if (ccz0.eq.1) then
			w(:,:,1) = 0.
		endif

		!No-slip condition
		if (ccz0.eq.2) then
			w(:,:,1) = 0.
		endif

		!Velocidade prescrita para manufaturada
		if (ccz0.eq.3) then
			w(:,:,1) = bzz1(:,:) + dpdz(:,:,1)
		endif

		!Superfície Livre (k = nz ou nz1)
		!Sempre na superfície livre será free-slip(em teste)
		if (cczf.eq.1) then
			w(:,:,nz1)   = 0.
		endif
		
		!Sempre na superfície livre será no-slip(em teste)
		if (cczf.eq.2) then
			w(:,:,nz1)   = 0.
		endif

		!Velocidade prescrita para manufaturada
		if (cczf.eq.3) then
			w(:,:,nz1)   = bzzf(:,:)    + dpdz(:,:,nz1)
		endif
	


elseif ((nlock == 3) .and. (obst_t .ne. 0)) then
		!Contorno para Level Set, o ideal é que o objeto seja representado por pelo menos dois grid por direção
		do j = 1,ny
		do i = 1,nx
		do k = kw(i,j)-1, 1 

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
		enddo
		enddo
		enddo

		CALL heaviside()
endif




END SUBROUTINE contorno

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

IMPLICIT NONE

real(8) :: x, y, x0, y0, a, d, sigx, sigy, tgaux, aux1, aux2, dz_v, erro, raio		!x mostra a localização atual
real(8),dimension(-1:nx1*2+2,-1:ny1*2+2) :: auxx
integer :: i,j, nxx, nyy

nxx = nx1*2
nyy = ny1*2

!ku, kv, kw indicam até que altura as velocidades tem que ser zeradas (até qual índice k)

!Desenhar um único objeto e fazer a divisão na hora de transformar em inteiro

if (obst_t == 0) then !Dunas (Daniel, Leonardo)
	write(*,*) "Sem obstáculo."
return

elseif (obst_t == 1) then !Dunas (Daniel, Leonardo)
	do j = -1,nyy+2
	do i = -1, nxx+2
		x = (i-1.)*(dx*0.5)
		auxx(i,j) = elev + amp*(1.+sin(fase+2.*pi*x/comp)) !Número de velocidades U zeradas, dependente da funçao altura da duna na parede esquerda do diferencial de volume
		!auxy(i,j) = elev + amp*(1.+sin(fase+2.*pi*x1/comp))-0.5*dz1 !Número de velocidades V zeradas
		!auxz(i,j) = elev + amp*(1.+sin(fase+2.*pi*x1/comp)) !Número de velocidades W zeradas
	enddo	
	enddo	

elseif (obst_t == 2) then !Duna de Yue et al. (2003) (Luisa, Leonardo)

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1.)*(dx*0.5)

	if (x.le.0.01) then
		auxx(i,j) = 0.02

	elseif (x.le.0.05) then
		tgaux=-tan(26.*pi/180.)
		auxx(i,j)=tgaux*x+0.02

	elseif (x.le.0.075) then
		auxx(i,j)=0.0005

	elseif (x.le.0.1375) then
		tgaux=tan(1.8*pi/180.)
		auxx(i,j)=tgaux*x+0.0005
		aux1=auxx(i,j)

	elseif (x.le.0.325) then
		tgaux=tan(5.*pi/180.)
		auxx(i,j)=tgaux*x+aux1
		aux2=auxx(i,j)

	elseif (x.le.0.39) then
	   tgaux=tan(1.8*pi/180.)
	   auxx(i,j)=tgaux*x+aux2
	else
		auxx(i,j)=0.02

	endif
	enddo
	enddo
	
elseif (obst_t == 3) then !gaussbump3D (Luisa, Leonardo)

	d=0.2    !4.!0.1/dz1+1. !Base parte "plana"

	x0 =3.       !3./dx!2.*nx/3. !Centro em x, do pico
	y0 =ny/2.*dy !ny/2.          !Centro em y, do pico
	a  = 2.0     !2./dz1!nz/5.   !Altura do pico
	sigx = 0.5   !0.5/dx!nx/5.   !Largura em x do pico
	sigy = 0.5   !0.5/dy!ny/5.   !Largura em y do pico

	do j=-1,nyy+2
		y = (j-1.)*(dy*0.5)
		aux1=(y-y0)*(y-y0)/(2.*sigy*sigy) 

	do i=-1,nxx+2
		x = (i-1.)*(dx*0.5)
		aux2=(x-x0)*(x-x0)/(2.*sigx*sigx)
		aux2=-aux2-aux1
		aux2=exp(aux2)
		auxx(i,j)=aux2*a +d
	enddo
	enddo

elseif (obst_t == 4) then !Beji e Battjes,(1994) (Leonardo)

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1)*(dx*0.5)

	if (x < 6.) then !Antes de subir
		auxx(i,j) = 0.

	elseif (x < 12.) then !Final da subida
		tgaux=tan(1./20.)
		auxx(i,j)=tgaux*(x-6.)

	elseif (x < 14.) then !Em cima do obstáculo
		auxx(i,j)= 0.30

	elseif (x < 17.) then !Descida
		tgaux=-tan(1./10.)
		auxx(i,j)=tgaux*(x-14.)+0.30
		
	else !Resto
		auxx(i,j) = 0.
	endif
	
	enddo
	enddo

elseif (obst_t == 5) then !(canal 0 - delft 1980)

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-2)*(dx*0.5)
	
	if (x <= 1.0) then !Raso plano
		auxx(i,j) = 0.2
		
	else !Ffundo plano
		auxx(i,j) = 0.
	endif
	enddo
	enddo

elseif (obst_t == 6) then !(canal 1_2 - delft 1980)

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1.)*dx*0.5

	if (x < 1.) then !Antes de subir
		auxx(i,j) = 0.2

	elseif (x < 1.4) then !Final da subida
		tgaux=-1./2.
		auxx(i,j)=tgaux*(x-1.) + 0.2

	elseif (x < 2.4) then !Em cima do obstáculo
		auxx(i,j)= 0.

	elseif (x < 2.8) then !Descida
		tgaux=1./2.
		auxx(i,j)=tgaux*(x-2.4)
		
	else !Resto
	   auxx(i,j) = 0.2
	endif
	enddo
	enddo

elseif (obst_t == 7) then !(SBRH - Buracos e calombos)

	raio = 0.12

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1.)*dx*0.5

	if (x < (1.-raio)) then !Antes de descer
		auxx(i,j) = 0.

	elseif (x < 1.) then !Final da descida

		auxx(i,j) = sqrt(-(x-1.)**2. +(raio**2.))

	else !Resto
	   auxx(i,j) = raio
 endif

	enddo
	enddo

elseif (obst_t == 8) then !(fennema, 1990)

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1.)*dx*0.5
	y = (j-1.)*dy*0.5

	if ((x >= 90.) .and. (x <= 100.)) then !Criando a parede
		auxx(i,j) = (nz+1)*dz
	endif
	
	if ((y > 95.) .and. (y < 170.)) then !Criando a parede
		auxx(i,j) = -dz
	endif
	
	enddo
	enddo

elseif (obst_t == 9) then !(aureli, 2008)

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1.)*dx*0.5
	y = (j-1.)*dy*0.5

	if ((x >= 0.8) .and. (x <= 0.825)) then !Criando a parede
		auxx(i,j) = (nz+1)*dz
	endif
	
	if ((y > 0.45) .and. (y < 0.75)) then !Criando a parede
		auxx(i,j) = -dz
	endif

	enddo
	enddo

elseif (obst_t == 10) then !(bd_koshizuka, 1995 e Kleefsman, 2005)

	auxx = -dz

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1.)*dx*0.5
	y = (j-1.)*dy*0.5

	if ((x >= (0.75-0.16*0.5) ) .and. (x <=(0.75+0.16*0.5) )) then !Criando a parede
		auxx(i,j) = 0.16!(nz+1)*dz
	endif
	
	if ((y > (0.5-0.4*0.5)) .and. (y < (0.5+0.4*0.5))) then !Criando a parede
		auxx(i,j) = -dz
	endif

	enddo
	enddo
	
elseif (obst_t == 12) then !rio Felipe

	auxx = -dz

	do j = -1,nyy+2
	do i = -1, nxx+2
	x = (i-1.)*dx*0.5
	y = (j-1.)*dy*0.5

	if ((y >= 1.5 ) .or. (y <=0.5 )) then !Criando a parede
		auxx(i,j) = 0.55!(nz+1)*dz
	else
		auxx(i,j) = 0.!(nz+1)*dz
	endif
	enddo
	enddo
	
endif


do j = 0, ny+1
do i = 0, nx1+1
	 ku(i,j) = nint(auxx(i*2-1,j*2)/dz+0.5)
	if (ku(i,j) > nz) ku(i,j) = nz
enddo
enddo

do j = 0, ny1+1
do i = 0, nx+1
	 kv(i,j) = nint(auxx(i*2,j*2-1)/dz+0.5)
	if (kv(i,j) > nz) kv(i,j) = nz
enddo
enddo

do j = 0, ny+1
do i = 0, nx+1
	 kw(i,j) = floor(auxx(i*2,j*2)/dz+1.)
	if (kw(i,j) > nz1) kw(i,j) = nz1
enddo
enddo

END SUBROUTINE obstaculo

SUBROUTINE sponge_layer(epis_z)

!Subrotina para representar a camada esponja, quando se quer um domínio com saída livre. 
!Não funciona com escoamentos, porque a velocidade é forçada a zero.

USE velpre
USE ls_param
IMPLICIT NONE

!Declarado também no programa

real(8), intent(out), dimension(nx,ny,nz1) :: epis_z

!Declarado apenas na rotina

!Profundidade do domínio
real(8), dimension(nx,ny,0:nz1) :: z_x
real(8), dimension(nx) :: zpfs
real(8), dimension(nx) :: xp,zpt

!Parâmetro da camada esponja (permeabilidade)
real(8) :: alfa, alfa0, alfa1, aux1, alt1, alt2

!Comprimento da camada esponja e posição do seu começo
real(8) :: l_sponge, x_inicial
integer :: i, j, k

!Inicialização do código
epis_z = 0.
z_x(:,:,0) = -dz
l_sponge  = 1. !Comprimento da camada esponja
x_inicial = nx*dx - l_sponge

alfa = 60. !A ser calibrado
alfa0 = 0. !A ser calibrado
alfa1 = 30. !2.*alfa !A ser calibrado

!Resolução do problema

if (esp_type == 1) then

	do j = 1, ny
	do i = nint(x_inicial/dx)+1, nx

	zpt(i)  = dz*nz  ! altura total

	do k = 1, nz-1
	
	if (ls(i,j,k)*ls(i,j,k+1) <= 0) then
		alt1 = dz*(k-0.5) + ls(i,j,k)
		alt2 = dz*(k+0.5) + ls(i,j,k+1)
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
	z_x(i,j,k) = z_x(i,j,k-1) + dz !Altura atual em w
	aux1 = (xp(i)-x_inicial)/l_sponge

	if (ls(i,j,k) >= 0.) then !Na água 
		epis_z(i,j,k) = alfa * aux1 * aux1 * (0.-z_x(i,j,k))     / (0.-zpfs(i)    )
	else
		epis_z(i,j,k) = alfa * aux1 * aux1 * (z_x(i,j,k)-zpt(i)) / (zpfs(i)-zpt(i))
	endif

	enddo
	enddo
	enddo

elseif (esp_type == 2) then

	do k = 1, nz
	do j = 1, ny
	do i = nint(x_inicial/dx)+1, nx
	xp(i) = real(i-0.5)*dx
	aux1 = (xp(i)-x_inicial)/l_sponge

	epis_z(i,j,k) = alfa0 + aux1 * (alfa1 - alfa0)

	enddo
	enddo
	enddo

elseif (esp_type == 3) then !Método da Tangente Hiperbólica

	do k = 1, nz
	do j = 1, ny
	do i = nint(x_inicial/dx)+1, nx
	xp(i) = real(i-0.5)*dx
	aux1 = (xp(i)-x_inicial)/l_sponge

	epis_z(i,j,k) = alfa* tanh(aux1*pi)

	enddo
	enddo
	enddo

endif

END SUBROUTINE sponge_layer

SUBROUTINE prd_corr(dpdx,dpdy,dpdz) !! arrumar rotina para eficiência!!
!Derivadas das pressões para adicionar nas condições de contorno (aproximar o valor em u^n+1 ...)

USE velpre
USE parametros

IMPLICIT NONE
!Declarado também no programa
real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: dpdx
real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: dpdy
real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: dpdz

real(8), dimension(nx1,ny,nz) :: rhox
real(8), dimension(nx,ny1,nz) :: rhoy
real(8), dimension(nx,ny,nz1) :: rhoz

integer :: i, j, k

call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
call interpz_cf(rho,nx,ny,nz,rhoz) !(nx,ny,nz1)


do k = 1, nz1
do j = 1, ny
do i = 1, nx

dpdz(i,j,k) = (prd(i,j,k) - prd(i,j,k-1))/ rhoz(i,j,k) *dt /dz

enddo
enddo
enddo

do k = 1, nz
do j = 1, ny1
do i = 1, nx

dpdy(i,j,k) = (prd(i,j,k) - prd(i,j-1,k))/ rhoy(i,j,k) *dt /dy

enddo
enddo
enddo

do k = 1, nz
do j = 1, ny
do i = 1, nx1

dpdx(i,j,k) = (prd(i,j,k) - prd(i-1,j,k))/ rhox(i,j,k) *dt /dx

enddo
enddo
enddo


dpdx(0,:,:) = dpdx(1,:,:) 
dpdy(0,:,:) = dpdy(1,:,:) 
dpdz(0,:,:) = dpdz(1,:,:) 

dpdx(nx1+1,:,:) = dpdx(nx1,:,:) 
dpdx(:,ny+1,:) = dpdx(:,ny,:) 
dpdx(:,:,nz+1) = dpdx(:,:,nz) 

dpdx(:,0,:) = dpdx(:,1,:) 
dpdy(:,0,:) = dpdy(:,1,:) 
dpdz(:,0,:) = dpdz(:,1,:) 

dpdy(nx+1,:,:) = dpdy(nx,:,:) 
dpdy(:,ny1+1,:) = dpdy(:,ny1,:) 
dpdy(:,:,nz+1) = dpdy(:,:,nz) 

dpdx(:,:,0) = dpdx(:,:,1) 
dpdy(:,:,0) = dpdy(:,:,1) 
dpdz(:,:,0) = dpdz(:,:,1) 

dpdz(nx+1,:,:) = dpdz(nx,:,:) 
dpdz(:,ny+1,:) = dpdz(:,ny,:) 
dpdz(:,:,nz1+1) = dpdz(:,:,nz1) 

END SUBROUTINE prd_corr

SUBROUTINE boundary_waves()

USE wave_c
USE parametros
USE velpre
USE ls_param

IMPLICIT NONE

integer :: i, j, k
real(8) :: aux1, aux2, aux3, aux4, aux5, h_fa, l_wa

real(8), dimension(0:nx) :: h_f

real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: u1
real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: v1
real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: w1

real(4),dimension(0:nx+1,0:ny+1,0:nz+1) :: ls1


!Reference: Coubilla, 2015 (Thesis)

ls1(1:nx,1:ny,1:nz) = ls(:,:,:)

ls1(0,:,:)    = ls1(1,:,:)
ls1(nx+1,:,:) = ls1(nx,:,:)

ls1(:,0,:)    = ls1(:,1,:)
ls1(:,ny+1,:) = ls1(:,ny,:)

ls1(:,:,0)    = ls1(:,:,1)
ls1(:,:,nz+1) = ls1(:,:,nz)


!Stokes I
if (wave_t == 1) then

	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, 1
	h_f(i)     = h0_f   + a_w*cos(n_w*dx*(i-0.5)-f_w*t)
	ls1(i,j,k) = h_f(i) - kp(k)

	if (ls1(i,j,k) >= 0.) then
		u1(i,j,k) = a_w * f_w * cosh(n_w*kp(k))/sinh(n_w*h0_f) * cos(n_w*dx*(i-1.)-f_w*t)
		w1(i,j,k) = a_w * f_w * sinh(n_w*(kp(k)-0.5*dz))/sinh(n_w*h0_f) * sin(n_w*dx*(i-0.5)-f_w*t)
	else
		u1(i,j,k) = 0.
		w1(i,j,k) = 0.	
	endif

	enddo
	enddo
	enddo

!Stokes II
elseif (wave_t == 2 ) then

	aux1 = tanh(n_w*h0_f) !Sigma

	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, 1

	aux2 = n_w*dx*(i-0.5)-f_w*t !Ômega para w e eta
	aux3 = n_w*dx*(i-1.)-f_w*t !Ômega para u

	h_f(i)     = h0_f + a_w*cos(aux2) +n_w*a_w*a_w*(3.-aux1*aux1)/(4.*aux1*aux1*aux1)*cos(2.*aux2)
	ls1(i,j,k) = h_f(i) - kp(k)

	if (ls1(i,j,k) >= 0.) then
		u1(i,j,k) = a_w*f_w*cosh(n_w*kp(k))/sinh(n_w*h0_f)*cos(aux3) &
		+3./4.*a_w*a_w*f_w*n_w*cosh(2.*n_w*kp(k))/((sinh(n_w*h0_f))**4.)*cos(2*aux3)

		w1(i,j,k) = a_w*f_w*sinh(n_w*(kp(k)-0.5*dz))/sinh(n_w*h0_f)*sin(aux2) &
		+3./4.*a_w*a_w*f_w*n_w*sinh(2.*n_w*(kp(k)-0.5*dz))/((sinh(n_w*h0_f))**4.)*sin(2.*aux3)
	else
		u1(i,j,k) = 0.
		w1(i,j,k) = 0.	
	endif
	
	enddo
	enddo
	enddo

!Stokes V
elseif (wave_t == 5 ) then

	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, 1

	aux2 = n_w*dx*(i-0.5)-f_w*t !Ômega para w e eta
	aux3 = n_w*dx*(i-1.)-f_w*t !Ômega para u

	aux4 = n_w*kp(k) !Para u
	aux5 = n_w*(kp(k)-0.5*dz) !Para w

	h_f(i) = h0_f +aeta1*cos(aux2) +aeta2*cos(2.*aux2) +aeta3*cos(3.*aux2) +aeta4*cos(4.*aux2) +aeta5*cos(5.*aux2)

	ls1(i,j,k) = h_f(i) - kp(k)

	if (ls1(i,j,k) >= 0.) then
		u1(i,j,k) = avel1*cosh(aux4)*cos(aux3) +avel2*cosh(2*aux4)*cos(2*aux3) +avel3*cosh(3*aux4)*cos(3*aux3) &
		+ avel4*cosh(4*aux4)*cos(4*aux3) +avel5*cosh(5*aux4)*cos(5*aux3)

		w1(i,j,k) = avel1*sinh(aux5)*sin(aux2) +avel2*sinh(2*aux5)*sin(2*aux2) +avel3*sinh(3*aux5)*sin(3*aux2) &
		+ avel4*sinh(4*aux5)*sin(4*aux2) +avel5*sinh(5*aux5)*sin(5*aux2)
	else
		u1(i,j,k) = 0.
		w1(i,j,k) = 0.	
	endif

	enddo
	enddo
	enddo

endif

ls(:,:,:) = ls1(1:nx,1:ny,1:nz)
bxx0(:,:) = u1(0,:,:)
bxy0(:,:) = v1(0,:,:)
bxz0(:,:) = w1(0,:,:)
bxx1(:,:) = u1(1,:,:)

!Streamfunction

!Solitary wave

END SUBROUTINE boundary_waves

SUBROUTINE waves_coef()
USE wave_c
USE parametros
IMPLICIT NONE

integer :: i, k, ii, nxii
real(8) :: aux1, aux2, aux3, aux4, aux5, h_fa, l_wa, lamb, lamb1, lamb2, erro, erro0, l_w1
real(8) :: s, c, a11, a13, a15, a22, a24, a33, a35, a44, a55, b22, b24, b33, b35, b44, b55, c1, c2

real(8), dimension(nx) :: h_f

!Reference: Coubilla, 2015 (Thesis)
do k = 0, nz+1
	kp(k) = (k-0.5)*dz
enddo

a_w  = 0.03 !Amplitude da onda
p_w  = 2.5  !Período da onda 
h0_f = 1.   !Profundidade do escoamento sem onda
h_f  = h0_f

l0_w = gz * p_w * p_w / (2.*pi)

!Stokes I e Stokes II
if (wave_t <= 2) then

	l_w = gz * p_w * p_w / (2.*pi)
	do i = 1, 1000
	l_w = gz * p_w * p_w / (2.*pi) * tanh(2.*pi*h0_f /l_w)
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
	s = sinh(2.*pi*h0_f/l_w1)
	c = cosh(2.*pi*h0_f/l_w1)

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
	lamb2 = sqrt((l_w1 /(l0_w*tanh(2.*pi*h0_f/l_w1)) -1.) / (c1 +lamb2**2.*c2))
	enddo

	erro = abs(lamb1 - lamb2)

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

	s = sinh(2.*pi*h0_f/l_w)
	c = cosh(2.*pi*h0_f/l_w)

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
	a55 = -(2880.*c**10.-72480.*c**8.+324000.*c**6.-432000.*c**4.+163470.*c*c-16245.)/(61440.*s**11.*(6.*c*c-1.)*(8.*c**4.-11.*c*c+3.))
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

