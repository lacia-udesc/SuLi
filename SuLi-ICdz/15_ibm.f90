!Subrotina definir IBM
!Referência: Auguste et al., 2019

!Implementação em 06/2023
!Leonardo Romero Monteiro

!############################################################

SUBROUTINE ibm_posicoes(id_ibm,obs_lsg,dimx,dimy,dimz)

!	USE disc
!	USE obst
	
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

	! limite na direção y
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


	! limite na direção z
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
	
END SUBROUTINE ibm_posicoes

!############################################################
SUBROUTINE ibm_images(obs_lsg,a1,dimx,dimy,dimz,id_ibm,der_t)

	USE disc, only : nx,ny,nz,dx,dy,dz,dzz,deltai

	IMPLICIT NONE
	
	integer :: i,j,k,l,ll,ii,jj,kk,i1,der_t
	integer :: dimx,dimy,dimz
	integer,dimension(dimx,dimy,dimz) :: id_ibm
	real(8),dimension(dimx,dimy,dimz) :: obs_lsg
	real(8),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: a1
	
	real(8),dimension(dimx,dimy,dimz) :: a1_image, norm_x, norm_y, norm_z, mod_ls
	
	real(8),dimension(dimx,dimy,dimz) :: a1_v, norm_x1, norm_y1, norm_z1, mod_lsa
	real(8),dimension(2) :: a1_imagep
	real(8),dimension(dimx) :: x
	real(8),dimension(dimy) :: y
	real(8),dimension(dimz) :: z
	real(8),dimension(8) :: lint
	real(8) :: v_sol, lg, li1, li2, sumlint,er, xl, yl, zl

	

	
	ll = 2 ! número de camadas imaginárias
	er = minval(deltai)*0.01 ! erro de aproximação da célula
	
	v_sol = 0. !velocidade alvo
	
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
	do k = 1, dimz
	do j = 1, dimy
	do i = 1, dimx

        ! inicia zerando tudo dentro do obstáculo !! não está parecendo ser necesssário
        if (obs_lsg(i,j,k) > 0. .and. id_ibm(i,j,k) .ne. 1 .and. der_t ==1) a1(i,j,k) = 0.
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
				xl = x(i) + (l*dx   +obs_lsg(i,j,k))*norm_x(i,j,k)
				yl = y(j) + (l*dy   +obs_lsg(i,j,k))*norm_y(i,j,k) 
				zl = z(k) + (l*dz(k)+obs_lsg(i,j,k))*norm_z(i,j,k)
				
				!arredondou para baixo. Assim sabemos que a média deve ser feita com o +1
				ii = floor((xl-x(i))/dx)+i
				jj = floor((yl-y(j))/dy)+j
				kk = floor((zl-z(k))/dz(k))+k
		    		! escolheu-se em fazer o cálculo da varipavel pelo método das distâncias invertidas! (Eq. 18 e texto posteriro no artigo)
		    		! cálculo da distância efetiva do ponto imaginário com as células em volta. A multiplicação verifica se está ou não dentro do ibm e usa a variável caso esteja fora do IBM.
				lint(1) = 1./(sqrt((xl-x(ii))**2.  +(yl-y(jj))**2.  +(zl-z(kk))**2.  ))*(1.-id_ibm(ii,jj,kk)      )
				lint(2) = 1./(sqrt((xl-x(ii+1))**2.+(yl-y(jj))**2.  +(zl-z(kk))**2.  ))*(1.-id_ibm(ii+1,jj,kk)    )
				lint(3) = 1./(sqrt((xl-x(ii))**2.  +(yl-y(jj+1))**2.+(zl-z(kk))**2.  ))*(1.-id_ibm(ii,jj+1,kk)    )
				lint(4) = 1./(sqrt((xl-x(ii))**2.  +(yl-y(jj))**2.  +(zl-z(kk+1))**2.))*(1.-id_ibm(ii,jj,kk+1)    )
				lint(5) = 1./(sqrt((xl-x(ii+1))**2.+(yl-y(jj+1))**2.+(zl-z(kk))**2.  ))*(1.-id_ibm(ii+1,jj+1,kk)  )
				lint(6) = 1./(sqrt((xl-x(ii+1))**2.+(yl-y(jj))**2.  +(zl-z(kk+1))**2.))*(1.-id_ibm(ii+1,jj,kk+1)  )
				lint(7) = 1./(sqrt((xl-x(ii))**2.  +(yl-y(jj+1))**2.+(zl-z(kk+1))**2.))*(1.-id_ibm(ii,jj+1,kk+1)  )
				lint(8) = 1./(sqrt((xl-x(ii+1))**2.+(yl-y(jj+1))**2.+(zl-z(kk+1))**2.))*(1.-id_ibm(ii+1,jj+1,kk+1))
			    
			    	! soma as distância para fazer a ponderação depois
				sumlint = sum(lint(:))
	
				! usa a distância e define a variável ! eq. 17 do artigo
				lint(1) = lint(1)*a1(ii,jj,kk)
				lint(2) = lint(2)*a1(ii+1,jj,kk)
				lint(3) = lint(3)*a1(ii,jj+1,kk)			
				lint(4) = lint(4)*a1(ii,jj,kk+1)
				lint(5) = lint(5)*a1(ii+1,jj+1,kk)            
				lint(6) = lint(6)*a1(ii+1,jj,kk+1)            
				lint(7) = lint(7)*a1(ii,jj+1,kk+1)
				lint(8) = lint(8)*a1(ii+1,jj+1,kk+1)
			     
				!interpolação para obter a imagem !fecha a equação 18.
				a1_imagep(l) = sum(lint(:))/sumlint

				if (sumlint == 0) a1_imagep(l) = 0.

			enddo
			
			!!@@## interpolador a (eq. 14 do artigo)
			li1 = (2.*deltai(k) - obs_lsg(i,j,k))/deltai(k) * (2.*obs_lsg(i,j,k))/(deltai(k) + obs_lsg(i,j,k))
			li2 = (obs_lsg(i,j,k) - deltai(k))/deltai(k) * (2.*obs_lsg(i,j,k))/(2.*deltai(k) + obs_lsg(i,j,k))
			lg  = (obs_lsg(i,j,k) - deltai(k))/(obs_lsg(i,j,k) + deltai(k))*(obs_lsg(i,j,k)-2.*deltai(k))/(obs_lsg(i,j,k)+2.*deltai(k))
		
			! cálculo das imagens que juntam as imagens 1 e 2 (eq. 12 do artigo)
			if (der_t == 1) a1_image(i,j,k) = (2.*lg*v_sol        + li1*a1_imagep(1) + li2*a1_imagep(2))/(1.+lg)
			if (der_t == 2) a1_image(i,j,k) = (2.*lg*a1_imagep(1) + li1*a1_imagep(1) + li2*a1_imagep(2))/(1.+lg)
			
			!!@@## interpolador b
			!li1 = (2.*deltai(k) - obs_lsg(i,j,k))/deltai(k)*obs_lsg(i,j,k)/deltai(k)
			!li2 = (obs_lsg(i,j,k) - deltai(k))/deltai(k) * (obs_lsg(i,j,k))/(2.*deltai(k))
			!lg  = (obs_lsg(i,j,k))/deltai(k) *(obs_lsg(i,j,k)-2.*deltai(k))/(2.*deltai(k))
		
			! cálculo das imagens que juntam as imagens 1 e 2
			!if (der_t == 1) a1_image(i,j,k) = lg*v_sol + li1*a1_imagep(1) + li2*a1_imagep(2)
			!if (der_t == 2) a1_image(i,j,k) = lg*a1_imagep(1) + li1*a1_imagep(1) + li2*a1_imagep(2)
			
			!!@@## interpolador c
			!li1 = 1./(deltai(k)    -obs_lsg(i,j,k))
			!li2 = 1./(2.*deltai(k) -obs_lsg(i,j,k))
			!lg  = 1./obs_lsg(i,j,k) 
		
			! cálculo das imagens que juntam as imagens 1 e 2
			!if (der_t == 1)	a1_image(i,j,k) = (lg*v_sol + li1*a1_imagep(1) + li2*a1_imagep(2))/(li1+li2+lg)
			!if (der_t == 2) a1_image(i,j,k) = (lg*a1_imagep(1) + li1*a1_imagep(1) + li2*a1_imagep(2))/(li1+li2+lg)
    		endif

		! eq. 21 do artigo
		if (der_t == 1) a1(i,j,k) =  2.*v_sol - a1_image(i,j,k)
		if (der_t == 2) a1(i,j,k) =  a1_image(i,j,k) !+ 2.*obs_lsg(i,j,k)*1.	! simplifiquei o Neumann para evitar o cálculo da derivada (eq. 21 do artigo do Auguste)
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
	
END SUBROUTINE ibm_images

