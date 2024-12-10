!Subrotina para calcular a pressão dinâmica pelo método do gradiente conjugado
!Referencia: Casulli (1992 e 1999)

!!! Implementação 15/04/2014
!Leonardo Romero Monteiro

!!! Modificações
!Leonardo Romero Monteiro em 13/01/2015

SUBROUTINE graddin()

	USE velpre
	USE param
	USE cond
	USE obst
	USE ls_param, only: rho_f1
	!$ USE omp_lib
	
	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8),save :: alfapr, alfamupr, alfadipr, betapr, betamupr
	real(8),dimension(nx1,ny,nz) :: matspri
	real(8),dimension(nx,ny1,nz) :: matsprj
	real(8),dimension(nx,ny,nz1) :: matsprk
	real(8),dimension(nx,ny,nz) :: matqpr
	real(8),dimension(0:nx1,0:ny1,0:nz1) :: matdpr, matepr, erropr, erroppr, mppr

	!contadores
	integer :: i, j, k, cont

	!auxiliares
	real(8),save :: aux1, aux2, aux3

	real(8),dimension(nx1,ny,nz) :: rhox, u1
	real(8),dimension(nx,ny1,nz) :: rhoy, v1
	real(8),dimension(nx,ny,nz1) :: rhoz, w1

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	    if (obst_t .ne. 0 .and. ibm_t == 2) then
		!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k) 
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
				if (obs_lsx(i,j,k) .ge. 0. ) then
					u1(i,j,k) = 0.
				else
					u1(i,j,k) = u(i,j,k)
				endif
			enddo
			enddo
			enddo
		!$OMP END PARALLEL DO
		
		!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k) 
			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
				if (obs_lsy(i,j,k) .ge. 0. ) then
					v1(i,j,k) = 0.
				else
					v1(i,j,k) = v(i,j,k)
				endif
			enddo
			enddo
			enddo
		!$OMP END PARALLEL DO
		
		!$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i,j,k) 
			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
				if (obs_lsz(i,j,k) .ge. 0. ) then
					w1(i,j,k) = 0.
				else
					w1(i,j,k) = w(i,j,k)
				endif
			enddo
			enddo
			enddo
		!$OMP END PARALLEL DO
	
	   endif
	   
   
	!%%%!-- Método do Gradiente Conjugado - Para Pressão Dinâmica --!%%%!

	!$OMP PARALLEL SECTIONS
	!$OMP  SECTION
	aux1 = dt / (dx*dx)
	!$OMP  SECTION
	aux2 = dt / (dy*dy) 
	!$OMP END PARALLEL SECTIONS


	if  (t_press .eq. 1) then !Projection Method
	
		call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
		call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
		call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)

	else !Implicit Pressure Algorithm (IPA) 
	!! obs.: seria possível reduzir o tempo do IPA deixando o matspr, matdpr e matapr... como constantes, calculando-os apenas uma vez (ver código do IPA do artigo). Caso precise de mais desempenho (estimado em 1%) fazer a modificação. Assim, optou-se em fazer uma mudança pouco eficiente, mas mais fácil de trabalhar.
	
		rhox=rho_f1
		rhoy=rho_f1
		rhoz=rho_f1
	
	endif

	!$OMP PARALLEL 
	!$OMP DO PRIVATE(i,j,k) 
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			matspri(i,j,k) = aux1/rhox(i,j,k)
		enddo
		enddo
		enddo
	!$OMP END DO NOWAIT

	!$OMP DO PRIVATE(i,j,k) 
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
			matsprj(i,j,k) = aux2/rhoy(i,j,k)
		enddo
		enddo
		enddo
	!$OMP END DO NOWAIT

	!$OMP DO PRIVATE(i,j,k) 
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			matsprk(i,j,k) = dt /(dz(k)*dz(k)*rhoz(i,j,k))
		enddo
		enddo
		enddo
	!$OMP END DO
	!$OMP END PARALLEL


	! Matrizes p e q
	!$OMP PARALLEL DO PRIVATE(i,j,k) 
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		matdpr(i,j,k) = matspri(i+1,j,k) + matspri(i,j,k) + matsprj(i,j+1,k) + matsprj(i,j,k) + matsprk(i,j,k+1) + matsprk(i,j,k)
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO 


	if (obst_t .ne. 0 .and. ibm_t == 2) then
		!$OMP PARALLEL DO PRIVATE(i,j,k) 
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			matqpr(i,j,k) = ( u1(i,j,k) - u1(i+1,j,k) )/dx + ( v1(i,j,k)-v1(i,j+1,k) )/dy  + (w1(i,j,k) - w1(i,j,k+1)) /dzz(k)
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO 
	else
		!$OMP PARALLEL DO PRIVATE(i,j,k) 
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			matqpr(i,j,k) = ( u(i,j,k) - u(i+1,j,k) )/dx + ( v(i,j,k)-v(i,j+1,k) )/dy  + (w(i,j,k) - w(i,j,k+1)) /dzz(k)
		enddo
		enddo
		enddo
		!$OMP END PARALLEL DO 
	endif


	!dentro do obstáculo, divergência nula
    
	!$OMP PARALLEL
	! Condições de contorno, von Neumann
    
	if (ccx0.eq.0) then  ! Condição periódica
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz1
		do j = 0, ny1
			matdpr(0,j,k) = matdpr(nx,j,k)
			matdpr(nx+1,j,k) = matdpr(1,j,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	else
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz1		
		do j = 0, ny1
			matdpr(0,j,k) = matdpr(1,j,k)
			matdpr(nx+1,j,k) = matdpr(nx,j,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	endif

	if (ccy0.eq.0) then  ! Condição periódica
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz1		
		do i = 0, nx1
			matdpr(i,0,k) = matdpr(i,ny,k)
			matdpr(i,ny+1,k) = matdpr(i,1,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	else
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz1		
		do i = 0, nx1
			matdpr(i,0,k)  = matdpr(i,1,k)
			matdpr(i,ny+1,k) = matdpr(i,ny,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	endif

		!$OMP DO PRIVATE(i,j,k) 
		do j = 0, ny1
		do i = 0, nx1
			matdpr(i,j,0)   = matdpr(i,j,1)
			matdpr(i,j,nz+1) = matdpr(i,j,nz)
		enddo
		enddo
		!$OMP END DO
		!$OMP END PARALLEL
	
	! Normalização da pressão dinâmica
	!$OMP PARALLEL
	!$OMP DO PRIVATE(i,j,k) 
	!SHARED(matspr,matdpr,matdpr)
	do k = 0, nz1
	do j = 0, ny1
	do i = 0, nx1
		matepr(i,j,k) = prd1(i,j,k) * sqrt(matdpr(i,j,k))
	enddo
	enddo
	enddo
	!$OMP END DO NOWAIT

	!Normalização das matrizes s e cálculo do primeiro erro
	
	!$OMP DO PRIVATE(i,j,k) 
	!SHARED(matspri,matdpr,matdpr)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		matspri(i,j,k) = matspri(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i-1,j,k))
	enddo
	enddo
	enddo
	!$OMP END DO NOWAIT

	!$OMP DO PRIVATE(i,j,k) 
	!SHARED(matsprj,matdpr,matdpr)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		matsprj(i,j,k) = matsprj(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i,j-1,k))
	enddo
	enddo
	enddo
	!$OMP END DO NOWAIT
	
	!$OMP DO PRIVATE(i,j,k) 
	!SHARED(matsprk,matdpr,matdpr)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		matsprk(i,j,k) = matsprk(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i,j,k-1))
	enddo
	enddo
	enddo
	!$OMP END DO
	!$OMP END PARALLEL
	
	!$OMP PARALLEL DO PRIVATE(i,j,k) 
	!SHARED(erropr,matepr,matspri,matsprj,matsprk,matqpr,matdpr)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		erropr(i,j,k) = matepr(i,j,k) - matspri(i+1,j,k) * matepr(i+1,j,k) - matspri(i,j,k) * matepr(i-1,j,k) &
			- matsprj(i,j+1,k) * matepr(i,j+1,k) - matsprj(i,j,k) * matepr(i,j-1,k) &
			- matsprk(i,j,k+1) * matepr(i,j,k+1) - matsprk(i,j,k) * matepr(i,j,k-1) - matqpr(i,j,k)/sqrt(matdpr(i,j,k))
	enddo
	enddo
	enddo
	!$OMP END PARALLEL DO
	
	!$OMP PARALLEL	
	
	!$OMP SINGLE
	alfamupr = 0.
	!$OMP END SINGLE
	
	!$OMP SINGLE
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		alfamupr = alfamupr + erropr(i,j,k) * erropr(i,j,k)
	enddo
	enddo
	enddo
	!$OMP END SINGLE
	
	!$OMP END PARALLEL
	
			
	!$OMP PARALLEL
	if (ccx0.eq.0) then  ! Condição periódica
		!$OMP  DO PRIVATE(i,j,k) 
		do k = 0, nz1
		do j = 0, ny1		
			erropr(0,j,k) = erropr(nx,j,k)
			erropr(nx+1,j,k) = erropr(1,j,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	else
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz1		
		do j = 0, ny1
			erropr(0,j,k) = erropr(1,j,k)
			erropr(nx+1,j,k) = erropr(nx,j,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	endif

	if (ccy0.eq.0) then  ! Condição periódica
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz1		
		do i = 0, nx1
			erropr(i,0,k) = erropr(i,ny,k)
			erropr(i,ny+1,k) = erropr(i,1,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	else
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz1		
		do i = 0, nx1
			erropr(i,0,k) = erropr(i,1,k)
			erropr(i,ny+1,k) = erropr(i,ny,k)
		enddo
		enddo
		!$OMP END DO NOWAIT
	endif

		!$OMP DO PRIVATE(i,j,k) 
		do j = 0, ny1
		do i = 0, nx1
			erropr(i,j,0) = erropr(i,j,1)
			erropr(i,j,nz+1) = erropr(i,j,nz)
		enddo
		enddo
		!$OMP END DO
		
		!$OMP DO PRIVATE(i,j,k) 
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, nx+1
			erroppr(i,j,k) = erropr(i,j,k)
		enddo
		enddo
		enddo		
		!$OMP END DO
				
	!$OMP END PARALLEL
	
	!$OMP PARALLEL SECTIONS

		!$OMP  SECTION
		cont = 0

		!$OMP  SECTION
		aux1 = abs(alfamupr)
		
		!$OMP  SECTION
		aux2 = 	(0.0001/(nx*ny*nz))	
		
	!$OMP END PARALLEL SECTIONS
	
	!%%%%%%%%%%%%%   loop da redução do erro   %%%%%%%%%%%%%%!
	do while ((aux1 > aux2) .and. (cont < 10000) )
!inicialização
	!$OMP PARALLEL
		!$OMP  WORKSHARE
			cont = cont +1	
			alfadipr = 0.
			betamupr = 0.
		!$OMP  END WORKSHARE
	!$OMP END PARALLEL

	 !$OMP PARALLEL DO PRIVATE(i,j,k)
			!SHARED(mppr,erroppr,matspri,matsprj,matsprk) 
			! Parâmetro mp e alfa
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				mppr(i,j,k) = erroppr(i,j,k) - erroppr(i+1,j,k)*matspri(i+1,j,k) - erroppr(i-1,j,k)*matspri(i,j,k) &
				- erroppr(i,j+1,k)*matsprj(i,j+1,k) - erroppr(i,j-1,k)*matsprj(i,j,k) & 
				- erroppr(i,j,k+1)*matsprk(i,j,k+1) - erroppr(i,j,k-1)*matsprk(i,j,k)
			enddo
			enddo
			enddo
	!$OMP END PARALLEL DO	

	if (omp_t == 0) then
	!$OMP PARALLEL	
		!$OMP SINGLE 
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				alfadipr = alfadipr + erroppr(i,j,k) * mppr(i,j,k)
			enddo
			enddo
			enddo
		!$OMP END SINGLE 
	!$OMP END PARALLEL


	else
	!$OMP PARALLEL	
		!$OMP DO PRIVATE(i,j,k) reduction (+:alfadipr)
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				alfadipr = alfadipr + erroppr(i,j,k) * mppr(i,j,k)
			enddo
			enddo
			enddo
		!$OMP END DO 
	!$OMP END PARALLEL
	endif

	!$OMP PARALLEL	
		!$OMP SINGLE			
			alfapr = alfamupr / alfadipr
		!$OMP END SINGLE
					
	!$OMP END PARALLEL

			! Recálculo das matrizes e, erro e parâmetro beta
	!$OMP PARALLEL
		!$OMP DO PRIVATE(i,j,k) 
			!SHARED(matepr,alfapr,erroppr)
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				matepr(i,j,k)  = matepr(i,j,k)  - alfapr * erroppr(i,j,k)
			enddo
			enddo
			enddo
		!$OMP END DO NOWAIT

		!$OMP DO PRIVATE(i,j,k) 
			!SHARED(erropr,alfapr,mppr)
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				erropr(i,j,k) = erropr(i,j,k) - alfapr * mppr(i,j,k)
			enddo
			enddo
			enddo
		!$OMP END DO
	!$OMP END PARALLEL

	if (omp_t == 0) then
	!$OMP PARALLEL	
		!$OMP SINGLE 
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				betamupr = betamupr + erropr(i,j,k) * erropr(i,j,k)
			enddo
			enddo
			enddo
		!$OMP END SINGLE 
	!$OMP END PARALLEL


	else
	!$OMP PARALLEL	
		!$OMP DO PRIVATE(i,j,k) reduction (+:betamupr)
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				betamupr = betamupr + erropr(i,j,k) * erropr(i,j,k)
			enddo
			enddo
			enddo
		!$OMP END DO 
	!$OMP END PARALLEL
	endif
			
	!$OMP PARALLEL
	!$OMP WORKSHARE
	betapr = betamupr/alfamupr
	alfamupr = betamupr
	!$OMP END WORKSHARE
	!$OMP END PARALLEL
							
			! Recálculo de erroppr
	!$OMP PARALLEL DO PRIVATE(i,j,k) 
			!SHARED(erroppr,erropr,betapr)
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				erroppr(i,j,k) = erropr(i,j,k) + betapr * erroppr(i,j,k)
			enddo
			enddo
			enddo
			!$OMP END PARALLEL DO

			! Condições de contorno
			
		!$OMP PARALLEL
			if (ccx0.eq.0) then  ! Condição periódica
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1
				do j = 0, ny1
					erroppr(0,j,k) = erroppr(nx,j,k)
					erroppr(nx+1,j,k) = erroppr(1,j,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			else
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1		
				do j = 0, ny1
					erroppr(0,j,k) = erroppr(1,j,k)
					erroppr(nx+1,j,k) = erroppr(nx,j,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			endif

			if (ccy0.eq.0) then  ! Condição periódica
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1		
				do i = 0, nx1
					erroppr(i,0,k) = erroppr(i,ny,k)
					erroppr(i,ny+1,k) = erroppr(i,1,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			else
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1		
				do i = 0, nx1
					erroppr(i,0,k) = erroppr(i,1,k)
					erroppr(i,ny+1,k) = erroppr(i,ny,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			endif

			!$OMP DO PRIVATE(i,j,k) 
				do j = 0, ny1
				do i = 0, nx1
					erroppr(i,j,0) = erroppr(i,j,1)
					erroppr(i,j,nz+1) = erroppr(i,j,nz)
				enddo
				enddo
			!$OMP END DO
				
			!$OMP SINGLE
			aux1 = abs(alfamupr)
			!$OMP END SINGLE
		
		
		!$OMP END PARALLEL
		enddo

		! se pular pressão ele vai avisar
		if (cont == 10000) write(*,*) "pulou pressão; ", "erro =", aux1

		! Condições de contorno
		
		!$OMP PARALLEL
			if (ccx0.eq.0) then  ! Condição periódica
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1
				do j = 0, ny1
					matepr(0,j,k) = matepr(nx,j,k)
					matepr(nx+1,j,k) = matepr(1,j,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			else
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1		
				do j = 0, ny1
					matepr(0,j,k) = matepr(1,j,k)
					matepr(nx+1,j,k) = matepr(nx,j,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			endif

			if (ccy0.eq.0) then  ! Condição periódica
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1		
				do i = 0, nx1
					matepr(i,0,k) = matepr(i,ny,k)
					matepr(i,ny+1,k) = matepr(i,1,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			else
			!$OMP DO PRIVATE(i,j,k) 
				do k = 0, nz1		
				do i = 0, nx1
					matepr(i,0,k) = matepr(i,1,k)
					matepr(i,ny+1,k) = matepr(i,ny,k)
				enddo
				enddo
			!$OMP END DO NOWAIT
			endif

			!$OMP DO PRIVATE(i,j,k) 
				do j = 0, ny1
				do i = 0, nx1
					matepr(i,j,0) = matepr(i,j,1)
					matepr(i,j,nz+1) = matepr(i,j,nz)			
				enddo
				enddo
			!$OMP END DO
		!$OMP END PARALLEL
		
		
		!write(*,*) "o erro do Gradiente Conjugado foi", alfamupr
		
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%
		! Desnormalização de matriz e para a pressão dinâmica
	!$OMP PARALLEL DO PRIVATE(i,j,k) 
		!SHARED(prd1,matepr,matdpr)
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, nx+1
			prd1(i,j,k) = matepr(i,j,k) / sqrt(matdpr(i,j,k))
		enddo
		enddo
		enddo
	!$OMP END PARALLEL DO

	!===============================================================================================================

	END SUBROUTINE graddin


!######################################################################################
SUBROUTINE pressh()

	USE velpre
	USE param

	IMPLICIT NONE

	integer :: i, j, k
	real(8),save :: intz

	!WORK Z-PENCILS

	   do k = nz,1 ,-1
	   do j = 1, ny
	   do i = 1, nx
	      if (k == nz) then !top boundary conditions
		 prd1(i,j,k)   = dz(k)*gz*rho(i,j,k)

	      elseif (k == nz-1) then !top boundary conditions

		 prd1(i,j,k)   = gz*(rho(i,j,k)*dz(k)+ rho(i,j,k+1)*dz(k+1))*0.5 + prd1(i,j,k+1)

	      elseif (k < nz-1) then

		 intz = (rho(i,j,k)*dz(k) + 4.*rho(i,j,k+1)*dz(k+1) + rho(i,j,k+2)*dz(k+2)) /3.
		 prd1(i,j,k)   = gz * intz + prd1(i,j,k+2) 

	      endif

	   enddo
	   enddo
	   enddo



END SUBROUTINE pressh

