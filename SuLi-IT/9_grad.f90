!Subrotina para calcular a pressão dinâmica pelo método do gradiente conjugado
!Referencia: Casulli (1992 e 1999)

!!! Implementação 15/04/2014
! Leonardo Romero Monteiro

!!! Modificações
! Leonardo Romero Monteiro - 13/01/2015

SUBROUTINE graddin()

	USE dzs
	USE velpre
	USE parametros
        USE cond

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8) :: alfapr, alfamupr, alfadipr, betapr, betamupr
	real(8), dimension(nx1,ny,nz) :: matspri
	real(8), dimension(nx,ny1,nz) :: matsprj
	real(8), dimension(nx,ny,nz1) :: matsprk
	real(8), dimension(nx,ny,nz) :: matqpr, matapripos, mataprineg, mataprjpos, mataprjneg, mataprkpos, mataprkneg
	real(8), dimension(0:nx1,0:ny1,0:nz1) :: matdpr, matepr, erropr, erroppr, mppr

	!contadores
	integer :: i, j, k

	!auxiliares
	real(8) :: aux1, aux2, aux3
	
	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================


	!%%%!-- Método do Gradiente Conjugado - Para Pressão Dinâmica --!%%%!

		aux1 = tetah * dt / (dx*dx)
		aux2 = tetah * dt / (dy*dy) 

		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			matspri(i,j,k) = aux1 * dzx(i,j,k)
		enddo

		do i = 1, nx
			matsprj(i,j,k) = aux2 * dzy(i,j,k) 
			matsprk(i,j,k) = tetah *dt /dzz(i,j,k)
		enddo
		enddo
	
		j = ny1
		do i = 1, nx
			matsprj(i,j,k) = aux2 * dzy(i,j,k) 
		enddo
		enddo

		matsprk(:,:,nz1) = 0.

		! Matrizes p e q
		do k = 1, nz-1
		do j = 1, ny
		do i = 1, nx
			matdpr(i,j,k) = matspri(i+1,j,k) + matspri(i,j,k) + matsprj(i,j+1,k) + matsprj(i,j,k) + matsprk(i,j,k+1) + matsprk(i,j,k)
			matqpr(i,j,k) = ( u(i,j,k)*dzx(i,j,k) - u(i+1,j,k)*dzx(i+1,j,k) )/dx &
			+ ( v(i,j,k)*dzy(i,j,k)-v(i,j+1,k)*dzy(i,j+1,k) )/dy  + w(i,j,k) - w(i,j,k+1)
		enddo
		enddo
		enddo

		k = nz
		do j = 1, ny
		do i = 1, nx
			! Descontinuidade da equação na superfície livre por causa da variação na Equação da Continuidade
			matdpr(i,j,nz) = matspri(i+1,j,nz) + matspri(i,j,nz) + matsprj(i,j+1,nz) + matsprj(i,j,nz) + matsprk(i,j,nz) + 1./(gz*tetah*dt) !+ 1./(gz*dt) !
			matqpr(i,j,nz) = ( u(i,j,nz)*dzx(i,j,nz) - u(i+1,j,nz)*dzx(i+1,j,nz) )/dx &
			+( v(i,j,nz)*dzy(i,j,nz) - v(i,j+1,nz)*dzy(i,j+1,nz) )/dy + w(i,j,nz) - (eta1(i,j)-deltah(i,j))/(tetah*dt) !- (eta1(i,j)-eta0(i,j))/dt !
		enddo
		enddo

		! Condições de contorno, von Neumann
		if (ccx0.eq.0) then  ! Condição periódica
			matdpr(0,:,:)   = matdpr(nx,:,:)
			matdpr(nx1,:,:) = matdpr(1,:,:)
		else
			matdpr(0,:,:)   = matdpr(1,:,:)
			matdpr(nx1,:,:) = matdpr(nx,:,:)
		endif

		if (ccy0.eq.0) then  ! Condição periódica
			matdpr(:,0,:)   = matdpr(:,ny,:)
			matdpr(:,ny1,:) = matdpr(:,1,:)
		else
			matdpr(:,0,:)   = matdpr(:,1,:)
			matdpr(:,ny1,:) = matdpr(:,ny,:)
		endif


		matdpr(:,:,0)   = matdpr(:,:,1)
		matdpr(:,:,nz1) = matdpr(:,:,nz)
		! Normalização da pressão dinâmica
		do k = 0, nz1
		do j = 0, ny1
		do i = 0, nx1
			matepr(i,j,k) = prd1(i,j,k) * sqrt(matdpr(i,j,k))
		enddo
		enddo
		enddo

		!Normalização das matrizes s e cálculo do primeiro erro
		erropr =   0.
		alfamupr = 0.
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx

			matapripos(i,j,k) = matspri(i+1,j,k)/sqrt(matdpr(i,j,k)*matdpr(i+1,j,k))
			mataprineg(i,j,k) = matspri(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i-1,j,k))

			mataprjpos(i,j,k) = matsprj(i,j+1,k)/sqrt(matdpr(i,j,k)*matdpr(i,j+1,k))
			mataprjneg(i,j,k) = matsprj(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i,j-1,k))

			mataprkpos(i,j,k) = matsprk(i,j,k+1)/sqrt(matdpr(i,j,k)*matdpr(i,j,k+1))
			mataprkneg(i,j,k) = matsprk(i,j,k)/sqrt(matdpr(i,j,k)*matdpr(i,j,k-1))


			erropr(i,j,k) = matepr(i,j,k) - matapripos(i,j,k) * matepr(i+1,j,k) - mataprineg(i,j,k) * matepr(i-1,j,k) &
			- mataprjpos(i,j,k) * matepr(i,j+1,k) - mataprjneg(i,j,k) * matepr(i,j-1,k) &
			- mataprkpos(i,j,k) * matepr(i,j,k+1) - mataprkneg(i,j,k) * matepr(i,j,k-1) - matqpr(i,j,k)/sqrt(matdpr(i,j,k))


			alfamupr = alfamupr + erropr(i,j,k) * erropr(i,j,k)
		enddo
		enddo
		enddo

		if (ccx0.eq.0) then  ! Condição periódica
			erropr(0,:,:)   = erropr(nx,:,:)
			erropr(nx1,:,:) = erropr(1,:,:)
		else
			erropr(0,:,:)   = erropr(1,:,:)
			erropr(nx1,:,:) = erropr(nx,:,:)
		endif

		if (ccy0.eq.0) then  ! Condição periódica
			erropr(:,0,:)   = erropr(:,ny,:)
			erropr(:,ny1,:) = erropr(:,1,:)
		else
			erropr(:,0,:)   = erropr(:,1,:)
			erropr(:,ny1,:) = erropr(:,ny,:)
		endif

		erropr(:,:,0)   = erropr(:,:,1)
		erropr(:,:,nz1) = erropr(:,:,nz)

		erroppr = erropr

		!%%%%%%%%%%%%%   loop da redução do erro   %%%%%%%%%%%%%%!
		do while (abs(alfamupr) > 0.00000000001)


			!inicialização
			alfapr   = 0.
			alfamupr = 0.
			alfadipr = 0.
			betapr   = 0.
			betamupr = 0.
			mppr     = 0.

			! Parâmetro mp e alfa

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				mppr(i,j,k) = erroppr(i,j,k) - erroppr(i+1,j,k) * matapripos(i,j,k) - erroppr(i-1,j,k) * mataprineg(i,j,k) &
				- erroppr(i,j+1,k) * mataprjpos(i,j,k) - erroppr(i,j-1,k) * mataprjneg(i,j,k) & 
				- erroppr(i,j,k+1) * mataprkpos(i,j,k) - erroppr(i,j,k-1) * mataprkneg(i,j,k)

				alfamupr = alfamupr + erropr(i,j,k) * erropr(i,j,k)
				alfadipr = alfadipr + erroppr(i,j,k) * mppr(i,j,k)
			enddo
			enddo
			enddo

			alfapr = alfamupr / alfadipr
	
			! Recálculo das matrizes e, erro e parâmetro beta


			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				matepr(i,j,k)  = matepr(i,j,k)  - alfapr * erroppr(i,j,k)
				erropr(i,j,k) = erropr(i,j,k) - alfapr * mppr(i,j,k)

				betamupr = betamupr + erropr(i,j,k) * erropr(i,j,k)
			enddo
			enddo
			enddo

			betapr = betamupr/alfamupr

			! Recálculo de errop
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				erroppr(i,j,k) = erropr(i,j,k) + betapr * erroppr(i,j,k)
			enddo
			enddo
			enddo

			! Condições de contorno

		if (ccx0.eq.0) then  ! Condição periódica
			matepr(0,:,:)   = matepr(nx,:,:)
			matepr(nx1,:,:) = matepr(1,:,:)
			erroppr(0,:,:)   = erroppr(nx,:,:)
			erroppr(nx1,:,:) = erroppr(1,:,:)
		else
			matepr(0,:,:)   = matepr(1,:,:)
			matepr(nx1,:,:) = matepr(nx,:,:)
			erroppr(0,:,:)   = erroppr(1,:,:)
			erroppr(nx1,:,:) = erroppr(nx,:,:)
		endif

		if (ccy0.eq.0) then  ! Condição periódica
			matepr(:,0,:)   = matepr(:,ny,:)
			matepr(:,ny1,:) = matepr(:,1,:)
			erroppr(:,0,:)   = erroppr(:,ny,:)
			erroppr(:,ny1,:) = erroppr(:,1,:)
		else
			matepr(:,0,:)   = matepr(:,1,:)
			matepr(:,ny1,:) = matepr(:,ny,:)
			erroppr(:,0,:)   = erroppr(:,1,:)
			erroppr(:,ny1,:) = erroppr(:,ny,:)
		endif


			matepr(:,:,0)   = matepr(:,:,1)
			matepr(:,:,nz1) = matepr(:,:,nz)


			erroppr(:,:,0)   = erroppr(:,:,1)
			erroppr(:,:,nz1) = erroppr(:,:,nz)

		enddo
		!%%%%%%%%%%%%%%%%%%%%%%%%%%%%

		! Desnormalização de matriz e para a pressão dinâmica
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, nx+1
			prd1(i,j,k) = matepr(i,j,k) / sqrt(matdpr(i,j,k))
		enddo
		enddo
		enddo



END SUBROUTINE graddin

!===============================================================================================================
		
!!! Subrotina para calcular o desnível não corrigido pelo método do gradiente conjugado
!Referencia: Casulli (1992)

!!! Implementação 15/04/2014
! Leonardo Romero Monteiro

!!! Modificações
! Leonardo Romero Monteiro - 13/01/2015


SUBROUTINE gradeta(parte1, parte2, parte3, parte4)

	USE velpre
	USE velpre
	USE parametros
        USE cond

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA

	real(8), intent(in), dimension(nx1,ny) :: parte1, parte2
	real(8), intent(in), dimension(nx,ny1) :: parte3, parte4

	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)
	real(8) :: alfa, alfamu, alfadi, beta, betamu
	real(8), dimension(nx1,ny) :: matrizsi
	real(8), dimension(nx,ny1) :: matrizsj
	real(8), dimension(nx,ny) :: matrizq, matrizadipos, matrizadineg, matrizadjpos, matrizadjneg
	real(8), dimension(0:nx1,0:ny1) :: matrizd, matrize, erro, errop, mp

	!contadores
	integer :: i,j

	!auxiliares
	real(8) :: aux1, aux2

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	CALL deltah_()


	!===============!%%%!-- Método do Gradiente Conjugado - Para Pressão Hidrostática --!%%%!===============

		! Matrizes s

	aux1 = gz * dt*dt/(dx*dx) * tetah*tetah
	aux2 = gz * dt*dt/(dy*dy) * tetah*tetah

	do j = 1, ny
	do i = 1, nx1
		matrizsi(i,j) = aux1 * parte1(i,j)
	enddo

	do i = 1, nx
		matrizsj(i,j) = aux2 * parte3(i,j)
	enddo
	enddo

	j = ny1
	do i = 1, nx
		matrizsj(i,j) = aux2 * parte3(i,j)
	enddo

	! Matrizes d, q e normalização do desnível para matriz e
	aux1 = tetah * dt/dx
	aux2 = tetah * dt/dy

	do j = 1, ny
	do i = 1, nx
		matrizd(i,j) = 1 + matrizsi(i+1,j) + matrizsi(i,j) + matrizsj(i,j+1) + matrizsj(i,j)
		matrizq(i,j) = deltah(i,j) - aux1 * (parte2(i+1,j) - parte2(i,j)) - aux2 * (parte4(i,j+1) - parte4(i,j))
	enddo
	enddo

	if (ccx0.eq.0) then  ! Condição periódica
		matrizd(0,:) =   matrizd(nx,:)
		matrizd(nx1,:) = matrizd(1,:)
	else
		matrizd(0,:) =   matrizd(1,:)
		matrizd(nx1,:) = matrizd(nx,:)
	endif

	if (ccy0.eq.0) then  ! Condição periódica
		matrizd(:,0) =   matrizd(:,ny)
		matrizd(:,ny1) = matrizd(:,1)
	else
		matrizd(:,0) =   matrizd(:,1)
		matrizd(:,ny1) = matrizd(:,ny)
	endif



	do j = 0, ny+1
	do i = 0, nx+1
		matrize(i,j) = eta1(i,j) * sqrt(matrizd(i,j))
	enddo
	enddo

	!Normalização das matrizes s e cálculo do primeiro erro
	erro = 0
	alfamu = 0
	do j = 1, ny
	do i = 1, nx
		matrizadipos(i,j) = matrizsi(i+1,j)/sqrt(matrizd(i,j)*matrizd(i+1,j))
		matrizadineg(i,j) = matrizsi(i,j)/sqrt(matrizd(i,j)*matrizd(i-1,j))

		matrizadjpos(i,j) = matrizsj(i,j)/sqrt(matrizd(i,j)*matrizd(i,j+1))
		matrizadjneg(i,j) = matrizsj(i,j+1)/sqrt(matrizd(i,j)*matrizd(i,j-1))


		erro(i,j) = matrize(i,j) - matrizadipos(i,j) * matrize(i+1,j) - matrizadineg(i,j) * matrize(i-1,j) &
		- matrizadjpos(i,j) * matrize(i,j+1) - matrizadjneg(i,j) * matrize(i,j-1) - matrizq(i,j)/sqrt(matrizd(i,j))


		alfamu = alfamu + erro(i,j) * erro(i,j)
	enddo
	enddo

	! Condições de contorno
	if (ccx0.eq.0) then  ! Condição periódica
		erro(0,:)   = erro(nx,:)
		erro(nx1,:) = erro(1,:)
	else
		erro(0,:)   = erro(1,:)
		erro(nx1,:) = erro(nx,:)
	endif

	if (ccy0.eq.0) then  ! Condição periódica
		erro(:,0)   = erro(:,ny)
		erro(:,ny1) = erro(:,1)
	else
		erro(:,0)   = erro(:,1)
		erro(:,ny1) = erro(:,ny)
	endif

	errop = erro

	!%%%%%    início do loop para a redução do erro   %%%%%%!
	do while (alfamu > 0.00000000001)

		!inicialização
		alfa   = 0.
		alfamu = 0.
		alfadi = 0.
		beta   = 0.
		betamu = 0.
		mp     = 0.

		! Parâmetro mp e alfa
		do j = 1, ny
		do i = 1, nx
			mp(i,j) = errop(i,j) - errop(i+1,j)*matrizadipos(i,j) - errop(i-1,j)*matrizadineg(i,j) &
			-errop(i,j+1)*matrizadjpos(i,j)-errop(i,j-1)*matrizadjneg(i,j) 

			alfamu = alfamu + erro(i,j) * erro(i,j)
			alfadi = alfadi + errop(i,j) * mp(i,j)
		enddo
		enddo

		alfa = alfamu / alfadi

		! Recálculo das matrizes e, erro e parâmetro beta
		do j = 1, ny
		do i = 1, nx
			matrize(i,j)  = matrize(i,j)  - alfa * errop(i,j)

			erro(i,j) = erro(i,j) - alfa * mp(i,j)

			betamu = betamu + erro(i,j) * erro(i,j)
		enddo
		enddo

		beta = betamu/alfamu

		! Recálculo de errop
		do j = 1, ny
		do i = 1, nx
			errop(i,j) = erro(i,j) + beta * errop(i,j)
		enddo
		enddo

			! Condições de contorno
		if (ccx0.eq.0) then  ! Condição periódica
			matrize(0,:)   = matrize(nx,:)
			matrize(nx1,:) = matrize(1,:)
			errop(0,:)   = errop(nx,:)
			errop(nx1,:) = errop(1,:)
		else
			matrize(0,:)   = matrize(1,:)
			matrize(nx1,:) = matrize(nx,:)
			errop(0,:)   = errop(1,:)
			errop(nx1,:) = errop(nx,:)
		endif

		if (ccy0.eq.0) then  ! Condição periódica
			matrize(:,0)   = matrize(:,ny)
			matrize(:,ny1) = matrize(:,1)
			errop(:,0)   = errop(:,ny)
			errop(:,ny1) = errop(:,1)
		else
			matrize(:,0)   = matrize(:,1)
			matrize(:,ny1) = matrize(:,ny)
			errop(:,0)   = errop(:,1)
			errop(:,ny1) = errop(:,ny)
		endif



		enddo
		!%%%%%    fim do loop para a redução do erro   %%%%%%!

		! Desnormalização de matriz e para desnível
		do j = 0, ny+1
		do i = 0, nx+1
			eta1(i,j) = matrize(i,j) / sqrt(matrizd(i,j))
		enddo
		enddo

!===============================================================================================================

END SUBROUTINE gradeta

!#############################################################################################################################

SUBROUTINE deltah_()

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA

	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8), dimension(nx1,ny) :: deltahi
	real(8), dimension(nx,ny1) :: deltahj

	!contadores
	integer :: i,j,k

	!auxiliares
	real(8) :: aux1, aux2


	!%%%!-- Cálculo do Deltah --!%%%!
	deltahi = 0.
	deltahj = 0.
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		deltahi(i,j) = deltahi(i,j) + dzx(i,j,k) * u(i,j,k)
	enddo
	enddo
	do j = 1, ny1	
	do i = 1, nx
		deltahj(i,j) = deltahj(i,j) + dzy(i,j,k) * v(i,j,k)
	enddo
	enddo
	enddo

	aux1 = (1.-tetah)*dt/dx
	aux2 = (1.-tetah)*dt/dy

	do j = 1, ny
	do i = 1, nx
		deltah(i,j) = eta0(i,j) - aux1*(deltahi(i+1,j)-deltahi(i,j)) - aux2 *(deltahj(i,j+1)-deltahj(i,j))
	enddo
	enddo


	END SUBROUTINE deltah_


