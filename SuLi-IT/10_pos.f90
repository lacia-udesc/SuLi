SUBROUTINE poseta(matrizgi,matrizgj,matrizgk,matrizainvi,matrizainvj,matrizainvk)

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE

	real(8), dimension(nx1,ny,nz,nz),intent(in)  :: matrizainvi
	real(8), dimension(nx,ny1,nz,nz),intent(in)  :: matrizainvj
	real(8), dimension(nx,ny,nz1,nz1),intent(in) :: matrizainvk

	real(8), dimension(nx1,ny,nz),intent(in) :: matrizgi
	real(8), dimension(nx,ny1,nz),intent(in) :: matrizgj
	real(8), dimension(nx,ny,nz1),intent(in) :: matrizgk

	integer :: i, j, k, kk
	real(8) :: aux1, aux2

	!Definindo a matriz de resolução de u, v e w
	u = 0.
	v = 0.
	w = 0.

	aux1 = tetah*gz*dt/dx
	aux2 = tetah*gz*dt/dy

	do k = 1, nz
	do kk = 1, nz
	do j = 1, ny
	do i = 1, nx1
		u(i,j,k) = u(i,j,k) + (matrizgi(i,j,kk) - aux1*(eta1(i,j)-eta1(i-1,j))*dzx(i,j,kk))*matrizainvi(i,j,kk,k)
	enddo
	enddo

	do j = 1, ny1
	do i = 1, nx
		v(i,j,k) = v(i,j,k) + (matrizgj(i,j,kk) - aux2*(eta1(i,j)-eta1(i,j-1))*dzy(i,j,kk))*matrizainvj(i,j,kk,k)
	enddo
	enddo
	enddo
	enddo

	do k = 1, nz1
	do kk = 1, nz1
	do i = 1, nx
	do j = 1, ny
		w(i,j,k) = w(i,j,k) + matrizgk(i,j,kk)*matrizainvk(i,j,kk,k)
	enddo
	enddo
	enddo
	enddo


if (mms_t == 2) then


	! Cálculo do desnível corrigido
	dz(:,:,nz) = dz(:,:,nz) + (eta1(:,:) - eta0(:,:))


	! Reajuste dos dzs (só na superfície livre) 
	CALL dz_it()

endif


END SUBROUTINE poseta

!#####################################################################################################################


SUBROUTINE posdin()

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE
	integer :: i, j, k
	real(8) :: aux1, aux2

	! Cálculo do desnível corrigido
	do j = 0, ny+1
	do i = 0, nx+1
		eta1(i,j) = eta1(i,j) + prd1(i,j,nz)/gz
		dz(i,j,nz) = dz(i,j,nz) + (eta1(i,j) - eta0(i,j))
	enddo
	enddo

	! Reajuste dos dzs (só na superfície livre) 
	CALL dz_it()

!===============================================================================================================
!PARTE 5	!%%%!-- Correção das Velocidades e do Desnível --!%%%!

	! velocidades corrigidas

	aux1 = tetah * dt / dx 
	aux2 = tetah * dt / dy 
	
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		u(i,j,k) = u(i,j,k) -aux1 * ( prd1(i,j,k)-prd1(i-1,j,k) ) 
	enddo
	enddo

	do j = 1, ny1
	do i = 1, nx
		v(i,j,k) = v(i,j,k) -aux2 * ( prd1(i,j,k)-prd1(i,j-1,k) )
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		w(i,j,k) = w(i,j,k) - tetah*dt/dzz(i,j,k) * ( prd1(i,j,k)-prd1(i,j,k-1) )
	enddo
	enddo
	enddo

	! pela conservação de massa
	!w(:,:,1) = 0. !velocidade vertical no fundo é zero
	!do k = 1, nz
	!do j = 1, ny
	!do i = 1, nx
	!	w(i,j,k+1) = w(i,j,k)-( u(i+1,j,k)*dzx(i+1,j,k)-u(i,j,k)*dzx(i,j,k) )/dx - ( v(i,j+1,k)*dzy(i,j+1,k)-v(i,j,k)*dzy(i,j,k) )/dy
	!enddo
	!enddo
	!enddo


	prd = prd1
	!prd1 = prd1 + prd0 !arrumar erro da pressão quando prd0 =/ 0


END SUBROUTINE posdin
