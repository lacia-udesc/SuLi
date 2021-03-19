!Subrotina definir as os diferentes dzs

!!! Implementação 12/12/2014
! Leonardo Romero Monteiro

!!! Modificações
! Leonardo Romero Monteiro - 13/01/2015

SUBROUTINE dz_i()

	USE dzs
        USE cond

	IMPLICIT NONE

	integer :: i, j, k

	do k = 0, nz+1
	do j = 0, ny+1
	do i = 1, nx1
		dzx(i,j,k) = (dz(i,j,k) + dz(i-1,j,k))*0.5
	enddo

	if (ccx0.eq.0) then  ! Condição periódica
		i = 0
		dzx(i,j,k) =  (dz(i,j,k) + dz(nx,j,k))*0.5 
		i = nx1+1
		dzx(i,j,k) =  (dz(1,j,k) + dz(i-1,j,k))*0.5  
	else
		i = 0
		dzx(i,j,k) = dzx(i+1,j,k) + (dzx(i+1,j,k) - dzx(i+2,j,k))
		i = nx1+1
		dzx(i,j,k) = dzx(i-1,j,k) + (dzx(i-1,j,k) - dzx(i-2,j,k))
	endif
	enddo
	enddo


	do k = 0, nz+1
	do i = 0, nx+1
	do j = 1, ny1
		dzy(i,j,k) = (dz(i,j,k) + dz(i,j-1,k))*0.5
	enddo



	if (ccy0.eq.0) then  ! Condição periódica
		j = 0
		dzy(i,j,k) = (dz(i,j,k) + dz(i,ny,k))*0.5 
		j = ny1+1
		dzy(i,j,k) = (dz(i,1,k) + dz(i,j-1,k))*0.5 
	else

		j = 0
		dzy(i,j,k) = dzy(i,j+1,k) + (dzy(i,j+1,k) - dzy(i,j+2,k))
		j = ny1+1
		dzy(i,j,k) = dzy(i,j-1,k) + (dzy(i,j-1,k) - dzy(i,j-2,k))
	endif

	enddo
	enddo


	do j = 0, ny+1
	do i = 0, nx+1
	do k = 1, nz1
		dzz(i,j,k) = (dz(i,j,k) + dz(i,j,k-1))*0.5
	enddo
		k = 0
		dzz(i,j,k) = dzz(i,j,k+1) + (dzz(i,j,k+1) - dzz(i,j,k+2))
		k = nz1+1
		dzz(i,j,k) = dzz(i,j,k-1) + (dzz(i,j,k-1) - dzz(i,j,k-2))
	enddo
	enddo
	!===============================================================================================================

END SUBROUTINE dz_i


SUBROUTINE dz_it()

	USE dzs
        USE cond

	IMPLICIT NONE

	integer :: i, j, k

	do k = nz, nz+1
	do j = 0, ny+1
	do i = 1, nx1
		dzx(i,j,k) = (dz(i,j,k) + dz(i-1,j,k))*0.5
	enddo
	if (ccx0.eq.0) then  ! Condição periódica
		i = 0
		dzx(i,j,k) =  (dz(i,j,k) + dz(nx,j,k))*0.5 
		i = nx1+1
		dzx(i,j,k) =  (dz(1,j,k) + dz(i-1,j,k))*0.5  
	else
		i = 0
		dzx(i,j,k) = dzx(i+1,j,k) + (dzx(i+1,j,k) - dzx(i+2,j,k))
		i = nx1+1
		dzx(i,j,k) = dzx(i-1,j,k) + (dzx(i-1,j,k) - dzx(i-2,j,k))
	endif
	enddo
	enddo

	do k = nz, nz+1
	do i = 0, nx+1
	do j = 1, ny1
		dzy(i,j,k) = (dz(i,j,k) + dz(i,j-1,k))*0.5
	enddo
	if (ccy0.eq.0) then  ! Condição periódica
		j = 0
		dzy(i,j,k) = (dz(i,j,k) + dz(i,ny,k))*0.5 
		j = ny1+1
		dzy(i,j,k) = (dz(i,1,k) + dz(i,j-1,k))*0.5 
	else
		j = 0
		dzy(i,j,k) = dzy(i,j+1,k) + (dzy(i,j+1,k) - dzy(i,j+2,k))
		j = ny1+1
		dzy(i,j,k) = dzy(i,j-1,k) + (dzy(i,j-1,k) - dzy(i,j-2,k))
	endif
	enddo
	enddo

	do j = 0, ny+1
	do i = 0, nx+1
		k = nz
		dzz(i,j,k) = (dz(i,j,k) + dz(i,j,k-1))*0.5

		k = nz1
		dzz(i,j,k) = (dz(i,j,k) + dz(i,j,k-1))*0.5

		k = nz1+1
		dzz(i,j,k) = dzz(i,j,k-1) + (dzz(i,j,k-1) - dzz(i,j,k-2))
	enddo
	enddo


END SUBROUTINE dz_it


!===============================================================================================================
!===============================================================================================================
!===============================================================================================================

SUBROUTINE zs(zx,zy,zz)

USE dzs

IMPLICIT NONE
!==================================================================================================================
integer :: i, j, k


real(8), dimension(0:nx1+1,0:ny+1,0:nz+1)   :: zx
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1)   :: zy
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1)   :: zz


	zx(:,:,0) = -dz1*0.5
	zy(:,:,0) = -dz1*0.5
	zz(:,:,0) = -dz1
	do k = 1, nz+1
	do j = 0, ny+1
	do i = 0, nx+1
		zx(i,j,k) = (dzx(i,j,k)+dzx(i,j,k-1))*0.5 + zx(i,j,k-1)
		zy(i,j,k) = (dzy(i,j,k)+dzy(i,j,k-1))*0.5 + zy(i,j,k-1)
		zz(i,j,k) = zz(i,j,k-1) + dz(i,j,k-1)
	enddo
	enddo
	enddo

	i = nx1+1
	do k = 1, nz+1
	do j = 0, ny+1
		zx(i,j,k) = (dzx(i,j,k)+dzx(i,j,k-1))*0.5 + zx(i,j,k-1)
	enddo
	enddo

	j = ny1+1
	do k = 1, nz+1
	do i = 0, nx+1
		zy(i,j,k) = (dzy(i,j,k)+dzy(i,j,k-1))*0.5 + zy(i,j,k-1)
	enddo
	enddo

	k = nz1+1
	do j = 0, ny+1
	do i = 0, nx+1
		zz(i,j,k) = zz(i,j,k-1) + dz(i,j,k-1)
	enddo
	enddo
	
END SUBROUTINE zs


!#####################################################################################################################
!#####################################################################################################################
!subrotinas criadas para encontrar a declividade ideal quando se possui o número de manning e a velocidade média.

SUBROUTINE analise_i()

USE parametros
USE analis
USE disc
	IMPLICIT NONE

	!!adaptação da declividade
	decliv1 = decliv
	ut_m0 = uinicial

END SUBROUTINE analise_i

!#####################################################################################################################

SUBROUTINE analise_f()

USE parametros
USE analis
	USE dzs
	USE velpre
	IMPLICIT NONE

real(8) :: ud

ut_m = sum(u(1:nx1,1:ny,1:nz)) / (nx1*ny*nz)

ud = ut_m - ut_m0


if (ut_m > uinicial + 0.00000001 .and. ud > -0.0001) then
decliv1 = decliv1 - 0.000000001
elseif (ut_m < uinicial - 0.00000001 .and. ud < 0.0001) then
decliv1 = decliv1 + 0.000000001
endif

gx = 9.80665 * sin(atan(decliv1))

gz = 9.80665 * cos(atan(decliv1))

ut_m0 = ut_m

if(mod(it, ceiling(dt_frame/dt)).eq.0) then
 write(*,*) "declividade",decliv1
endif
END SUBROUTINE analise_f


!#####################################################################################################################
!#####################################################################################################################
!Subrotinas para calcular matrizes A**-1

!!! Implementação 15/04/2014
! Leonardo Romero Monteiro

!!! Modificações
! Leonardo Romero Monteiro - 13/01/2015

SUBROUTINE matrizainv(matrizainvi,matrizainvj,matrizainvk)

	USE dzs
	USE parametros
	USE obst
	USE smag

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA

	real(8), intent(out), dimension(nx1,ny,nz,nz)  :: matrizainvi
	real(8), intent(out), dimension(nx,ny1,nz,nz)  :: matrizainvj
	real(8), intent(out), dimension(nx,ny,nz1,nz1) :: matrizainvk

	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)
	
	real(8), dimension(nx1,ny,nz1) :: dzzx, xznut
	real(8), dimension(nx,ny1,nz1) :: dzzy, yznut
	real(8), dimension(nz,3)       :: matrizi, matrizj
	real(8), dimension(nz1,3)      :: matrizk

	real(8), dimension(nz,nz)      :: invmatrizi, invmatrizj
	real(8), dimension(nz1,nz1)    :: invmatrizk

	!contadores
	integer :: i, j, k

	!auxiliares
	real(8) :: aux1

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA

	! Definição das matrizes zzs
	do k = 1, nz1
	do i = 1, nx1
	do j = 1, ny
		dzzx(i,j,k) = (dz(i,j,k) + dz(i-1,j,k) + dz(i,j,k-1) + dz(i-1,j,k-1))*0.25
	enddo
	enddo

	do j = 1, ny1
	do i = 1, nx
		dzzy(i,j,k) = (dz(i,j,k) + dz(i,j-1,k) + dz(i,j,k-1) + dz(i,j-1,k-1))*0.25
	enddo
	enddo
	enddo

	call interpz_cf(xnut,nx1,ny,nz,xznut) !(nx1,ny,nz1)
	call interpz_cf(ynut,nx,ny1,nz,yznut) !(nx,ny1,nz1)

	!%%%!--Matrizes A (Matrix)-!%%%!
	! Matriz Ai**-1

	do j = 1, ny
	do i = 1, nx1
	! Reinício destas variáveis para cada variação horizontal
		matrizi = 0.
		invmatrizi = 0.

		matrizi(1,2) = - xznut(i,j,2) * dt / dzzx(i,j,2)
		matrizi(1,3) = matrizi(1,2)
		matrizi(1,1) = dzx(i,j,1) - matrizi(1,2) 

		do k = 2, (nz-1)
			matrizi(k,2) = - xznut(i,j,k+1) * dt / dzzx(i,j,k+1)
			matrizi(k,3) = matrizi(k,2)
			matrizi(k,1) = dzx(i,j,k) - matrizi(k,2) - matrizi(k-1,3)
		enddo 

		matrizi(nz,1) = dzx(i,j,nz) - matrizi(nz-1,3)
		
		! Inversão da matriz A
	    	CALL findinv(matrizi, invmatrizi, nz) 

		! Gravar a matriz A**-1 produzida para as variações horizontais
		matrizainvi (i,j,:,:) = invmatrizi (:,:)
	enddo
	enddo


	! Matriz Aj**-1
	do j = 1, ny1
	do i = 1, nx
		matrizj = 0.
		invmatrizj = 0.

		matrizj(1,2) = - yznut(i,j,2) * dt / dzzy(i,j,2)
		matrizj(1,3) = matrizj(1,2)
		matrizj(1,1) = dzy(i,j,1) - matrizj(1,2)

		do k = 2, (nz-1)
			matrizj(k,2) = - yznut(i,j,k+1) * dt / dzzy(i,j,k+1)
			matrizj(k,3) = matrizj(k,2)
			matrizj(k,1) = dzy(i,j,k) - matrizj(k,2) - matrizj(k-1,3)
		enddo 

		matrizj(nz,1) = dzy(i,j,nz) - matrizj(nz-1,3)

		! Inversão da matriz A
	    	 CALL findinv(matrizj, invmatrizj, nz) 

		! Gravar a matriz A**-1 produzida para as variações horizontais
		matrizainvj (i,j,:,:) = invmatrizj (:,:)
	enddo
	enddo



	! Matriz Ak**-1 !! É invertido o nz1 do resto por causa do posicionamento!
	do j = 1, ny
	do i = 1, nx
		matrizk = 0.
		invmatrizk = 0.

		matrizk(1,2) = - nut(i,j,1) * dt / dz(i,j,1)
		matrizk(1,3) = matrizk(1,2)
		matrizk(1,1) = dzz(i,j,1) - matrizk(1,2) !- matrizk(1,2)

		do k = 2, (nz1-1)
			matrizk(k,2) = - nut(i,j,k) * dt / dz(i,j,k)
			matrizk(k,3) = matrizk(k,2)
			matrizk(k,1) = dzz(i,j,k) - matrizk(k,2) - matrizk(k-1,3)
		enddo 

		matrizk(nz1,1) = dzz(i,j,nz1) - matrizk(nz1-1,3)

	    	 CALL findinv(matrizk, invmatrizk, nz1) 

		matrizainvk (i,j,:,:) = invmatrizk (:,:)
	enddo
	enddo


	! comentários
	!! nos contornos os "matrizk" estão x2 para considerar as condições de contorno virtuais

	END SUBROUTINE matrizainv



!########################################################################################

SUBROUTINE findinv(matriz, invmatriz, nz)
! algoritmo de inversão de matrizes tridiagonais !
	implicit none
	!Declarations
	integer, intent(in) :: nz
	real(8), intent(in),  dimension(nz,3) :: matriz
	real(8), intent(out), dimension(nz,nz) :: invmatriz

	integer :: k, kk, j
	real(8), dimension(0:nz) :: tetha
	real(8), dimension(nz+1) :: phi
	real(8) :: multmat

        multmat = 1.

	tetha(0) = 1.
	tetha(1) = matriz(1,1)
	do k = 2, nz
		tetha(k) = matriz(k,1)*tetha(k-1) - matriz(k-1,2)*matriz(k-1,3)*tetha(k-2)
	enddo

	phi(nz+1) = 1.
	phi(nz) = matriz(nz,1)
	do k = nz-1, 1, -1
		phi(k) = matriz(k,1)*phi(k+1) - matriz(k,2)*matriz(k,3)*phi(k+2)
	enddo


	! para i < j, triangulo superior
	do k = 1, nz ! i
	do kk = k+1, nz ! j
		do j = k, kk-1
		multmat = matriz(j,2) * multmat
		enddo
	   invmatriz(k,kk) = (-1.)**(k+kk) * multmat * tetha(k-1) * phi(kk+1) / tetha(nz)
           multmat = 1.
	enddo
	enddo

	! para i > j, triangulo inferior
	do k = 1, nz
	do kk = k-1, 1, -1
		do j = kk, k-1
		multmat = matriz(j,3) * multmat
		enddo
	   invmatriz(k,kk) = (-1.)**(k+kk) * multmat * tetha(kk-1) * phi(k+1) / tetha(nz)
           multmat = 1.
	enddo
	enddo

	! para j = i, centro
	do k = 1, nz
	   invmatriz(k,k) = tetha(k-1) * phi(k+1) / tetha(nz)
	enddo




END SUBROUTINE findinv

