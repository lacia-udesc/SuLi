SUBROUTINE posdin()

	USE velpre
	USE parametros
	USE mms_m
	USE ls_param
	IMPLICIT NONE
	integer :: i, j, k
	real(8) :: aux1, aux2
	real(8), dimension(nx1,ny,nz) :: rhox, hs_x
	real(8), dimension(nx,ny1,nz) :: rhoy, hs_y
	real(8), dimension(nx,ny,nz1) :: rhoz, hs_z

!===============================================================================================================
!PARTE 5	!%%%!-- Correção das Velocidades e do Desnível --!%%%!

	! velocidades corrigidas

	aux1 = dt / dx 
	aux2 = dt / dy 

	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,rhoz) !(nx,ny,nz1)

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		u(i,j,k) = u(i,j,k) -aux1 * ( prd1(i,j,k)-prd1(i-1,j,k) ) / rhox(i,j,k)! + dt*tf_px(i,j,k)/rhox(i,j,k)
	enddo
	enddo

	do j = 1, ny1
	do i = 1, nx
		v(i,j,k) = v(i,j,k) -aux2 * ( prd1(i,j,k)-prd1(i,j-1,k) ) / rhoy(i,j,k)! + dt*tf_py(i,j,k)/rhoy(i,j,k)
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		w(i,j,k) = w(i,j,k) - dt / dz * ( prd1(i,j,k)-prd1(i,j,k-1) ) / rhoz(i,j,k)! + dt*tf_pz(i,j,k)/rhoz(i,j,k)

	enddo
	enddo
	enddo

!	w(:,:,1) = 0. !velocidade vertical no fundo é zero
!	do k = 1, nz
!	do j = 1, ny
!	do i = 1, nx
!		w(i,j,k+1) = w(i,j,k)-(u(i+1,j,k)-u(i,j,k))*dz/dx - (v(i,j+1,k)-v(i,j,k))*dz/dy
!	enddo
!	enddo
!	enddo

!	call interpx_cf(hs,nx,ny,nz,hs_x) !(nx1,ny,nz)
!	call interpy_cf(hs,nx,ny,nz,hs_y) !(nx,ny1,nz)
!	call interpz_cf(hs,nx,ny,nz,hs_z) !(nx,ny,nz1)


!	do k = 1, nz
!	do j = 1, ny
!	do i = 1, nx1
!	 if (hs_x(i,j,k) < 0.0001) then
!	  u(i,j,k) = 0.
!	 endif
!	enddo
!	enddo
!	enddo

!	do k = 1, nz
!	do j = 1, ny1
!	do i = 1, nx
!	 if (hs_y(i,j,k) < 0.0001) then
!	  v(i,j,k) = 0.
!	 endif
!	enddo
!	enddo
!	enddo

!	do k = 1, nz1
!	do j = 1, ny
!	do i = 1, nx
!	 if (hs_z(i,j,k) < 0.0001) then
!	  w(i,j,k) = 0.
!	 endif
!	enddo
!	enddo
!	enddo

	prd = prd1
	!prd1 = prd1 + prd0 !arrumar erro da pressão quando prd0 =/ 0


END SUBROUTINE posdin
