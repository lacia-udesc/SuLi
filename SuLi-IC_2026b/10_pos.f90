SUBROUTINE posdin()

	USE velpre
	USE param
	USE ls_param
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	real(8) :: aux1, aux2
	real(8),dimension(nx1,ny,nz) :: rhox, hs_x
	real(8),dimension(nx,ny1,nz) :: rhoy, hs_y
	real(8),dimension(nx,ny,nz1) :: rhoz, hs_z

!PARTE 5	!%%%!-- Correção das Velocidades e do Desnível --!%%%!

	! velocidades corrigidas
	aux1 = dt * inv_dx
	aux2 = dt * inv_dy

	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)


    !$omp parallel
    !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		u(i,j,k) = u(i,j,k) -aux1 * ( prd1(i,j,k)-prd1(i-1,j,k) ) / rhox(i,j,k)! + dt*tf_px(i,j,k)/rhox(i,j,k)
	enddo
	enddo
	enddo
    !$omp end do nowait

    !$omp do collapse(3)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		v(i,j,k) = v(i,j,k) -aux2 * ( prd1(i,j,k)-prd1(i,j-1,k) ) / rhoy(i,j,k)! + dt*tf_py(i,j,k)/rhoy(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel

	if (t_press.eq.1) then

        !$omp parallel do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			w(i,j,k) = w(i,j,k) - dt * inv_dzz(k) * ( prd1(i,j,k)-prd1(i,j,k-1) ) / rhoz(i,j,k)! + dt*tf_pz(i,j,k)/rhoz(i,j,k)

		enddo
		enddo
		enddo
        !$omp end parallel do

	else

		w(:,:,1) = 0. !velocidade vertical no fundo é zero

		! Recorrência sequencial em k (w(i,j,k+1) depende de w(i,j,k)):
		! paraleliza por coluna (j,i) e mantém k seriado dentro de cada coluna.
        !$omp parallel
        !$omp do
		do j = 1, ny
		do i = 1, nx
			do k = 1, nz
				w(i,j,k+1) = w(i,j,k)-(u(i+1,j,k)-u(i,j,k))*dz(k)*inv_dx - (v(i,j+1,k)-v(i,j,k))*dz(k)*inv_dy
			enddo
		enddo
		enddo
        !$omp end do
        !$omp end parallel

		! anula a velocidade no ar para puramente hidrostatico
		call interpx_cf(hs,nx,ny,nz,hs_x) !(nx1,ny,nz)
		call interpy_cf(hs,nx,ny,nz,hs_y) !(nx,ny1,nz)
		call interpz_cf(hs,nx,ny,nz,dz,hs_z) !(nx,ny,nz1)


        !$omp parallel
        !$omp do 
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		 if (hs_x(i,j,k) < 0.0001) then
		  u(i,j,k) = 0.
		 endif
		enddo
		enddo
		enddo
		!$omp end do nowait

        !$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		 if (hs_y(i,j,k) < 0.0001) then
		  v(i,j,k) = 0.
		 endif
		enddo
		enddo
		enddo
		!$omp end do nowait

        !$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		 if (hs_z(i,j,k) < 0.0001) then
		  w(i,j,k) = 0.
		 endif
		enddo
		enddo
		enddo
        !$omp end do
        !$omp end parallel


	endif
	prd = prd1
	!prd1 = prd1 + prd0 !arrumar erro da pressão quando prd0 =/ 0


END SUBROUTINE posdin
