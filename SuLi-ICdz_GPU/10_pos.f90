SUBROUTINE posdin()

	USE velpre
	USE param
	USE ls_param
	!$ USE omp_lib

	IMPLICIT NONE
	integer :: i, j, k
	real(8),save :: aux1, aux2
	real(8),dimension(nx1,ny,nz) :: rhox, hs_x
	real(8),dimension(nx,ny1,nz) :: rhoy, hs_y
	real(8),dimension(nx,ny,nz1) :: rhoz, hs_z

!PARTE 5	!%%%!-- Correção das Velocidades e do Desnível --!%%%!

	! velocidades corrigidas
	aux1 = dt / dx 
	aux2 = dt / dy 

	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)

#ifdef USE_GPU
        !$omp target map(to: aux1, aux2, rhox, rhoy) map(tofrom: u, v)
#endif
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
#ifdef USE_GPU
        !$omp end target
#endif


	if (t_press.eq.1) then

#ifdef USE_GPU
                !$omp target map(to: rhoz, dzz, prd1) map(tofrom: w)
                !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
                !$omp parallel do collapse(3) private(i,j,k)
#endif
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			w(i,j,k) = w(i,j,k) - dt / dzz(k) * ( prd1(i,j,k)-prd1(i,j,k-1) ) / rhoz(i,j,k)! + dt*tf_pz(i,j,k)/rhoz(i,j,k)

		enddo
		enddo
		enddo
#ifdef USE_GPU
                !$omp end teams distribute parallel do
                !$omp end target
#else
                !$omp end parallel do
#endif

	else

		w(:,:,1) = 0. !velocidade vertical no fundo é zero

#ifdef USE_GPU
                !$omp target map(to: u, v, dx, dy, dz) map(tofrom: w)
                !$omp teams distribute parallel do collapse(3) private(i,j,k)
#else
                !$omp parallel do collapse(3) private(i,j,k)
#endif
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			w(i,j,k+1) = w(i,j,k)-(u(i+1,j,k)-u(i,j,k))*dz(k)/dx - (v(i,j+1,k)-v(i,j,k))*dz(k)/dy
		enddo
		enddo
		enddo
#ifdef USE_GPU
                !$omp end teams distribute parallel do
                !$omp end target
#else
                !$omp end parallel do
#endif

		! anula a velocidade no ar para puramente hidrostatico
		call interpx_cf(hs,nx,ny,nz,hs_x) !(nx1,ny,nz)
		call interpy_cf(hs,nx,ny,nz,hs_y) !(nx,ny1,nz)
		call interpz_cf(hs,nx,ny,nz,dz,hs_z) !(nx,ny,nz1)

#ifdef USE_GPU
                !$omp target map(to: hs_x, hs_y, hs_z) map(from: u, v, w)
#endif
                !$omp parallel

                !$omp do collapse(3)
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

                !$omp do collapse(3)
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

                !$omp do collapse(3)
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
#ifdef USE_GPU
                !$omp end target
#endif

	endif
	prd = prd1
	!prd1 = prd1 + prd0 !arrumar erro da pressão quando prd0 =/ 0


END SUBROUTINE posdin
