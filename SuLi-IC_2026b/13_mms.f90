!Subrotina de verificação do código utilizando soluções manufaturadas

!!! Implementação 15/05/2017
! Leonardo Romero Monteiro


! obs.:1 ao utilizar esta subrotina, deve-se adiconar os termos fontes para as equações de NS e continuidade.

! obs.:2 os casos permanentes só são rodados uma vez, no final do código.




!!!!!!!! PARTE DO CÓDIGO DESATUALIZADA, SEM DZ VARIÁVEL !!!!!!!!!!!


SUBROUTINE mms_i()

	USE velpre
	USE param
	USE cond
	use mms_m
	use ls_param
	!$ USE omp_lib

	IMPLICIT NONE
	!==================================================================================================================
	integer :: i, j, k
	real(8),save :: x,y,z,h,ha, hpi, dpi
	real(8),dimension(nx,ny)     :: h_m
	real(8),dimension(nx1,ny)    :: h_mx
	real(8),dimension(nx,ny1)    :: h_my
	real(8),dimension(nx1,ny,nz) :: lsx
	real(8),dimension(nx,ny1,nz) :: lsy
	real(8),dimension(nx,ny,nz1) :: lsz

	hpi=0.5*pi
	dpi=2.0*pi

	erro_t = 0.

	open (unit=2031989, action= 'write', file= 'dados//erros_mms.txt', status= 'unknown')

	! cálculo das condições iniciais para verificar com as soluções manufaturadas
	if (mms_t == 1) then ! Weng, 2009 - Caso permanente


		!$omp parallel do
                do k = 1, nz
                do j = 1, ny
                do i = 1, nx
                        x = (i-0.5)*dx
			y = (j-0.5)*dy
			z = (k-0.5)*dz

			ls(i,j,k) =  a*dsin(x)*dsin(y)+h0-z
                enddo
                enddo
                enddo
                !$omp end parallel do



                !$omp parallel
                !$omp do
		do j = 1, ny
		do i = 1, nx
		x = (i-0.5)*dx
		y = (j-0.5)*dy
		h_m(i,j) = a*dsin(x)*dsin(y) + h0
		enddo
		enddo
		!$omp end do

                !$omp do
		do j = 1, ny
		do i = 1, nx1
		        x = (i-1. )*dx
		        y = (j-0.5)*dy

		        h_mx(i,j) = a*dsin(x)*dsin(y) + h0
		enddo
		enddo
		!$omp end do

                !$omp do
		do j = 1, ny1
		do i = 1, nx
		        x = (i-0.5)*dx
		        y = (j-1. )*dy

		        h_my(i,j) = a*dsin(x)*dsin(y) + h0
		enddo
		enddo
		!$omp end do
		!$omp end parallel


		!$omp parallel do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny
		y = (j-0.5)*dy
		do i = 1, nx1
		x = (i-1.0)*dx

		u(i,j,k) =  dsin(x)*dcos(y)*dsin(pi*0.5*z/h_mx(i,j)) - dcos(x)*dsin(y)*(dcos(pi*2.*z/h_mx(i,j))-1.)
		enddo
		enddo
		enddo
                !$omp end parallel do



                !$omp parallel do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny1
		y = (j-1.0)*dy
		do i = 1, nx
		x = (i-0.5)*dx

		v(i,j,k) = -dcos(x)*dsin(y)*dsin(pi*0.5*z/h_my(i,j)) + dsin(x)*dcos(y)*(dcos(pi*2.*z/h_my(i,j))-1.)
		enddo
		enddo
		enddo
                !$omp end parallel do


                !$omp parallel do
		do k = 1, nz1
		z = (k-1.0)*dz
		do j = 1, ny
		y = (j-0.5)*dy
		do i = 1, nx
		x = (i-0.5)*dx

		w(i,j,k) = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/h_m(i,j)) &
			     -2.*pi*z/h_m(i,j)*dcos(2.*pi*z/h_m(i,j))+2.*pi)/(2.*pi)
		enddo
		enddo
		enddo
                !$omp end parallel do


	elseif (mms_t == 2) then ! Weng, 2009 - Caso não permanente


		!$omp parallel do
                do k = 1, nz
                do j = 1, ny
                do i = 1, nx
			x = (i-0.5)*dx-pi
			y = (j-0.5)*dy-pi
			z = (k-0.5)*dz

			ls(i,j,k) =  a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0-z
                enddo
                enddo
                enddo
                !$omp end parallel do


                !$omp parallel
                !$omp do
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 h_m(i,j) = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx1
		x = (i-1. )*dx-pi

		 h_mx(i,j) = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do j = 1, ny1
		y = (j-1. )*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 h_my(i,j) = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		enddo
		enddo
		!$omp end do
		!$omp end parallel


                !$omp parallel do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx1
		x = (i-1.0)*dx-pi

		 u(i,j,k) =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/h_mx(i,j)))*dsin(t)
		enddo
		enddo
		enddo
                !$omp end parallel do


                !$omp parallel do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny1
		y = (j-1.0)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 v(i,j,k) =  -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/h_my(i,j)))*dsin(t)
		enddo
		enddo
		enddo
                !$omp end parallel do

                !$omp parallel do
		do k = 1, nz1
		z = (k-1.0)*dz
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 w(i,j,k) =  -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/h_m(i,j)))*dsin(t)
		enddo
		enddo
		enddo
	        !$omp end parallel do

                !$omp parallel do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 prd1(i,j,k) = -coef*dcos(x)*dcos(y)*dcos(hpi*z/h_m(i,j))
		enddo
		enddo
		enddo
                !$omp end parallel do

	endif


END SUBROUTINE mms_i

!==================================================================================================================
!==================================================================================================================

SUBROUTINE mms()

	USE velpre
	USE param
	use mms_m
	use ls_param
	!$ USE omp_lib

	IMPLICIT NONE
	!==================================================================================================================
	integer :: i, j, k, nxc, nyc, nzc, npc

	real(8),dimension(nx1,ny,nz) :: u_m
	real(8),dimension(nx,ny1,nz) :: v_m
	real(8),dimension(nx,ny,nz1) :: w_m
	real(8),dimension(nx,ny,nz)  :: p_m
	real(8),dimension(nx,ny)     :: h_m
	real(8),dimension(nx1,ny)    :: h_mx
	real(8),dimension(nx,ny1)    :: h_my
	real(8),dimension(nx1,ny,nz) :: lsx
	real(8),dimension(nx,ny1,nz) :: lsy
	real(8),dimension(nx,ny,nz1) :: lsz

	real(8),save :: x,y,z,h, hpi

	real(8),save :: erro_u1, erro_v1, erro_w1, erro_p1

	hpi = pi * 0.5

	! cálculo das soluções manufaturadas
	if (mms_t == 1) then ! Weng, 2009 - Caso permanente


                !$omp parallel
                !$omp do
		do j = 1, ny
		y = (j-0.5)*dy
		do i = 1, nx
		x = (i-0.5)*dx

		h_m(i,j) = a*dsin(x)*dsin(y) + h0
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do j = 1, ny
		y = (j-0.5)*dy
		do i = 1, nx1
		x = (i-1. )*dx

		h_mx(i,j) = a*dsin(x)*dsin(y) + h0
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do j = 1, ny1
		y = (j-1. )*dy
		do i = 1, nx
		x = (i-0.5)*dx

		h_my(i,j) = a*dsin(x)*dsin(y) + h0
		enddo
		enddo
		!$omp end do
		!$omp end parallel



                !$omp parallel
                !$omp do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny
		y = (j-0.5)*dy
		do i = 1, nx1
		x = (i-1.0)*dx

		u_m(i,j,k) =  dsin(x)*dcos(y)*dsin(pi*0.5*z/h_mx(i,j)) - dcos(x)*dsin(y)*(dcos(pi*2.*z/h_mx(i,j))-1.)
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny1
		y = (j-1.0)*dy
		do i = 1, nx
		x = (i-0.5)*dx

		v_m(i,j,k) = -dcos(x)*dsin(y)*dsin(pi*0.5*z/h_my(i,j)) + dsin(x)*dcos(y)*(dcos(pi*2.*z/h_my(i,j))-1.)
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz1
		z = (k-1.0)*dz
		do j = 1, ny
		y = (j-0.5)*dy
		do i = 1, nx
		x = (i-0.5)*dx

		w_m(i,j,k) = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/h_m(i,j)) &
			     -2.*pi*z/h_m(i,j)*dcos(2.*pi*z/h_m(i,j))+2.*pi)/(2.*pi)
		enddo
		enddo
		enddo
		!$omp end do
		!$omp end parallel

	elseif (mms_t == 2) then ! Weng, 2009 - Caso não permanente


                !$omp parallel
                !$omp do
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 h_m(i,j) = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx1
		x = (i-1. )*dx-pi

		 h_mx(i,j) = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do j = 1, ny1
		y = (j-1. )*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 h_my(i,j) = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		enddo
		enddo
		!$omp end do
		!$omp end parallel

                !$omp parallel
                !$omp do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx1
		x = (i-1.0)*dx-pi

		 u_m(i,j,k) =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/h_mx(i,j)))*dsin(t)
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny1
		y = (j-1.0)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 v_m(i,j,k) =  -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/h_my(i,j)))*dsin(t)
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz1
		z = (k-1.0)*dz
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 w_m(i,j,k) =  -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/h_m(i,j)))*dsin(t)
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz
		z = (k-0.5)*dz
		do j = 1, ny
		y = (j-0.5)*dy-pi
		do i = 1, nx
		x = (i-0.5)*dx-pi

		 p_m(i,j,k) =  -coef*dcos(x)*dcos(y)*dcos(hpi*z/h_m(i,j))
		enddo
		enddo
		enddo
		!$omp end do
		!$omp end parallel

	endif

		call interpx_cf(ls,nx,ny,nz,lsx) !(nx1,ny,nz)
		call interpy_cf(ls,nx,ny,nz,lsy) !(nx,ny1,nz)
		call interpz_cf(ls,nx,ny,nz,lsz) !(nx,ny,nz1)

	! cálculo dos erros
	erro_u = 0. 
	erro_v = 0.
	erro_w = 0.
	nxc = 0
	nyc = 0
	nzc = 0
	npc = 0


		!$omp parallel
		!$omp do
                do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		if (lsx(i,j,k) > 0.) then
			erro_u(i,j,k) = (u_m(i,j,k)-u(i,j,k))**2.
			nxc = nxc + 1
		endif
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		if (lsy(i,j,k) > 0.) then
			erro_v(i,j,k) = (v_m(i,j,k)-v(i,j,k))**2.
			nyc = nyc + 1
		endif
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		if (lsz(i,j,k) > 0.) then
			erro_w(i,j,k) = (w_m(i,j,k)-w(i,j,k))**2.
			nzc = nzc + 1
		endif
		enddo
		enddo
		enddo
		!$omp end do nowait

                !$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
		if (ls(i,j,k) > 0.) then
			erro_p(i,j,k) = (p_m(i,j,k)-prd1(i,j,k))**2.
			npc = npc + 1
		endif
		enddo
		enddo
		enddo
		!$omp end do
		!$omp end parallel

		erro_u1 = dsqrt(sum(erro_u)/nxc ) 
		erro_v1 = dsqrt(sum(erro_v)/nyc ) 
		erro_w1 = dsqrt(sum(erro_w)/nzc ) 
		erro_p1 = dsqrt(sum(erro_p)/npc ) 

	! plotagem dos erros
	if(dmod(it, ceiling(dt_frame/dt)).eq.0) then
	write(2031989,*) t, dx, erro_u1, dy, erro_v1, dz, erro_w1, 'e_p=', erro_p1
	endif

	erro_t = (erro_u1 + erro_v1 + erro_w1)/3. + erro_t


	!if (it == ts) write(2031989,*) "dt = ", dt, " ",  "erro total acumulado=", erro_t/ts
	if (it == ts) write(*,*) "dt = ", dt, " ",  "erro total acumulado=", erro_t/ts

END SUBROUTINE mms

!==================================================================================================================
!==================================================================================================================

SUBROUTINE mms_bc()

	USE velpre
	USE param
	USE cond
	use mms_m
	use ls_param
	!!$ USE omp_lib

	IMPLICIT NONE
	!==================================================================================================================
	integer :: i,j,k
	real(8),save :: x,y,z,ha,hpi,dpi


	hpi=0.5*pi
	dpi=2.0*pi

	if ((mms_t == 1) .and. (t == 0.)) then ! Weng, 2009 - Caso permanente

                !!$omp parallel

                !!$omp do
		!! c.c. bx
		do k = 0, nz+1
		do j = 0, ny+1
		z = (k-0.5)*dz
		y = (j-0.5)*dy

		i = 0
		x = (i-1.0)*dx
		ha = a*dsin(x)*dsin(y) + h0
		bxx0(j,k)   = dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)

		i = 1
		x = (i-1.0)*dx
		ha = a*dsin(x)*dsin(y) + h0
		bxx1(j,k)   = dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)

		i = nx1
		x = (i-1.0)*dx
		ha = a*dsin(x)*dsin(y) + h0
		bxxf(j,k)   = dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)

		i = nx1+1
		x = (i-1.0)*dx
		ha = a*dsin(x)*dsin(y) + h0
		bxxf1(j,k)  = dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)


		i = 0
		x = (i-.5)*dx
		y = (j-1.)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bxy0(j,k)  = -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)
		x = (i-.5)*dx
		y = (j-.5)*dy
		z = (k-1.)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bxz0(j,k)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			     -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)
		i = nx+1
		x = (i-.5)*dx
		y = (j-1.)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bxyf(j,k)  = -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)
		x = (i-.5)*dx
		y = (j-.5)*dy
		z = (k-1.)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bxzf(j,k)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			     -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)

		enddo
		enddo
                !!$omp end do

                !!$omp do
		!! c.c. by
		do k = 0, nz+1
		do i = 0, nx+1
		z = (k-0.5)*dz
		x = (i-0.5)*dx

		j = 0
		y = (j-1.0)*dy
		ha = a*dsin(x)*dsin(y) + h0
		byy0(i,k)  =  -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)
		j = 1
		y = (j-1.0)*dy
		ha = a*dsin(x)*dsin(y) + h0
		byy1(i,k)  =  -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)
		j = ny1
		y = (j-1.0)*dy
		ha = a*dsin(x)*dsin(y) + h0
		byyf(i,k)  =  -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)
		j = ny1+1
		y = (j-1.0)*dy
		ha = a*dsin(x)*dsin(y) + h0
		byyf1(i,k) =  -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)


		j = 0
		x = (i-1.)*dx
		y = (j-.5)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		byx0(i,k)  = dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)
		x = (i-.5)*dx
		y = (j-.5)*dy
		z = (k-1.)*dz
		ha = a*dsin(x)*dsin(y) + h0
		byz0(i,k)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			     -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)
		j = ny+1
		x = (i-1.)*dx
		y = (j-.5)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		byxf(i,k)  =  dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)
		x = (i-.5)*dx
		y = (j-.5)*dy
		z = (k-1.)*dz
		ha = a*dsin(x)*dsin(y) + h0
		byzf(i,k)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			     -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)
		enddo
		enddo
                !!$omp end do

                !!$omp do
		!! c.c. bz
		do j = 0, ny+1
		do i = 0, nx+1
		y = (j-0.5)*dy
		x = (i-0.5)*dx
		ha = a*dsin(x)*dsin(y) + h0

		k = 0
		z = (k-1.)*dz
		bzz0(i,j)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			     -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)
		k = 1
		z = (k-1.)*dz
		bzz1(i,j)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			     -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)
		k = nz1
		z = (k-1.)*dz
		bzzf(i,j)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			     -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)
		k = nz1+1
		z = (k-1.)*dz
		bzzf1(i,j)  = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(2.*pi*z/ha) &
			      -2.*pi*z/ha*dcos(2.*pi*z/ha)+2.*pi)/(2.*pi)


		k= 0
		x = (i-1.)*dx
		y = (j-.5)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bzx0(i,j)  =  dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)
		x = (i-.5)*dx
		y = (j-1.)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bzy0(i,j)  = -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)

		k = nz+1
		x = (i-1.)*dx
		y = (j-.5)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bzxf(i,j)  =  dsin(x)*dcos(y)*dsin(pi*0.5*z/ha) - dcos(x)*dsin(y)*(dcos(pi*2.*z/ha)-1.)
		x = (i-.5)*dx
		y = (j-1.)*dy
		z = (k-.5)*dz
		ha = a*dsin(x)*dsin(y) + h0
		bzyf(i,j)  = -dcos(x)*dsin(y)*dsin(pi*0.5*z/ha) + dsin(x)*dcos(y)*(dcos(pi*2.*z/ha)-1.)

		enddo
		enddo
		!!$omp end do
		!!$omp end parallel

	elseif (mms_t == 2) then ! Weng, 2009 - Caso não permanente

                !!$omp parallel

                !!$omp do
		!! c.c. bx
		do k = 0, nz+1
		do j = 0, ny+1
		z = (k-0.5)*dz
		y = (j-0.5)*dy-pi

		i = 0
		x = (i-1.0)*dx-pi


		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxx0(j,k)   = dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		i = 1
		x = (i-1.0)*dx-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxx1(j,k)   = dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		i = nx1
		x = (i-1.0)*dx-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxxf(j,k)   = dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		i = nx1+1
		x = (i-1.0)*dx-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxxf1(j,k)  = dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)


		i = 0
		x = (i-.5)*dx-pi
		y = (j-1.)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxy0(j,k)  = -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		x = (i-.5)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-1.)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxz0(j,k)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		i = nx+1
		x = (i-.5)*dx-pi
		y = (j-1.)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxyf(j,k)  = -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		x = (i-.5)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-1.)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bxzf(j,k)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		enddo
		enddo
                !!$omp end do




                !!$omp do
		!! c.c. by
		do k = 0, nz+1
		do i = 0, nx+1
		z = (k-.5)*dz
		x = (i-0.5)*dx-pi

		j = 0
		y = (j-1.0)*dy-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byy0(i,k)  =  -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		j = 1
		y = (j-1.0)*dy-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byy1(i,k)  =  -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		j = ny1
		y = (j-1.0)*dy-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byyf(i,k)  =  -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		j = ny1+1
		y = (j-1.0)*dy-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byyf1(i,k) =  -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)


		j = 0
		x = (i-1.)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byx0(i,k)  = dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		x = (i-.5)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-1.)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byz0(i,k)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		j = ny+1
		x = (i-1.)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byxf(i,k)  =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		x = (i-.5)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-1.)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		byzf(i,k)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		enddo
		enddo
                !!$omp end do





                !!$omp do
		!! c.c. bz
		do j = 0, ny+1
		do i = 0, nx+1
		y = (j-0.5)*dy-pi
		x = (i-0.5)*dx-pi
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0

		k = 0
		z = (k-1.)*dz
		bzz0(i,j)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		k = 1
		z = (k-1.)*dz
		bzz1(i,j)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		k = nz1
		z = (k-1.)*dz
		bzzf(i,j)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		k = nz1+1
		z = (k-1.)*dz
		bzzf1(i,j)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)


		k= 0
		x = (i-1.)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bzx0(i,j)  =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		x = (i-.5)*dx-pi
		y = (j-1.)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bzy0(i,j)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		k = nz+1
		x = (i-1.)*dx-pi
		y = (j-.5)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bzxf(i,j)  =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)
		x = (i-.5)*dx-pi
		y = (j-1.)*dy-pi
		z = (k-.5)*dz
		ha = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		bzyf(i,j)  = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-z/ha))*dsin(t)

		enddo
		enddo
		!!$omp end do
		!!$omp end parallel
	endif


END SUBROUTINE mms_bc

!==================================================================================================================
!==================================================================================================================

SUBROUTINE termo_fonte1()

	! x, y and z are location of nodes in a 3D domain
	! h0 is the reference surface height 
	! a  is the surface wave amplitude
	! sn is the sign of h(x,y): sn=1.0 or -1.0
	! L is the number of vertical mesh and LMAX is the maximum 

	use mms_m
	!!$ USE omp_lib

	IMPLICIT NONE

	integer :: i,j,k

	real(8),save :: x,y,z,hpi,dpi,h1,h2,h3,zh1,h4,h5,dhdx,dhdy,d2hdx2,d2hdy2,d2hdxy, visc

	real(8),save :: ux, vx, wx
	real(8),save :: uy, vy, wy
	real(8),save :: uz, vz, wz

	real(8),save :: dudx,dudy,dudz,d2udx2,d2udy2,d2udz2,d2udxy,d2udxz,d2udyz
	real(8),save :: dvdx,dvdy,dvdz,d2vdx2,d2vdy2,d2vdz2,d2vdxy,d2vdxz,d2vdyz 
	real(8),save :: dwdx,dwdy,dwdz,d2wdx2,d2wdy2,d2wdz2,d2wdxy,d2wdxz,d2wdyz

	hpi=0.5*pi
	dpi=2.0*pi

	visc = 0.00100798/998.
	
	!!$omp parallel
	!-------- derivadas para o termo fonte Fu
	!!$omp do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		x = (i-1.0)*dx
		y = (j-0.5)*dy
		z = (k-0.5)*dz

		dhdx=   a*dcos(x)*dsin(y)
		dhdy=   a*dsin(x)*dcos(y)
		d2hdx2=-a*dsin(x)*dsin(y)	
		d2hdy2=-a*dsin(x)*dsin(y)
		d2hdxy= a*dcos(x)*dcos(y)

		h1 = a*dsin(x)*dsin(y) + h0
		h2 = h1*h1
		h3 = h2*h1
		h4 = h3*h1
		h5 = h4*h1

		zh1 = z/h1

		ux =  dsin(x)*dcos(y)*dsin(hpi*zh1) - dcos(x)*dsin(y)*(dcos(dpi*zh1)-1.)
		vx = -dcos(x)*dsin(y)*dsin(hpi*zh1) + dsin(x)*dcos(y)*(dcos(dpi*zh1)-1.)
		wx = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(dpi*zh1) &
			     -dpi*zh1*dcos(dpi*zh1)+dpi)/(dpi)


		dudx= dcos(x)*dcos(y)*dsin(hpi*zh1)+dsin(x)*dsin(y)*(dcos(dpi*zh1)-1.) &
		     -dsin(x)*dcos(y)*(hpi*z/h2)*dcos(hpi*zh1)*dhdx &
		     -dcos(x)*dsin(y)*(dpi*z/h2)*dsin(dpi*zh1)*dhdx
		
		dudy=-dsin(x)*dsin(y)*dsin(hpi*zh1)-dcos(x)*dcos(y)*(dcos(dpi*zh1)-1.) &
		     +dsin(x)*dcos(y)*(-hpi*z/h2)*dcos(hpi*zh1)*dhdy &
		     +dcos(x)*dsin(y)*(-dpi*z/h2)*dsin(dpi*zh1)*dhdy

		dudz=(hpi/h1)*dsin(x)*dcos(y)*dcos(hpi*zh1)+(dpi/h1)*dcos(x)*dsin(y)*dsin(dpi*zh1)
	  
		d2udx2=-ux &
		       	+2.*dcos(x)*dcos(y)*(-hpi*z/h2)*dcos(hpi*zh1)*dhdx &
			+2.*dsin(x)*dsin(y)*( dpi*z/h2)*dsin(dpi*zh1)*dhdx &
			+dsin(x)*dcos(y)*(-(hpi*z/h2)**2*dsin(hpi*zh1)+ (pi*z/h3)*dcos(hpi*zh1))*dhdx**2 &
			+dcos(x)*dsin(y)*( (dpi*z/h2)**2*dcos(dpi*zh1)+4.*pi*z/h3*dsin(dpi*zh1))*dhdx**2 &
		       -(dsin(x)*dcos(y)*hpi*z/h2*dcos(hpi*zh1)+dcos(x)*dsin(y)*(dpi*z/h2)*dsin(dpi*zh1))*d2hdx2


		d2udy2=-ux &
			-2.*dsin(x)*dsin(y)*(-hpi*z/h2)*dcos(hpi*zh1)*dhdy &
			-2.*dcos(x)*dcos(y)*( dpi*z/h2)*dsin(dpi*zh1)*dhdy &
			+dsin(x)*dcos(y)*(-(hpi*z/h2)**2*dsin(hpi*zh1)+ (pi*z/h3)*dcos(hpi*zh1))*dhdy**2 &
			+dcos(x)*dsin(y)*( (dpi*z/h2)**2*dcos(dpi*zh1)+4.*pi*z/h3*dsin(dpi*zh1))*dhdy**2 &
			-(dsin(x)*dcos(y)*hpi*z/h2*dcos(hpi*zh1)-dcos(x)*dsin(y)*(dpi*z/h2)*dsin(dpi*zh1))*d2hdy2

	       
		d2udz2=-(hpi/h1)**2.*dsin(x)*dcos(y)*dsin(hpi*zh1) &
			+(dpi/h1)**2.*dcos(x)*dsin(y)*dcos(dpi*zh1)
		

		d2vdxy=dsin(x)*dcos(y)*dsin(hpi*zh1)-dcos(x)*dsin(y)*(dcos(dpi*zh1)-1.0) &
			+dsin(x)*dsin(y)*(-hpi*z/h2)*dcos(hpi*zh1)*dhdy &
			+dcos(x)*dcos(y)*( dpi*z/h2)*dsin(dpi*zh1)*dhdy &
			+dcos(x)*dcos(y)*( hpi*z/h2)*dcos(hpi*zh1)*dhdx &
			-dsin(x)*dsin(y)*( dpi*z/h2)*dsin(dpi*zh1)*dhdx &
			+dcos(x)*dsin(y)*( (hpi*z/h2)**2*dsin(hpi*zh1)- (pi*z/h3)*dcos(hpi*zh1))*dhdx*dhdy &
			+dsin(x)*dcos(y)*(-(dpi*z/h2)**2*dcos(dpi*zh1)-(4.*pi*z/h3)*dsin(dpi*zh1))*dhdx*dhdy &
			+(dcos(x)*dsin(y)*(hpi*z/h2)*dcos(hpi*zh1)+dsin(x)*dcos(y)*dpi*z/h2*dsin(dpi*zh1))*d2hdxy

		d2wdxz=-a/dpi*dsin(2.*x)*(dpi*dpi*z/h2)*dsin(dpi*zh1) &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2) &
			*((8.*pi*pi*z/h3)*dsin(dpi*zh1)+(8.*pi*pi*pi*z*z/h4)*dcos(dpi*zh1))*dhdx

		tf_u(i,j,k) = ux*dudx + vx*dudy + wx*dudz &
			     -visc*(d2udx2+d2udy2+d2udz2+d2udx2+d2vdxy+d2wdxz)

	enddo
	enddo
	enddo
	!!$omp end do


	!-------- derivadas para o termo fonte Fv
	!!$omp do
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		x = (i-0.5)*dx
		y = (j-1.0)*dy
		z = (k-0.5)*dz

		dhdx=   a*dcos(x)*dsin(y)
		dhdy=   a*dsin(x)*dcos(y)
		d2hdx2=-a*dsin(x)*dsin(y)	
		d2hdy2=-a*dsin(x)*dsin(y)
		d2hdxy= a*dcos(x)*dcos(y)

		h1 = a*dsin(x)*dsin(y) + h0
		h2 = h1*h1
		h3 = h2*h1
		h4 = h3*h1
		h5 = h4*h1

		zh1 = z/h1

		uy =  dsin(x)*dcos(y)*dsin(hpi*zh1) - dcos(x)*dsin(y)*(dcos(dpi*zh1)-1.)
		vy = -dcos(x)*dsin(y)*dsin(hpi*zh1) + dsin(x)*dcos(y)*(dcos(dpi*zh1)-1.)
		wy = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(dpi*zh1) &
			     -dpi*zh1*dcos(dpi*zh1)+dpi)/(dpi)

		   dvdx= dsin(x)*dsin(y)*dsin(hpi*zh1)+dcos(x)*dcos(y)*(dcos(dpi*zh1)-1.0) &
			+dcos(x)*dsin(y)*(hpi*z/h2)*dcos(hpi*zh1)*dhdx &
			+dsin(x)*dcos(y)*(dpi*z/h2)*dsin(dpi*zh1)*dhdx

		   dvdy=-dcos(x)*dcos(y)*dsin(hpi*zh1)-dsin(x)*dsin(y)*(dcos(dpi*zh1)-1.0) &
			+dcos(x)*dsin(y)*(hpi*z/h2)*dcos(hpi*zh1)*dhdy &
			+dsin(x)*dcos(y)*(dpi*z/h2)*dsin(dpi*zh1)*dhdy

		   dvdz=-(hpi/h1)*dcos(x)*dsin(y)*dcos(hpi*zh1)-dsin(x)*dcos(y)*(dpi/h1)*dsin(dpi*zh1)

		   d2vdx2=-vy &
			+2.*dsin(x)*dsin(y)*(-hpi*z/h2)*dcos(hpi*zh1)*dhdx &
			+2.*dcos(x)*dcos(y)*( dpi*z/h2)*dsin(dpi*zh1)*dhdx &
			-dcos(x)*dsin(y)*(-(hpi*z/h2)**2*dsin(hpi*zh1)+   pi*z/h3*dcos(hpi*zh1))*dhdx**2 &
			-dsin(x)*dcos(y)*( (dpi*z/h2)**2*dcos(dpi*zh1)+4.*pi*z/h3*dsin(dpi*zh1))*dhdx**2 &
			+(dcos(x)*dsin(y)*hpi*z/h2*dcos(hpi*zh1)+dsin(x)*dcos(y)*dpi*z/h2*dsin(dpi*zh1))*d2hdx2

		   d2vdy2=-vy &
			-2.*dcos(x)*dcos(y)*(-hpi*z/h2)*dcos(hpi*zh1)*dhdy &
			-2.*dsin(x)*dsin(y)*( dpi*z/h2)*dsin(dpi*zh1)*dhdy &
			+(dcos(x)*dsin(y)*hpi*z/h2*dcos(hpi*zh1)+dsin(x)*dcos(y)*(dpi*z/h2)*dsin(dpi*zh1))*d2hdy2 &  
			+ dcos(x)*dsin(y)*((hpi*z/h2)**2*dsin(hpi*zh1)-   pi*z/h3*dcos(hpi*zh1))*dhdy*dhdy &
			- dsin(x)*dcos(y)*((dpi*z/h2)**2*dcos(dpi*zh1)+4.*pi*z/h3*dsin(dpi*zh1))*dhdy*dhdy
			
		   d2vdz2=(hpi/h1)**2*dcos(x)*dsin(y)*dsin(hpi*zh1)-dsin(x)*dcos(y)*(dpi/h1)**2*dcos(dpi*zh1)

		   d2udxy=-dcos(x)*dsin(y)*dsin(hpi*zh1)+dsin(x)*dcos(y)*(dcos(dpi*zh1)-1.0) &
			+(dsin(x)*dsin(y)*dhdx-dcos(x)*dcos(y)*dhdy)*(hpi*z/h2)*dcos(hpi*zh1) &
			+(dsin(x)*dsin(y)*dhdy-dcos(x)*dcos(y)*dhdx)*(dpi*z/h2)*dsin(dpi*zh1) &
			-(dsin(x)*dcos(y)*(hpi*z/h2)*dcos(hpi*zh1)+dcos(x)*dsin(y)*(dpi*z/h2)*dsin(dpi*zh1))*d2hdxy &
			+ dsin(x)*dcos(y)*(-(hpi*z/h2)**2*dsin(hpi*zh1)+   pi*z/h3*dcos(hpi*zh1))*dhdx*dhdy &
			+ dcos(x)*dsin(y)*( (dpi*z/h2)**2*dcos(dpi*zh1)+4.*pi*z/h3*dsin(dpi*zh1))*dhdx*dhdy

		   d2wdyz= a/dpi*dsin(2.*y)*(dpi*dpi*z/h2)*dsin(dpi*zh1) &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2) &
			*((8.*pi*pi*z/h3)*dsin(dpi*zh1)+(8.*pi*pi*pi*z*z/h4)*dcos(dpi*zh1))*dhdy	

		tf_v(i,j,k) = uy*dvdx + vy*dvdy + wy*dvdz &
			     -visc*(d2vdx2+d2vdy2+d2vdz2+d2udxy+d2vdy2+d2wdyz) 

	enddo
	enddo
	enddo
	!!$omp end do

	!-------- derivatives of w
	!!$omp do
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		x = (i-0.5)*dx
		y = (j-0.5)*dy
		z = (k-1.0)*dz

		dhdx=   a*dcos(x)*dsin(y)
		dhdy=   a*dsin(x)*dcos(y)
		d2hdx2=-a*dsin(x)*dsin(y)	
		d2hdy2=-a*dsin(x)*dsin(y)
		d2hdxy= a*dcos(x)*dcos(y)

		h1 = a*dsin(x)*dsin(y) + h0
		h2 = h1*h1
		h3 = h2*h1
		h4 = h3*h1
		h5 = h4*h1

		zh1 = z/h1

		uz =  dsin(x)*dcos(y)*dsin(pi*0.5*zh1) - dcos(x)*dsin(y)*(dcos(pi*2.*zh1)-1.)
		vz = -dcos(x)*dsin(y)*dsin(pi*0.5*zh1) + dsin(x)*dcos(y)*(dcos(pi*2.*zh1)-1.)
		wz = -a*(dsin(x)**2.*dcos(y)**2.-dcos(x)**2.*dsin(y)**2.)*(dsin(dpi*zh1) &
			     -dpi*zh1*dcos(dpi*zh1)+dpi)/(dpi)



		   dwdx=-a/dpi*dsin(2.*x)*(dsin(dpi*zh1)-(dpi*zh1)*dcos(dpi*zh1)+dpi) &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2)*(dpi*dpi*z*z/h3)*dsin(dpi*zh1)*dhdx

		   dwdy= a/dpi*dsin(2.*y)*(dsin(dpi*zh1)-(dpi*zh1)*dcos(dpi*zh1)+dpi) &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2)*(dpi*dpi*z*z/h3)*dsin(dpi*zh1)*dhdy

		   dwdz=-a*dpi*z/h2*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2)*dsin(dpi*zh1)

		   d2wdx2=-a/pi*dcos(2.*x)*(dsin(dpi*zh1)-(dpi*zh1)*dcos(dpi*zh1)+dpi) &
			+a/pi*dsin(2.*x)*(dpi*dpi*z*z/h3)*dsin(dpi*zh1)*dhdx &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2)*(-(dpi**3.*z**3./h5) &
			*dcos(dpi*zh1)-(12.*pi*pi*z*z/h4)*dsin(dpi*zh1))*dhdx*dhdx &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2)*(4.*pi*pi*z*z/h3)*dsin(dpi*zh1)*d2hdx2

		   d2wdy2= a/pi*dcos(2.*y)*(dsin(dpi*zh1)-(dpi*zh1)*dcos(dpi*zh1)+dpi) &
			-a/pi*dsin(2.*y)*(dpi*dpi*z*z/h3)*dsin(dpi*zh1)*dhdy &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2)*(-(dpi**3.*z**3./h5) &
			*dcos(dpi*zh1)-(12.*pi*pi*z*z/h4)*dsin(dpi*zh1))*dhdy*dhdy &
			+a/dpi*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2)*(dpi*dpi*z*z/h3)*dsin(dpi*zh1)*d2hdy2  
			
		   d2wdz2=-a*(dsin(x)**2*dcos(y)**2-dcos(x)**2*dsin(y)**2) &
			*(dpi*dpi*z/h3*dcos(dpi*zh1)+(dpi/h2)*dsin(dpi*zh1))
			
		   d2vdyz=-dcos(x)*dcos(y)*( hpi/h1)*dcos(hpi*zh1)+dsin(x)*dsin(y)*(dpi/h1)*dsin(dpi*zh1) &
			+dcos(x)*dsin(y)*(-pi*pi*z/4./h3*dsin(hpi*zh1)+hpi/h2*dcos(hpi*zh1))*dhdy &
			+dsin(x)*dcos(y)*( 4.*pi*pi*z/h3*dcos(dpi*zh1)+dpi/h2*dsin(dpi*zh1))*dhdy

		   d2udxz= dcos(x)*dcos(y)*(hpi/h1)*dcos(hpi*zh1)-dsin(x)*dsin(y)*(dpi/h1)*dsin(dpi*zh1) &
			+dsin(x)*dcos(y)*(pi*pi*z/4./h3*dsin(hpi*zh1)-hpi/h2*dcos(hpi*zh1))*dhdx &
			-dcos(x)*dsin(y)*(4.*pi*pi*z/h3*dcos(dpi*zh1)+dpi/h2*dsin(dpi*zh1))*dhdx

		tf_w(i,j,k) = uz*dwdx + vz*dwdy + wz*dwdz &
			     -visc*(d2wdx2+d2wdy2+d2wdz2+d2udxz+d2vdyz+d2wdz2) 

	enddo
	enddo
	enddo
	!!$omp end do
	!!$omp end parallel


END SUBROUTINE termo_fonte1

!==================================================================================================================
!==================================================================================================================

SUBROUTINE termo_fonte2()


	! x, y and z are location of nodes in a 3D domain
	! h0 is the reference surface height 
	! a  is the surface wave amplitude
	! sn is the sign of h(x,y): sn=1.0 or -1.0
	! L is the number of vertical mesh and LMAX is the maximum 

	use mms_m
	!!$ USE omp_lib

	IMPLICIT NONE

	integer :: i,j,k

	real(8),save :: c,s,x,y,hpi,dpi,h1,h2,h3,zh1,h4,h5,dhdx,dhdy,d2hdx2,d2hdy2,d2hdxy,z,visc

	real(8),save :: dcdx,dcdy,dcdz,dsdx,dsdy,dsdz,d2cdx2,d2cdy2,d2cdz2,d2cdxy,d2cdxz,d2cdyz

	real(8),save :: dpdx,dpdy,dpdz

	real(8),save :: ux, vx, wx, dudt
	real(8),save :: uy, vy, wy, dvdt
	real(8),save :: uz, vz, wz, dwdt

	real(8),save :: dudx,dudy,dudz,d2udx2,d2udy2,d2udz2,d2udxy,d2udxz,d2udyz
	real(8),save :: dvdx,dvdy,dvdz,d2vdx2,d2vdy2,d2vdz2,d2vdxy,d2vdxz,d2vdyz 
	real(8),save :: dwdx,dwdy,dwdz,d2wdx2,d2wdy2,d2wdz2,d2wdxy,d2wdxz,d2wdyz

	hpi=0.5*pi
	dpi=2.0*pi

	visc = 0.00100798/998.

        !!$omp parallel
	!-------- derivadas para o termo fonte Fu
	!!$omp do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		x = (i-1.0)*dx-pi
		y = (j-0.5)*dy-pi
		z = (k-0.5)*dz

		dhdx=   -a*0.5*dsin(x*0.5)*dcos(y*0.5)*dcos(t)
		dhdy=   -a*0.5*dcos(x*0.5)*dsin(y*0.5)*dcos(t)
		d2hdx2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdy2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdxy= -a*0.25*dsin(x*0.5)*dsin(y*0.5)*dcos(t)

		h1 = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		h2 = h1*h1
		h3 = h2*h1

		zh1 = z/h1


		c = dcos(hpi*(1.0-zh1))
		s = dsin(hpi*(1.0-zh1)) 

		dcdx = -hpi*z/h2*s*dhdx
		dcdy = -hpi*z/h2*s*dhdy
		dcdz =  hpi/h1*s

		dsdx =  hpi*z/h2*c*dhdx 	
		dsdy =  hpi*z/h2*c*dhdy	
		dsdz = -hpi/h1*c

		d2cdx2 = -hpi*z/h2*dsdx*dhdx +pi*z/h3*s*dhdx*dhdx -hpi*z/h2*s*d2hdx2
		d2cdy2 = -hpi*z/h2*dsdy*dhdy +pi*z/h3*s*dhdy*dhdy -hpi*z/h2*s*d2hdy2
		d2cdz2 = -(hpi/h1)**2.*c
		d2cdxy = -hpi*z/h2*dsdy*dhdx +pi*z/h3*s*dhdx*dhdy -hpi*z/h2*s*d2hdxy
		d2cdxz =  hpi*z/h2*hpi/h1*c*dhdx-hpi/h2*s*dhdx
		d2cdyz =  hpi*z/h2*hpi/h1*c*dhdy-hpi/h2*s*dhdy

		ux =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-zh1))*dsin(t)
		vx = -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-zh1))*dsin(t)
		wx = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-zh1))*dsin(t)

		dudt = dsin(y)*dcos(0.5*x)*dcos(0.5*x)*c*dcos(t)  

		dudx= -0.5*dsin(x)*dsin(y)*c*dsin(t)+dsin(y)*dcos(0.5*x)*dcos(0.5*x)*dcdx*dsin(t)
		
		dudy= dcos(0.5*x)*dcos(0.5*x)*dcos(y)*c*dsin(t)+dsin(y)*dcos(0.5*x)*dcos(0.5*x)*dcdy*dsin(t)

		dudz= dcos(0.5*x)*dcos(0.5*x)*dsin(y)*dcdz*dsin(t)
	  
		d2udx2= -0.5*dcos(x)*dsin(y)*c*dsin(t)-dsin(x)*dsin(y)*dcdx*dsin(t)+dsin(y)*dcos(0.5*x)*dcos(0.5*x)*d2cdx2*dsin(t)

		d2udy2= -dcos(0.5*x)*dcos(0.5*x)*dsin(y)*c*dsin(t) +2.0*dcos(y)*dcos(0.5*x)*dcos(0.5*x)*dcdy*dsin(t) &
			+dsin(y)*dcos(0.5*x)*dcos(0.5*x)*d2cdy2*dsin(t)

	       	d2udz2= dcos(0.5*x)*dcos(0.5*x)*dsin(y)*d2cdz2*dsin(t)
		
		d2vdxy= 0.5*dsin(y)*dcos(x)*c*dsin(t)+0.5*dsin(x)*dsin(y)*dcdx*dsin(t) &
		       -dcos(0.5*y)*dcos(0.5*y)*dcos(x)*dcdy*dsin(t)-dsin(x)*dcos(0.5*y)*dcos(0.5*y)*d2cdxy*dsin(t)

		d2wdxz= 0.5*a*dsin(0.5*x)*dcos(0.5*y)*dcdz*dsin(t)-a*dcos(0.5*x)*dcos(0.5*y)*d2cdxz*dsin(t)

		tf_u(i,j,k) = dudt + ux*dudx + vx*dudy + wx*dudz &
			     -visc*(d2udx2+d2udy2+d2udz2+d2udx2+d2vdxy+d2wdxz)


		dpdx=-dsin(x)*dcos(y)*dcos(hpi*zh1)+dcos(x)*dcos(y)*dsin(hpi*zh1)*hpi*zh1/h1*dhdx
		tf_px(i,j,k) =  dpdx * coef

	enddo
	enddo
	enddo
	!!$omp end do


	!-------- derivadas para o termo fonte Fv
	!!$omp do
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		x = (i-0.5)*dx-pi
		y = (j-1.0)*dy-pi
		z = (k-0.5)*dz

		dhdx=   -a*0.5*dsin(x*0.5)*dcos(y*0.5)*dcos(t)
		dhdy=   -a*0.5*dcos(x*0.5)*dsin(y*0.5)*dcos(t)
		d2hdx2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdy2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdxy= -a*0.25*dsin(x*0.5)*dsin(y*0.5)*dcos(t)

		h1 = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		h2 = h1*h1
		h3 = h2*h1

		zh1 = z/h1

		c = dcos(hpi*(1.0-zh1))
		s = dsin(hpi*(1.0-zh1)) 

		dcdx = -hpi*z/h2*s*dhdx
		dcdy = -hpi*z/h2*s*dhdy
		dcdz =  hpi/h1*s

		dsdx =  hpi*z/h2*c*dhdx 	
		dsdy =  hpi*z/h2*c*dhdy	
		dsdz = -hpi/h1*c

		d2cdx2 = -hpi*z/h2*dsdx*dhdx +pi*z/h3*s*dhdx*dhdx -hpi*z/h2*s*d2hdx2
		d2cdy2 = -hpi*z/h2*dsdy*dhdy +pi*z/h3*s*dhdy*dhdy -hpi*z/h2*s*d2hdy2
		d2cdz2 = -(hpi/h1)**2.*c
		d2cdxy = -hpi*z/h2*dsdy*dhdx +pi*z/h3*s*dhdx*dhdy -hpi*z/h2*s*d2hdxy
		d2cdxz =  hpi*z/h2*hpi/h1*c*dhdx-hpi/h2*s*dhdx
		d2cdyz =  hpi*z/h2*hpi/h1*c*dhdy-hpi/h2*s*dhdy

		uy =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-zh1))*dsin(t)
		vy = -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-zh1))*dsin(t)
		wy = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-zh1))*dsin(t)

		  dvdt=-dsin(x)*dcos(0.5*y)*dcos(0.5*y)*c*dcos(t)   
	   
		   dvdx= -dcos(x)*dcos(0.5*y)*dcos(0.5*y)*c*dsin(t)-dsin(x)*dcos(0.5*y)*dcos(0.5*y)*dcdx*dsin(t)

		   dvdy= 0.5*dsin(x)*dsin(y)*c*dsin(t)-dsin(x)*dcos(0.5*y)*dcos(0.5*y)*dcdy*dsin(t)

		   dvdz= -dcos(0.5*y)*dcos(0.5*y)*dsin(x)*dcdz*dsin(t)

		   d2vdx2= dsin(x)*dcos(0.5*y)*dcos(0.5*y)*c*dsin(t) -2.0*dcos(x)*dcos(0.5*y)*dcos(0.5*y)*dcdx*dsin(t) &
			  -dsin(x)*dcos(0.5*y)*dcos(0.5*y)*d2cdx2*dsin(t)

		   d2vdy2= 0.5*dsin(x)*dcos(y)*c*dsin(t)+dsin(x)*dsin(y)*dcdy*dsin(t)-dsin(x)*dcos(0.5*y)*dcos(0.5*y)*d2cdy2*dsin(t)
			
		   d2vdz2= -dsin(x)*dcos(0.5*y)*dcos(0.5*y)*d2cdz2*dsin(t)

		   d2udxy= -0.5*dsin(x)*dcos(y)*c*dsin(t)-0.5*dsin(x)*dsin(y)*dcdy*dsin(t) &
			   +dcos(0.5*x)*dcos(0.5*x)*dcos(y)*dcdx*dsin(t)+dsin(y)*dcos(0.5*x)*dcos(0.5*x)*d2cdxy*dsin(t)

		   d2wdyz= 0.5*a*dcos(0.5*x)*dsin(0.5*y)*dcdz*dsin(t)-a*dcos(0.5*x)*dcos(0.5*y)*d2cdyz*dsin(t)	

		tf_v(i,j,k) = dvdt + uy*dvdx + vy*dvdy + wy*dvdz &
			     -visc*(d2vdx2+d2vdy2+d2vdz2+d2udxy+d2vdy2+d2wdyz) 

		dpdy=-dcos(x)*dsin(y)*dcos(hpi*zh1)+dcos(x)*dcos(y)*dsin(hpi*zh1)*hpi*zh1/h1*dhdy
		tf_py(i,j,k) =  dpdy * coef
	enddo
	enddo
	enddo
	!!$omp end do

	!-------- derivatives of w
	!!$omp do
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		x = (i-0.5)*dx-pi
		y = (j-0.5)*dy-pi
		z = (k-1.0)*dz

		dhdx=   -a*0.5*dsin(x*0.5)*dcos(y*0.5)*dcos(t)
		dhdy=   -a*0.5*dcos(x*0.5)*dsin(y*0.5)*dcos(t)
		d2hdx2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdy2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdxy= -a*0.25*dsin(x*0.5)*dsin(y*0.5)*dcos(t)

		h1 = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		h2 = h1*h1
		h3 = h2*h1

		zh1 = z/h1

		c = dcos(hpi*(1.0-zh1))
		s = dsin(hpi*(1.0-zh1)) 

		dcdx = -hpi*z/h2*s*dhdx
		dcdy = -hpi*z/h2*s*dhdy
		dcdz =  hpi/h1*s

		dsdx =  hpi*z/h2*c*dhdx 	
		dsdy =  hpi*z/h2*c*dhdy	
		dsdz = -hpi/h1*c

		d2cdx2 = -hpi*z/h2*dsdx*dhdx +pi*z/h3*s*dhdx*dhdx -hpi*z/h2*s*d2hdx2
		d2cdy2 = -hpi*z/h2*dsdy*dhdy +pi*z/h3*s*dhdy*dhdy -hpi*z/h2*s*d2hdy2
		d2cdz2 = -(hpi/h1)**2.*c
		d2cdxy = -hpi*z/h2*dsdy*dhdx +pi*z/h3*s*dhdx*dhdy -hpi*z/h2*s*d2hdxy
		d2cdxz =  hpi*z/h2*hpi/h1*c*dhdx-hpi/h2*s*dhdx
		d2cdyz =  hpi*z/h2*hpi/h1*c*dhdy-hpi/h2*s*dhdy

		uz =  dsin(y)*dcos(x*0.5)*dcos(x*0.5)*dcos(hpi*(1.-zh1))*dsin(t)
		vz = -dsin(x)*dcos(y*0.5)*dcos(y*0.5)*dcos(hpi*(1.-zh1))*dsin(t)
		wz = -a*dcos(x*0.5)*dcos(y*0.5)*dcos(hpi*(1.-zh1))*dsin(t)

		  dwdt=  -a*dcos(0.5*x)*dcos(0.5*y)*c*dcos(t)

		   dwdx= 0.5*a*dsin(0.5*x)*dcos(0.5*y)*c*dsin(t)-a*dcos(0.5*x)*dcos(0.5*y)*dcdx*dsin(t)

		   dwdy= 0.5*a*dcos(0.5*x)*dsin(0.5*y)*c*dsin(t)-a*dcos(0.5*x)*dcos(0.5*y)*dcdy*dsin(t)

		   dwdz= -a*dcos(0.5*x)*dcos(0.5*y)*dcdz*dsin(t)

		   d2wdx2= 0.25*a*dcos(0.5*x)*dcos(0.5*y)*c*dsin(t)+a*dsin(0.5*x)*dcos(0.5*y)*dcdx*dsin(t) &
			  -a*dcos(0.5*x)*dcos(0.5*y)*d2cdx2*dsin(t)

		   d2wdy2= 0.25*a*dcos(0.5*x)*dcos(0.5*y)*c*dsin(t)+a*dcos(0.5*x)*dsin(0.5*y)*dcdy*dsin(t) &
			  -a*dcos(0.5*x)*dcos(0.5*y)*d2cdy2*dsin(t)
			
		   d2wdz2= -a*dcos(0.5*x)*dcos(0.5*y)*d2cdz2*dsin(t)
			
		   d2vdyz= 0.5*dsin(x)*dsin(y)*dcdz*dsin(t)-dsin(x)*dcos(0.5*y)*dcos(0.5*y)*d2cdyz*dsin(t)

		   d2udxz=-0.5*dsin(x)*dsin(y)*dcdz*dsin(t)+dsin(y)*dcos(0.5*x)*dcos(0.5*x)*d2cdxz*dsin(t)

		tf_w(i,j,k) = dwdt + uz*dwdx + vz*dwdy + wz*dwdz &
			     -visc*(d2wdx2+d2wdy2+d2wdz2+d2udxz+d2vdyz+d2wdz2) 

		dpdz=-dcos(x)*dcos(y)*dsin(hpi*zh1)*hpi/h1
		tf_pz(i,j,k) =  dpdz * coef

	enddo
	enddo
	enddo
	!!$omp end do
	!!$omp end parallel


END SUBROUTINE termo_fonte2


SUBROUTINE termo_fontep2()

	use mms_m
	!!$ USE omp_lib

	IMPLICIT NONE

	integer :: i,j,k

	real(8),save :: x,y,hpi,dpi,h1,h2,h3,zh1,h4,h5,dhdx,dhdy,d2hdx2,d2hdy2,d2hdxy,z

	real(8),save :: d2pdx2,d2pdy2,d2pdz2


	hpi=0.5*pi
	dpi=2.0*pi



	!-------- derivatives of pressure
	!!$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		x = (i-0.5)*dx-pi
		y = (j-0.5)*dy-pi
		z = (k-0.5)*dz

		dhdx=   -a*0.5*dsin(x*0.5)*dcos(y*0.5)*dcos(t)
		dhdy=   -a*0.5*dcos(x*0.5)*dsin(y*0.5)*dcos(t)
		d2hdx2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdy2= -a*0.25*dcos(x*0.5)*dcos(y*0.5)*dcos(t)
		d2hdxy= -a*0.25*dsin(x*0.5)*dsin(y*0.5)*dcos(t)

		h1 = a*dcos(x*0.5)*dcos(y*0.5)*dcos(t) + h0
		h2 = h1*h1

		zh1 = z/h1


		d2pdx2= -dcos(x)*dcos(y)*dcos(hpi*zh1)-2.0*dsin(x)*dcos(y)*dsin(hpi*zh1)*hpi*zh1/h1*dhdx &
			+dcos(x)*dcos(y)*(-dcos(hpi*zh1)*(hpi*zh1/h1)**2*dhdx**2-dsin(hpi*zh1)*pi*zh1/h2*dhdx**2+dsin(hpi*zh1)*hpi*zh1/h1*d2hdx2)

		d2pdy2= -dcos(x)*dcos(y)*dcos(hpi*zh1)-2.0*dcos(x)*dsin(y)*dsin(hpi*zh1)*hpi*zh1/h1*dhdy &
			+dcos(x)*dcos(y)*(-dcos(hpi*zh1)*(hpi*zh1/h1)**2*dhdy**2-dsin(hpi*zh1)*pi*zh1/h2*dhdy**2+dsin(hpi*zh1)*hpi*zh1/h1*d2hdy2)

		d2pdz2= -dcos(x)*dcos(y)*dcos(hpi*zh1)*(hpi/h1)**2




		tf_p(i,j,k) = (d2pdx2 + d2pdy2 + d2pdz2) * coef

	enddo
	enddo
	enddo
	!!$omp end parallel do


END SUBROUTINE termo_fontep2

	! zerando os termos no ar !
!	call interpx_cf(ls,nx,ny,nz,lsx) !(nx1,ny,nz)
!	call interpy_cf(ls,nx,ny,nz,lsy) !(nx,ny1,nz)
!	call interpz_cf(ls,nx,ny,nz,lsz) !(nx,ny,nz1)


!	   do k = 1, nz
!	   do j = 1, ny
!	   do i = 1, nx1
!		if ( lsx(i,j,k) < 0.) then
		! u(i,j,k) = 0.
!		endif
!	   enddo
!	   enddo
!	   enddo

!	   do k = 1, nz
!	   do j = 1, ny
!		if ( lsx(1,j,k) < 0.) then
!		 bxx0(j,k)  = 0.
!		 bxx1(j,k)  = 0.
!		 bxy0(j,k)  = 0.
!		 bxz0(j,k)  = 0.
!		endif

!		if ( lsx(nx1,j,k) < 0.) then
!		 bxxf(j,k)  = 0.
!		 bxxf1(j,k) = 0.
!		 bxyf(j,k)  = 0.
!		 bxzf(j,k)  = 0.
!		endif
!	   enddo
!	   enddo

!	   do k = 1, nz
!	   do j = 1, ny1
!	   do i = 1, nx
!		if ( lsy(i,j,k) < 0.) then
		! v(i,j,k) = 0.
!		endif
!	   enddo
!	   enddo
!	   enddo

!	   do k = 1, nz
!	   do i = 1, nx
!		if ( lsy(i,1,k) < 0.) then
!		 byy0(i,k)  = 0.
!		 byy1(i,k)  = 0.
!		 byx0(i,k)  = 0.
!		 byz0(i,k)  = 0.
!		endif

!		if ( lsy(i,ny1,k) < 0.) then
!		 byyf(i,k)  = 0.
!		 byyf1(i,k) = 0.
!		 byxf(i,k)  = 0.
!		 byzf(i,k)  = 0.
!		endif
!	   enddo
!	   enddo


!	   do k = 1, nz1
!	   do j = 1, ny
!	   do i = 1, nx
!		if ( lsz(i,j,k) < 0. ) then
		! w(i,j,k) = 0.
!		endif
!	   enddo
!	   enddo
!	   enddo

!	   do j = 1, ny
!	   do i = 1, nx
!		if ( lsz(i,j,1) < 0.) then
!		 bzz0(i,j)  = 0.
!		 bzz1(i,j)  = 0.
!		 bzx0(i,j)  = 0.
!		 bzy0(i,j)  = 0.
!		endif

!		if ( lsz(i,j,nz1) < 0.) then
!		 bzzf(i,j)  = 0.
!		 bzzf1(i,j) = 0.
!		 bzxf(i,j)  = 0.
!		 bzyf(i,j)  = 0.
!		endif
!	   enddo
!	   enddo
