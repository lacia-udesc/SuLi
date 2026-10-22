!Subrotina visco --> viscosidade por modelo de turbulência

!!! Implementação 19/07/16
! Luísa Vieira Lucchese

!!! Modificação 
! Leonardo Romero Monteiro em 18/09/16
! Luísa V. Lucchese em 31/12/16
! Bruna Fernanda Soares em 15/08/23

SUBROUTINE visco()

	USE param
	USE tempo_hora_data
	USE les
	USE velpre
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k

	real(8) :: p1, p2, p3, lcar, sij
	real(8) :: ydist, nu_loc, k_loc, omega_loc, Pk, gamma_loc, beta_loc, arg1

	real(8),dimension(nx,ny,nz) :: Fka, diss, des, e0, fe0, fe1
	real(8),dimension(nx,ny,nz) :: Fom, e0o, fe0o, fe1o
	real(8),dimension(nx,ny,nz) :: dkdx, dkdy, dkdz, dodx, dody, dodz
	real(8),dimension(nx,ny,nz) :: S2, Smag, F1, F2
	real(8),dimension(nx,ny,nz) :: sigk, sigom

	real(8),dimension(nx1,ny,nz) :: dudx_i, dudy_i, dudz_i
	real(8),dimension(nx,ny1,nz) :: dvdx_i, dvdy_i, dvdz_i
	real(8),dimension(nx,ny,nz1) :: dwdx_i, dwdy_i, dwdz_i

	real(8),dimension(nx,ny,nz) :: dudx, dudy, dudz
    real(8),dimension(nx,ny,nz) :: dvdx, dvdy, dvdz
    real(8),dimension(nx,ny,nz) :: dwdx, dwdy, dwdz

	real(8),dimension(nx1,ny,nz) :: dudx_x, dudy_x, dudz_x, dwdx_x, dvdx_x
	real(8),dimension(nx,ny1,nz) :: dudy_y, dvdx_y, dwdy_y, dvdz_y, dvdy_y
	real(8),dimension(nx,ny,nz1) :: dwdy_z, dvdz_z, dwdz_z, dwdx_z, dudz_z
	real(8),dimension(nx1,ny1,nz) :: dvdx_a, dudy_a
	real(8),dimension(nx1,ny,nz1) :: dwdx_a, dudz_a 
   	real(8),dimension(nx,ny1,nz1) :: dwdy_a, dvdz_a

	real(8),parameter :: zero = 0.0d0
	real(8),parameter :: one  = 1.0d0
	real(8),parameter :: tiny = 1.0d-14
	real(8),parameter :: omega_min = 1.0d-6

	! Constantes SST (Menter)
	real(8),parameter :: betaStar = 0.09d0
	real(8),parameter :: a1_sst = 0.31d0
	real(8),parameter :: sigma_k1 = 0.85d0, sigma_om1 = 0.5d0
	real(8),parameter :: sigma_k2 = 1.0d0,  sigma_om2 = 0.856d0
	real(8),parameter :: beta1 = 0.075d0, beta2 = 0.0828d0
	real(8),parameter :: gamma1 = 5.d0/9.d0, gamma2 = 0.44d0

	!===================================================================================================================
	if (m_turb == 0) then ! 0 = sem modelo de turbulência

		mut = 0.
		ka  = 0.

	else

		!computando as derivadas
		call derivax(u,nx1,ny,nz,dudx_i); call derivax(v,nx,ny1,nz,dvdx_i); call derivax(w,nx,ny,nz1,dwdx_i)
		call derivay(u,nx1,ny,nz,dudy_i); call derivay(v,nx,ny1,nz,dvdy_i); call derivay(w,nx,ny,nz1,dwdy_i)
  		call derivaz(u,nx1,ny,nz,inv_2dz,dudz_i); call derivaz(v,nx,ny1,nz,inv_2dz,dvdz_i); call derivaz(w,nx,ny,nz1,inv_2dzz,dwdz_i)

 		! as interpolações abaixo são interpolando tudo para o centro da celula

		! interpolação para u
		call interpx_fc(dudx_i,nx1,ny,nz,dudx)!dudx(i,j,k)=(dudx_i(i+1,j,k)+dudx_i(i,j,k))*0.5
		call interpx_fc(dudy_i,nx1,ny,nz,dudy)!dudy(i,j,k)=(dudy_i(i+1,j,k)+dudy_i(i,j,k))*0.5
		call interpx_fc(dudz_i,nx1,ny,nz,dudz)!dudz(i,j,k)=(dudz_i(i+1,j,k)+dudz_i(i,j,k))*0.5

		! interpolação para v
		call interpy_fc(dvdx_i,nx,ny1,nz,dvdx)!dvdx(i,j,k)=(dvdx_i(i,j+1,k)+dvdx_i(i,j,k))*0.5
		call interpy_fc(dvdy_i,nx,ny1,nz,dvdy)!dvdy(i,j,k)=(dvdy_i(i,j+1,k)+dvdy_i(i,j,k))*0.5
		call interpy_fc(dvdz_i,nx,ny1,nz,dvdz)!dvdz(i,j,k)=(dvdz_i(i,j+1,k)+dvdz_i(i,j,k))*0.5

		! interpolação para w
		call interpz_fc(dwdx_i,nx,ny,nz1,dzz,dwdx)!dwdx(i,j,k)=(dwdx_i(i,j,k+1)+dwdx_i(i,j,k))*0.5
		call interpz_fc(dwdy_i,nx,ny,nz1,dzz,dwdy)!dwdy(i,j,k)=(dwdy_i(i,j,k+1)+dwdy_i(i,j,k))*0.5
		call interpz_fc(dwdz_i,nx,ny,nz1,dzz,dwdz)!dwdz(i,j,k)=(dwdz_i(i,j,k+1)+dwdz_i(i,j,k))*0.5

		!===================================================================================================================
		if (m_turb == 1) then ! 1 = LES Smagorinsky-Lilly clássico

			ka  = 0.

			!$omp parallel do collapse(3) private(p1,p2,p3,sij)
			do k=1,nz
			do j=1,ny
			do i=1,nx
	   			p1=dudx(i,j,k)+0.5*(dudy(i,j,k)+dvdx(i,j,k)+dudz(i,j,k)+dwdx(i,j,k)) ! p1+p2+p3 = Sij
	   			p2=dvdy(i,j,k)+0.5*(dvdx(i,j,k)+dudy(i,j,k)+dvdz(i,j,k)+dwdy(i,j,k))
	   			p3=dwdz(i,j,k)+0.5*(dwdy(i,j,k)+dvdz(i,j,k)+dudz(i,j,k)+dwdx(i,j,k))
	   			sij = p1+p2+p3
	   			!mut(i,j,k)=delta*delta*csmag*csmag*dsqrt(2.*(p1+p2+p3)*(p1+p2+p3))*rho(i,j,k) 
	   			mut(i,j,k)=deltai(k)*deltai(k)*csmag*csmag*dsqrt(2.d0*sij*sij)*rho(i,j,k) ! viscosidade dinâmica modelada
			enddo
			enddo
			enddo
			!$omp end parallel do

		!===================================================================================================================
		elseif (m_turb == 5) then ! LES Smagorinsky dinâmico (Germano-Lilly)

			ka = 0.d0
			call visco_smag_dinamico(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)

		!===================================================================================================================
		elseif (m_turb == 4) then ! 4 = k-omega SST (Menter)

			! Tensor taxa de deformação (centro da célula)
			! S2 = 2*Sij*Sij (magnitude ao quadrado)
			!$omp parallel do collapse(3) private(p1,p2,p3)
			do k=1,nz
			do j=1,ny
			do i=1,nx
				p1 = dudx(i,j,k)
				p2 = dvdy(i,j,k)
				p3 = dwdz(i,j,k)
				S2(i,j,k) = 2.d0*(p1*p1 + p2*p2 + p3*p3 + 2.d0*( &
					(0.5d0*(dudy(i,j,k)+dvdx(i,j,k)))**2 + &
					(0.5d0*(dudz(i,j,k)+dwdx(i,j,k)))**2 + &
					(0.5d0*(dvdz(i,j,k)+dwdy(i,j,k)))**2 ))
				Smag(i,j,k) = dsqrt(max(S2(i,j,k), zero))
			enddo
			enddo
			enddo
			!$omp end parallel do

			do tt = 1, ntt
				dt = a_dt(tt)

				CALL contorno_sst()

				! Garantias de positividade/robustez no interior (max() não elimina NaN)
				!$omp parallel do collapse(3)
				do k = 1, nz
				do j = 1, ny
				do i = 1, nx
					if (ka(i,j,k) /= ka(i,j,k)) ka(i,j,k) = 0.d0
					if (omg(i,j,k) /= omg(i,j,k)) omg(i,j,k) = omega_min
					ka(i,j,k)  = max(ka(i,j,k), 0.d0)
					omg(i,j,k) = max(omg(i,j,k), omega_min)
				enddo
				enddo
				enddo
				!$omp end parallel do

				! Gradientes de k e omega para termo de difusão cruzada
				call derivax(ka,nx,ny,nz,dkdx)
				call derivay(ka,nx,ny,nz,dkdy)
				call derivaz(ka,nx,ny,nz,inv_2dz,dkdz)
				call derivax(omg,nx,ny,nz,dodx)
				call derivay(omg,nx,ny,nz,dody)
				call derivaz(omg,nx,ny,nz,inv_2dz,dodz)

				! Blending + sigmas + nut/mut
				!$omp parallel do collapse(3) private(ydist,nu_loc,k_loc,omega_loc,arg1)
				do k = 1, nz
				do j = 1, ny
				do i = 1, nx
					ydist = max(distw(i,j,k), 1.0d-12)                         ! y (distance-to-wall)
					nu_loc = max(ls_mu(i,j,k)/rho(i,j,k), 1.0d-12)             ! nu (cinemática)
					k_loc = max(ka(i,j,k), zero)                               ! k
					omega_loc = max(omg(i,j,k), omega_min)                     ! omega

					! CDkw (forma cinemática)
					diss(i,j,k) = 2.d0*sigma_om2*(dkdx(i,j,k)*dodx(i,j,k) + dkdy(i,j,k)*dody(i,j,k) + dkdz(i,j,k)*dodz(i,j,k)) / omega_loc
					diss(i,j,k) = max(diss(i,j,k), 1.0d-10)

					! F1
					des(i,j,k) = dsqrt(k_loc)/(betaStar*omega_loc*ydist)      ! term1
					F2(i,j,k)  = 500.d0*nu_loc/(ydist*ydist*omega_loc)        ! term2 (usa F2 como temporário)
					F1(i,j,k)  = 4.d0*sigma_om2*k_loc/(diss(i,j,k)*ydist*ydist)! term3 (usa F1 como temporário)
					arg1 = min(max(des(i,j,k), F2(i,j,k)), F1(i,j,k))          ! arg1
					F1(i,j,k) = dtanh((arg1*arg1)*(arg1*arg1))

					! F2
					arg1 = max( &
						2.d0*dsqrt(k_loc)/(betaStar*omega_loc*ydist), &
						500.d0*nu_loc/(ydist*ydist*omega_loc) )
					F2(i,j,k) = dtanh(arg1*arg1)

					sigk(i,j,k)  = F1(i,j,k)*sigma_k1  + (one - F1(i,j,k))*sigma_k2
					sigom(i,j,k) = F1(i,j,k)*sigma_om1 + (one - F1(i,j,k))*sigma_om2

					nut(i,j,k) = a1_sst*k_loc/max(a1_sst*omega_loc, Smag(i,j,k)*F2(i,j,k) + tiny)
					! Limiter (evita viscosidade turbulenta patológica por omega muito baixo)
					nut(i,j,k) = min(nut(i,j,k), 1.0d4*max(ls_mu(i,j,k)/rho(i,j,k), 1.0d-12))
					mut(i,j,k) = nut(i,j,k)*rho(i,j,k)
				enddo
				enddo
				enddo
				!$omp end parallel do

				! k: convecção + difusão
				CALL convdiff_sst_k(Fka, sigk)
				!$omp parallel do collapse(3) private(k_loc,omega_loc,Pk)
				do k = 1, nz
				do j = 1, ny
				do i = 1, nx
					k_loc = max(ka(i,j,k), zero)
					omega_loc = max(omg(i,j,k), omega_min)
					Pk = nut(i,j,k)*S2(i,j,k)
					Pk = min(Pk, 10.d0*betaStar*k_loc*omega_loc)
					Fka(i,j,k) = Fka(i,j,k) + Pk - betaStar*k_loc*omega_loc
				enddo
				enddo
				enddo
				!$omp end parallel do

				CALL tempo_esc(Fka,nx,ny,nz,ka(1:nx,1:ny,1:nz),e0,fe0,fe1)
				!$omp parallel do collapse(3)
				do k = 1, nz
				do j = 1, ny
				do i = 1, nx
					if (ka(i,j,k) /= ka(i,j,k)) then
						ka(i,j,k) = 0.d0
					else
						ka(i,j,k) = max(ka(i,j,k), zero)
					endif
				enddo
				enddo
				enddo
				!$omp end parallel do

				CALL contorno_sst()

				! Recalcula gradientes (k mudou) para cross-diffusion em omega
				call derivax(ka,nx,ny,nz,dkdx)
				call derivay(ka,nx,ny,nz,dkdy)
				call derivaz(ka,nx,ny,nz,inv_2dz,dkdz)
				call derivax(omg,nx,ny,nz,dodx)
				call derivay(omg,nx,ny,nz,dody)
				call derivaz(omg,nx,ny,nz,inv_2dz,dodz)

				! omega: convecção + difusão
				CALL convdiff_sst_omega(Fom, sigom)
				!$omp parallel do collapse(3) private(k_loc,omega_loc,Pk,gamma_loc,beta_loc)
				do k = 1, nz
				do j = 1, ny
				do i = 1, nx
					k_loc = max(ka(i,j,k), zero)
					omega_loc = max(omg(i,j,k), omega_min)
					Pk = nut(i,j,k)*S2(i,j,k)
					Pk = min(Pk, 10.d0*betaStar*k_loc*omega_loc)

					! gamma e beta (omega eq)
					diss(i,j,k) = F1(i,j,k)*gamma1 + (one - F1(i,j,k))*gamma2
					des(i,j,k)  = F1(i,j,k)*beta1  + (one - F1(i,j,k))*beta2
					gamma_loc = diss(i,j,k)
					beta_loc  = des(i,j,k)

					! produção e dissipação
					if (nut(i,j,k) > 1.0d-12) then
						Fom(i,j,k) = Fom(i,j,k) + gamma_loc*Pk/max(nut(i,j,k), 1.0d-12) - beta_loc*omega_loc*omega_loc
					else
						Fom(i,j,k) = Fom(i,j,k) - beta_loc*omega_loc*omega_loc
					endif

					! difusão cruzada (SST)
					Fom(i,j,k) = Fom(i,j,k) + 2.d0*(one - F1(i,j,k))*sigma_om2 * &
						(dkdx(i,j,k)*dodx(i,j,k) + dkdy(i,j,k)*dody(i,j,k) + dkdz(i,j,k)*dodz(i,j,k)) / omega_loc
				enddo
				enddo
				enddo
				!$omp end parallel do

				CALL tempo_esc(Fom,nx,ny,nz,omg(1:nx,1:ny,1:nz),e0o,fe0o,fe1o)
				!$omp parallel do collapse(3)
				do k = 1, nz
				do j = 1, ny
				do i = 1, nx
					if (omg(i,j,k) /= omg(i,j,k)) then
						omg(i,j,k) = omega_min
					else
						omg(i,j,k) = max(omg(i,j,k), omega_min)
					endif
				enddo
				enddo
				enddo
				!$omp end parallel do

				CALL contorno_sst()
			enddo

		!===================================================================================================================
		elseif (m_turb == 2 .or. m_turb == 3) then ! 2 = LES - energia cinemática turbulenta, 3 = DES (Heinz, 2020)

		! comprimento característico	
		lcar = max(prof,ampl)

			do tt = 1, ntt

				dt   = a_dt(tt)
				
				!$omp parallel do collapse(3) private(p1,p2)
				do k = 1, nz 
				do j = 1, ny
				do i = 1, nx
						
					p1 = dsqrt(max(ka(i,j,k), zero))   ! sqrt(k)
					p2 = max(ka(i,j,k), zero) * p1      ! k^1.5 (clip: evita 0/0 abaixo)
					nut(i,j,k)  = cmu*p1*deltai(k)
					diss(i,j,k) = max(p2/lcar, 1.0d-20) ! taxa de dissipação > denominador = comprimento caracteristico
						! sem obstáculo > comprimento caracteristico = prof
						! com obstáculo > comprimento caracteristico = 0.2m altura do degrau (backward facing step/T1 - Delft, 1980)
					des(i,j,k) = p2/(cdes*deltai(k)*diss(i,j,k))	! proteção: diss>0 evita divisão por zero (FPE)
					
				enddo
				enddo
				enddo
				!$omp end parallel do
				
				if (m_turb == 3) then
					! modificação da equação k para DES
					!$omp parallel do collapse(3)
					do k = 1, nz
					do j = 1, ny
					do i = 1, nx
						des(i,j,k) = max(one,des(i,j,k))
					enddo
					enddo
					enddo
					!$omp end parallel do
				else
					!$omp parallel do collapse(3)
					do k = 1, nz
					do j = 1, ny
					do i = 1, nx
						des(i,j,k) = one
					enddo
					enddo
					enddo
					!$omp end parallel do
				endif

				CALL convdiff_les(Fka)       ! termos convectivos e difusivos da energia cinética


				!$omp parallel do collapse(3) private(p1,p2,p3,sij)
				do k = 1, nz                 ! adição dos outros termos da energia cinética
				do j = 1, ny
				do i = 1, nx
					p1=dudx(i,j,k)+0.5*(dudy(i,j,k)+dvdx(i,j,k)+dudz(i,j,k)+dwdx(i,j,k)) ! p1+p2+p3 = Sij
	   				p2=dvdy(i,j,k)+0.5*(dvdx(i,j,k)+dudy(i,j,k)+dvdz(i,j,k)+dwdy(i,j,k))
	   				p3=dwdz(i,j,k)+0.5*(dwdy(i,j,k)+dvdz(i,j,k)+dudz(i,j,k)+dwdx(i,j,k))
	   				sij = p1+p2+p3
										
					Fka(i,j,k) = Fka(i,j,k) + nut(i,j,k)*2.d0*sij*sij - des(i,j,k)*diss(i,j,k)
				enddo
				enddo
				enddo
				!$omp end parallel do

				CALL tempo_esc(Fka,nx,ny,nz,ka(1:nx,1:ny,1:nz),e0,fe0,fe1) ! evolução temporal da energia cinética turbulenta (ka)

				CALl contorno_les


				!$omp parallel do collapse(3)
				do k = 1, nz
				do j = 1, ny
				do i = 1, nx				
					ka(i,j,k) = max(ka(i,j,k), zero)
				enddo
				enddo
				enddo
                !$omp end parallel do


			enddo

!========================================================================================
!		        nut = cmu*ka(1:nx,1:ny,1:nz)**(0.5)*delta ! viscosidade modelada
!			mut = nut*rho     ! viscosidade dinâmica modelada
!========================================================================================


		!$omp parallel do collapse(3) private(p1)
                do k = 1, nz
                do j = 1, ny
                do i = 1, nx
					p1 = dsqrt(max(ka(i,j,k), zero)) ! sqrt(k)
		         	nut(i,j,k) = cmu*p1*deltai(k) ! viscosidade modelada
			 		mut(i,j,k) = nut(i,j,k)*rho(i,j,k)     ! viscosidade dinâmica modelada
                enddo
                enddo
                enddo
                !$omp end parallel do

		endif
	endif

!===================================================================================================================


	!$omp parallel do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
	        mu(i,j,k) = ls_mu(i,j,k) + mut(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

	call interpx_cf(mu,nx,ny,nz,xmu)
	call interpy_cf(mu,nx,ny,nz,ymu)
	call interpz_cf(mu,nx,ny,nz,dz,zmu)

END SUBROUTINE visco


SUBROUTINE visco_smag_dinamico(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz)

	USE disc
	USE cond
	USE les
	USE velpre

	IMPLICIT NONE

	real(8),intent(in),dimension(nx,ny,nz) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz
	integer :: i,j,k,ii,jj,kk,di,dj,dk,a,b
	real(8) :: sumw,weight,uu,vv,ww,smag_local,smag_test,delta_test
	real(8) :: trace_l,trace_m,numerator,denominator,cs2
	real(8) :: sij(3,3),sbar(3,3),qbar(3,3),lten(3,3),mten(3,3),ubar(3),uij(3,3),uvec(3)
	real(8) :: qij(3,3),delta_local

	! Filtro teste caixa 3x3x3. Nas direções não periódicas, a janela é
	! normalizada com as células disponíveis; em z ela é ponderada por dz.
	! O coeficiente dinâmico é limitado para evitar backscatter instável.
	!!$omp parallel do collapse(3) default(none) &
	!!$omp& shared(dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,nx,ny,nz,ccx0,ccy0,u,v,w,dz,deltai,rho,mut) &
	!!$omp& private(i,j,k,ii,jj,kk,di,dj,dk,a,b,sumw,weight,uu,vv,ww,smag_local,smag_test,delta_test,trace_l,trace_m,numerator,denominator,cs2,sij,sbar,qbar,lten,mten,ubar,uij,uvec,qij,delta_local)
	do k=1,nz
	do j=1,ny
	do i=1,nx
		ubar = 0.d0
		uij = 0.d0
		sbar = 0.d0
		qbar = 0.d0
		sumw = 0.d0

		do dk=-1,1
			kk=k+dk
			if (kk < 1 .or. kk > nz) cycle
			do dj=-1,1
				jj=j+dj
				if (ccy0 == 0) then
					jj=modulo(jj-1,ny)+1
				elseif (jj < 1 .or. jj > ny) then
					cycle
				endif
				do di=-1,1
					ii=i+di
					if (ccx0 == 0) then
						ii=modulo(ii-1,nx)+1
					elseif (ii < 1 .or. ii > nx) then
						cycle
					endif

					weight=dz(kk)
					uu=0.5d0*(u(ii,jj,kk)+u(ii+1,jj,kk))
					vv=0.5d0*(v(ii,jj,kk)+v(ii,jj+1,kk))
					ww=0.5d0*(w(ii,jj,kk)+w(ii,jj,kk+1))
					uvec=(/uu,vv,ww/)

					sij=0.d0
					sij(1,1)=dudx(ii,jj,kk)
					sij(2,2)=dvdy(ii,jj,kk)
					sij(3,3)=dwdz(ii,jj,kk)
					sij(1,2)=0.5d0*(dudy(ii,jj,kk)+dvdx(ii,jj,kk))
					sij(2,1)=sij(1,2)
					sij(1,3)=0.5d0*(dudz(ii,jj,kk)+dwdx(ii,jj,kk))
					sij(3,1)=sij(1,3)
					sij(2,3)=0.5d0*(dvdz(ii,jj,kk)+dwdy(ii,jj,kk))
					sij(3,2)=sij(2,3)
					smag_local=dsqrt(2.d0*sum(sij*sij))
					delta_local=deltai(kk)
					qij=delta_local*delta_local*smag_local*sij

					sumw=sumw+weight
					ubar(1)=ubar(1)+weight*uu
					ubar(2)=ubar(2)+weight*vv
					ubar(3)=ubar(3)+weight*ww
					uij = uij + weight * spread(uvec,1,3) * spread(uvec,2,3)
					sbar = sbar + weight * sij
					qbar = qbar + weight * qij
				enddo
			enddo
		enddo

		if (sumw <= 0.d0) then
			mut(i,j,k)=0.d0
			cycle
		endif
		ubar=ubar/sumw
		uij=uij/sumw
		sbar=sbar/sumw
		qbar=qbar/sumw

		lten=uij
		do a=1,3
			do b=1,3
				lten(a,b)=lten(a,b)-ubar(a)*ubar(b)
			enddo
		enddo
		trace_l=(lten(1,1)+lten(2,2)+lten(3,3))/3.d0
		do a=1,3
			lten(a,a)=lten(a,a)-trace_l
		enddo

		smag_test=dsqrt(2.d0*sum(sbar*sbar))
		delta_test=2.d0*deltai(k)
		mten=2.d0*(qbar-delta_test*delta_test*smag_test*sbar)
		trace_m=(mten(1,1)+mten(2,2)+mten(3,3))/3.d0
		do a=1,3
			mten(a,a)=mten(a,a)-trace_m
		enddo

		numerator=sum(lten*mten)
		denominator=sum(mten*mten)
		if (denominator > 1.d-30) then
			cs2=max(0.d0,min(numerator/denominator,0.23d0*0.23d0))
		else
			cs2=0.d0
		endif

		sij=0.d0
		sij(1,1)=dudx(i,j,k)
		sij(2,2)=dvdy(i,j,k)
		sij(3,3)=dwdz(i,j,k)
		sij(1,2)=0.5d0*(dudy(i,j,k)+dvdx(i,j,k))
		sij(2,1)=sij(1,2)
		sij(1,3)=0.5d0*(dudz(i,j,k)+dwdx(i,j,k))
		sij(3,1)=sij(1,3)
		sij(2,3)=0.5d0*(dvdz(i,j,k)+dwdy(i,j,k))
		sij(3,2)=sij(2,3)
		smag_local=dsqrt(2.d0*sum(sij*sij))
		mut(i,j,k)=cs2*deltai(k)*deltai(k)*smag_local*rho(i,j,k)
	enddo
	enddo
	enddo
	!!$omp end parallel do

END SUBROUTINE visco_smag_dinamico
