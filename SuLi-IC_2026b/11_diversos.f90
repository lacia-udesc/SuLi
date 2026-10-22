!Subrotina para calcular diversos


!Implementação em 03/07/2017
! Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 13/07/2017

SUBROUTINE tempo()

	USE vartempo
	USE velpre
	!$ USE omp_lib

	IMPLICIT NONE

	!Contadores
	integer :: i, j, k
	real(8) :: aux1

	if ((t_press .eq. 2) .and. (mms_t .eq. 0)) then


        !$omp parallel
        !$omp do
        do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			Fu(i,j,k) = Fu(i,j,k) + px(i,j,k)
		enddo
		enddo
		enddo
        !$omp end do nowait

		!$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
			Fv(i,j,k) = Fv(i,j,k) + py(i,j,k)
		enddo
		enddo
		enddo
        !$omp end do nowait

        !$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			Fw(i,j,k) = Fw(i,j,k) + pz(i,j,k)
		enddo
		enddo
		enddo
		!$omp end do
		!$omp end parallel

	endif


	if ((t_tempo == 0) .or. ((t_tempo == 3) .or. (t_tempo == 4) .and. (it == 1))) then !Euler e primeiro do AB2 e AB3


	    !$omp parallel
	    !$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		u(i,j,k) = u(i,j,k) + dt*Fu(i,j,k)

		enddo
		enddo	
		enddo
		!$omp end do nowait

        !$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		v(i,j,k) = v(i,j,k) + dt*Fv(i,j,k)

		enddo
		enddo
		enddo
		!$omp end do nowait

        !$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		w(i,j,k) = w(i,j,k) + dt*Fw(i,j,k)

		enddo
		enddo
		enddo
		!$omp end do
		!$omp end parallel

		if (t_tempo == 3 .or. t_tempo == 4) then !Só para AB2
			! cópias explícitas (evita atribuições de array fora do OpenMP)
			!$omp parallel
			!$omp do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
				fu0(i,j,k) = Fu(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do nowait

			!$omp do
			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
				fv0(i,j,k) = Fv(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do nowait

			!$omp do
			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
				fw0(i,j,k) = Fw(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do
			!$omp end parallel
		endif
	elseif (t_tempo == 1) then !RK 2
		if (tt == 1) then
			aux1 = dt !Já está modificado para ser 0.5 do dt original

            !$omp parallel
            !$omp do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u0(i,j,k)  = u(i,j,k)
			u(i,j,k)   = u(i,j,k) + aux1*Fu(i,j,k)
			enddo
			enddo	
			enddo
			!$omp end do nowait

            !$omp do
			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v0(i,j,k)  = v(i,j,k)
			v(i,j,k)   = v(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do nowait

            !$omp do
			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w0(i,j,k)  = w(i,j,k)
			w(i,j,k)   = w(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do
			!$omp end parallel


		elseif (tt == 2) then
			aux1 = dt

            !$omp parallel
            !$omp do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u(i,j,k) = u0(i,j,k) +aux1 * Fu(i,j,k)	
			enddo
			enddo	
			enddo
            !$omp end do nowait

            !$omp do
			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v(i,j,k) = v0(i,j,k) +aux1 * Fv(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do nowait

            !$omp do
			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w(i,j,k) = w0(i,j,k) +aux1 * Fw(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do
			!$omp end parallel

		endif
	elseif (t_tempo == 2) then !RK 3
		if (tt == 1) then
			aux1 = dt !Já está modificado para ser 0.5 do dt original


            !$omp parallel
            !$omp do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u0(i,j,k)  = u(i,j,k)
			fu0(i,j,k) = Fu(i,j,k)		
			u(i,j,k)   = u(i,j,k) + aux1*Fu(i,j,k)
			enddo
			enddo	
			enddo
			!$omp end do nowait

            !$omp do
			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v0(i,j,k)  = v(i,j,k)
			fv0(i,j,k) = Fv(i,j,k)	
			v(i,j,k)   = v(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo
            !$omp end do nowait

            !$omp do
			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w0(i,j,k)  = w(i,j,k)
			fw0(i,j,k) = Fw(i,j,k)	
			w(i,j,k)   = w(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do
			!$omp end parallel

		elseif (tt == 2) then
			aux1 = 2.*dt


		    !$omp parallel
			!$omp do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			fu1(i,j,k) = Fu(i,j,k)		
			u(i,j,k)   = u0(i,j,k) - dt*fu0(i,j,k) + aux1*Fu(i,j,k)	
			enddo
			enddo	
			enddo
            !$omp end do nowait

            !$omp do
			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			fv1(i,j,k) = Fv(i,j,k)	
			v(i,j,k)   = v0(i,j,k) - dt*fv0(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo
            !$omp end do nowait

            !$omp do
			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			fw1(i,j,k) = Fw(i,j,k)	
			w(i,j,k)   = w0(i,j,k) - dt*fw0(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo
			!$omp end do
			!$omp end parallel

		elseif (tt == 3) then
			aux1 = dt/6.

            !$omp parallel
            !$omp do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u(i,j,k) = u0(i,j,k) +aux1 * (fu0(i,j,k) + 4.*fu1(i,j,k) + Fu(i,j,k)) 	
			enddo
			enddo	
			enddo
            !$omp end do nowait

            !$omp do
			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v(i,j,k) = v0(i,j,k) +aux1 * (fv0(i,j,k) + 4.*fv1(i,j,k) + Fv(i,j,k)) 	
			enddo
			enddo
			enddo
            !$omp end do nowait

            !$omp do
			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w(i,j,k) = w0(i,j,k) +aux1 * (fw0(i,j,k) + 4.*fw1(i,j,k) + Fw(i,j,k)) 	
			enddo
			enddo
			enddo
			!$omp end do
			!$omp end parallel


		endif
	elseif (( (t_tempo == 3) .and. (it .ne. 1)).or. (t_tempo == 4) .and. (it .eq. 2)) then ! AB2 e segundo AB3

	    !$omp parallel
	    !$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		u(i,j,k) = u(i,j,k) + dt*0.5*(3*Fu(i,j,k)-fu0(i,j,k))
		fu0(i,j,k) = Fu(i,j,k)	
		enddo
		enddo	
		enddo
        !$omp end do nowait

        !$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		v(i,j,k) = v(i,j,k) + dt*0.5*(3*Fv(i,j,k)-fv0(i,j,k))
		fv0(i,j,k) = Fv(i,j,k)	
		enddo
		enddo
		enddo
        !$omp end do nowait

        !$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		w(i,j,k) = w(i,j,k) + dt*0.5*(3*Fw(i,j,k)-fw0(i,j,k))
		fw0(i,j,k) = Fw(i,j,k)	
		enddo
		enddo
		enddo		
        !$omp end do
        !$omp end parallel

		if (t_tempo == 4) then
			fu1 = fu0
			fv1 = fv0
			fw1 = fw0
		endif
elseif ( (t_tempo == 4) .and. (it .ge. 3)) then ! AB3

		aux1 = dt/12.

        !$omp parallel
        !$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		u(i,j,k) = u(i,j,k) + aux1*(23.*Fu(i,j,k)-16.*fu0(i,j,k)+5.*fu1(i,j,k))
		fu1(i,j,k) = fu0(i,j,k)
		fu0(i,j,k) = Fu(i,j,k)
		enddo
		enddo	
		enddo
        !$omp end do nowait

        !$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		v(i,j,k) = v(i,j,k) + aux1*(23.*Fv(i,j,k)-16.*fv0(i,j,k)+5.*fv1(i,j,k))
		fv1(i,j,k) = fv0(i,j,k)	
		fv0(i,j,k) = Fv(i,j,k)	
		enddo
		enddo
		enddo
        !$omp end do nowait

        !$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		w(i,j,k) = w(i,j,k) + aux1*(23.*Fw(i,j,k)-16.*fw0(i,j,k)+5.*fw1(i,j,k))
		fw1(i,j,k) = fw0(i,j,k)	
		fw0(i,j,k) = Fw(i,j,k)	
		enddo
		enddo
		enddo
		!$omp end do
        !$omp end parallel


endif

END SUBROUTINE tempo

!##################################################################################################################

SUBROUTINE tempo_esc(Fe,dimx,dimy,dimz,e,e0,fe0,fe1) ! Evolução temporal para variáveis escalares centralizadas

	USE vartempo
	USE velpre
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k, dimx,dimy,dimz ! contadores
	real(8) :: aux1 ! auxiliar

	real(8),intent(OUT),dimension(dimx,dimy,dimz) :: e ! e = escalar
	real(8),dimension(dimx,dimy,dimz) :: e0,fe0,fe1,Fe

	! RK 3
		if (tt == 1) then
			aux1 = dt ! Já está modificado para ser 0.5 do dt original

            !$omp parallel do
			do k = 1, dimz
			do j = 1, dimy
			do i = 1, dimx
				e0(i,j,k)  = e(i,j,k)
				fe0(i,j,k) = Fe(i,j,k)		
				e(i,j,k)   = e(i,j,k) + aux1*Fe(i,j,k)
			enddo
			enddo
			enddo
            !$omp end parallel do

					
		elseif (tt == 2) then
			aux1 = 2.*dt


            !$omp parallel do
			do k = 1, dimz
			do j = 1, dimy
			do i = 1, dimx
				fe1(i,j,k) = Fe(i,j,k)		
				e(i,j,k)   = e0(i,j,k) - dt*fe0(i,j,k) + aux1*Fe(i,j,k)	
			enddo
			enddo	
			enddo
            !$omp end parallel do


		elseif (tt == 3) then
			aux1 = dt/6.


            !$omp parallel do
			do k = 1, dimz
			do j = 1, dimy
			do i = 1, dimx
				e(i,j,k)   = e0(i,j,k) + aux1*(fe0(i,j,k) + 4.*fe1(i,j,k) + Fe(i,j,k)) 	
			enddo
			enddo
			enddo
            !$omp end parallel do

		endif

END SUBROUTINE tempo_esc

!##################################################################################################################

SUBROUTINE build_distw()
	! Campo distance-to-wall (m) para uso no modelo k-omega SST.
	! Considera apenas paredes com condição de contorno tipo parede (no-slip/semi-slip)
	! e superfícies de obstáculos (obs_lss) quando ibm_t==1 ou ibm_t==2.

	USE disc
	USE cond
	USE obst
	USE les, only: distw

	IMPLICIT NONE

	integer :: i, j, k
	real(8) :: dxw, dyw, dzw, dmin
	real(8), parameter :: huge = 1.0d30
	real(8), parameter :: eps  = 1.0d-12

	!$omp parallel do private(dmin)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		dmin = huge

		! Paredes em X (somente no-slipou semi-slip)
		if (ccx0 == 2 .or. ccx0 == 4) dmin = min(dmin, max(x(i), eps))
		if (ccxf == 2 .or. ccxf == 4) dmin = min(dmin, max(lx - x(i), eps))

		! Paredes em Y (no-slip ou semi-slip)
		if (ccy0 == 2 .or. ccy0 == 4) dmin = min(dmin, max(y(j), eps))
		if (ccyf == 2 .or. ccyf == 4) dmin = min(dmin, max(ly - y(j), eps))

		! Fundo em Z (no-slip ou semi-slip)
		if (ccz0 == 2 .or. ccz0 == 4) dmin = min(dmin, max(z(k), eps))

		! Obstáculos (campo função-distância do IBM, disponível para ibm_t = 1 e 2)
		if ((ibm_t == 1 .or. ibm_t == 2) .and. obst_t /= 0) then
			dmin = min(dmin, max(dabs(obs_lss(i,j,k)), eps))
		endif

		! Se não houver nenhuma parede considerada, mantém um valor grande (evita divisão por zero)
		if (dmin >= huge*0.5d0) dmin = max(z(k), eps)
		
		distw(i,j,k) = dmin
	enddo
	enddo
	enddo
	!$omp end parallel do

END SUBROUTINE build_distw

!##################################################################################################################

SUBROUTINE complementos()

	USE velpre
	USE param
	USE ls_param
	USE vartempo
	USE mms_m
	USE obst
	USE les
	!$ USE omp_lib
	
	IMPLICIT NONE


	!Contadores
	integer :: i, j, k
	real(8),dimension(nx,ny,nz) :: aux
	real(8),dimension(nx1,ny,nz) :: dhsdx, xrho, vx, hsxx
	real(8),dimension(nx,ny1,nz) :: dhsdy, yrho, uy
	real(8),dimension(nx,ny,nz1) :: dhsdz, zrho, epis_z, hszz
	
	real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: karho, rhof

	!! tensão superficial, MMS, rugosidade, gravidade e camada esponja

	!reinicializando a subrotina
	epis_z = 0.
	
	
	!Camada esponja.
	if (esp_type > 0) then
		CALL sponge_layer(epis_z)
	endif

	! Tensão superficial (verificação t_tens adicionada para correção lógica)
	if (t_tens == 1) then
		aux = hsx*kurv/rho	
		call interpx_cf(aux,nx,ny,nz,dhsdx) !(nx1,ny,nz)
		
		aux = hsy*kurv/rho	
		call interpy_cf(aux,nx,ny,nz,dhsdy) !(nx,ny1,nz)
		
		aux = hsz*kurv/rho	
		call interpz_cf(aux,nx,ny,nz,dz,dhsdz) !(nx,ny,nz1)
	else
		! Zerar termo de tensão superficial quando desativado
		dhsdx = 0.
		dhsdy = 0.
		dhsdz = 0.
	endif


	call interpx_cf(hs,nx,ny,nz,hsxx) !(nx1,ny,nz)
	call interpz_cf(hs,nx,ny,nz,dz,hszz) !(nx,ny,nz1)
	
	!Interpolação velocidades	
	!CALL interpxy_ff(ub,nx,ny,nz,nx1,ny1,uy)
	!CALL interpyx_ff(vb,nx,ny,nz,nx1,ny1,vx)
	
	!Energia cinética
	call fantasmas_cf(rho,nx,ny,nz,rhof)
	karho = ka*rhof
	
	call interpx_cf(rho,nx,ny,nz,xrho) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,yrho) !(nx,ny1,nz)	
	call interpz_cf(rho,nx,ny,nz,dz,zrho) !(nx,ny,nz1)


        !$omp parallel
        !$omp do
        do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			Fu(i,j,k) = Fu(i,j,k) & 
				    + gx*hsxx(i,j,k)   &  ! força da gravidade na água!
			    	    - sigma*dhsdx(i,j,k) &  ! tensão superficial (zero se t_tens == 0)
			    	    + tf_u(i,j,k) &  ! MMS
			    	    - 2./3.*(karho(i,j,k)-karho(i-1,j,k))/(dx*xrho(i,j,k)) !& ! energia cinética
			    	    !- ((1./0.41)*log(dy/0.016)**(-2.) * sqrt( ub(i,j,k)*ub(i,j,k) + vx(i,j,k)*vx(i,j,k) ) * ub(i,j,k) ) /dz !*xrho(i,j,k) ! rugosidade do fundo (Ferziger et al., 2020)
			    	    !- ( gz / (chezy*chezy) * sqrt( ub(i,j,k)*ub(i,j,k) + vx(i,j,k)*vx(i,j,k) ) * u(i,j,k) ) / dz  ! rugosidade do fundo (Casulli, Zanolli, 2022)
		enddo
		enddo
		enddo
        !$omp end do nowait

        !$omp do
        do k = 1, nz
        do j = 1, ny1
		do i = 1, nx
			Fv(i,j,k) = Fv(i,j,k) &
				    - sigma*dhsdy(i,j,k) &  ! tensão superficial (zero se t_tens == 0)
				    + tf_v(i,j,k) &  ! MMS
				    - 2./3.*(karho(i,j,k)-karho(i,j-1,k))/(dy*yrho(i,j,k)) !& ! energia cinética
				    !- ((1./0.41)*log(dy/0.016)**(-2.) * sqrt( uy(i,j,k)*uy(i,j,k) + vb(i,j,k)*vb(i,j,k) ) * vb(i,j,k) ) /dz !*yrho(i,j,k) ! rugosidade do fundo (Ferziger et al., 2020)
				    !- ( gz / (chezy*chezy) * sqrt( uy(i,j,k)*uy(i,j,k) + vb(i,j,k)*vb(i,j,k) ) * v(i,j,k) ) / dz  ! rugosidade do fundo (Casulli, Zanolli, 2022)
		enddo
		enddo
		enddo
		!$omp end do nowait

        !$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			Fw(i,j,k) = Fw(i,j,k) &
				    - gz*hszz(i,j,k) &  ! força da gravidade na água!
				    - (w(i,j,k)-0.)*epis_z(i,j,k) & 
				    - sigma*dhsdz(i,j,k) &  ! tensão superficial (zero se t_tens == 0) 
			    	    + tf_w(i,j,k) &  ! MMS
			    	    - 2./3.*(karho(i,j,k)-karho(i,j,k-1))/(dz(k)*zrho(i,j,k)) ! energia cinética
		enddo 
		enddo
		enddo
		!$omp end do
		!$omp end parallel


	
END SUBROUTINE complementos


!##################################################################################################################

SUBROUTINE ruido_turb() 


	USE velpre
	USE les
	
	IMPLICIT NONE

	integer :: i, j, k
	real(8) :: aux1
	real(8),dimension(3):: rr
	

	if (m_turb .ge. 1) then
		if (uinicial>0) then
			aux1 = sqrt(iturb*iturb*(uinicial*uinicial+10.**(-6))*1.5) ! ka = turbulence kinetic energy ! ruído (10-6) adicionado para casos sem velocidade inicial.
			!$omp parallel do private(rr)
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				!ka(i,j,k) = aux1
				call random_number(rr)
				rr = (rr-0.5)
				u(i,j,k) = u(i,j,k)+ aux1*rr(1)
				v(i,j,k) = v(i,j,k)+ aux1*rr(2)
				w(i,j,k) = w(i,j,k)+ aux1*rr(3)		
			enddo
			enddo
			enddo
			!$omp end parallel do
		else
			!$omp parallel do private(rr)
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				!ka(i,j,k) = aux1
				call random_number(rr)
				rr = (rr-0.5)
				u(i,j,k) = u(i,j,k)*(1.d0 + iturb*rr(1))
				v(i,j,k) = v(i,j,k)*(1.d0 + iturb*rr(2))
				w(i,j,k) = w(i,j,k)*(1.d0 + iturb*rr(3))		
			enddo
			enddo
			enddo
			!$omp end parallel do
		endif
	endif
	
END SUBROUTINE ruido_turb



!##################################################################################################################


SUBROUTINE divergente(diver) ! apenas para IPA !

	USE velpre
	USE obst
	!$ USE omp_lib

	implicit none

	real(8),intent(out),dimension(nx,ny,nz)  :: diver
	real(8),dimension(nx1,ny,nz) :: u1
	real(8),dimension(nx,ny1,nz) :: v1
	real(8),dimension(nx,ny,nz1) :: w1
	
	integer :: i, j, k

 	if (obst_t .ne. 0 .and. (ibm_t == 1 .or. ibm_t == 2)) then


        !$omp parallel
        !$omp do
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
        !$omp end do nowait

        !$omp do
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
        !$omp end do nowait

        !$omp do
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
		!$omp end do
		!$omp end parallel

        !$omp parallel do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			diver(i,j,k) = ( u1(i,j,k) - u1(i+1,j,k) )*inv_dx + ( v1(i,j,k)-v1(i,j+1,k) )*inv_dy  + (w1(i,j,k) - w1(i,j,k+1)) /dzz(k)
		enddo
		enddo
		enddo
        !$omp end parallel do

	else

        !$omp parallel do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			diver(i,j,k) = ( u(i,j,k) - u(i+1,j,k) )*inv_dx + ( v(i,j,k)-v(i,j+1,k) )*inv_dy  + (w(i,j,k) - w(i,j,k+1)) /dzz(k)
		enddo
		enddo
		enddo
        !$omp end parallel do

	endif	

END SUBROUTINE divergente


!##################################################################################################################


SUBROUTINE exp_press() ! apenas para IPA !

USE velpre
!$ USE omp_lib

implicit none
	real(8),dimension(nx1,ny,nz) :: prdx, rhox, drhodx
	real(8),dimension(nx,ny1,nz) :: prdy, rhoy, drhody
	real(8),dimension(nx,ny,nz1) :: prdz, rhoz, drhodz
	real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: prdaux
	integer :: i, j, k
	
	prd0 = prd1 + prd0  !inscremento da pressão
	prdaux = prd1 + prd0  !variável para alterar px,py e pz, com um incremento penalizado
	
	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)

	!$omp parallel
	!$omp do
        do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		px(i,j,k) = -(prdaux(i,j,k)-prdaux(i-1,j,k))*inv_dx/rhox(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do nowait

	!$omp do
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		py(i,j,k) = -(prdaux(i,j,k)-prdaux(i,j-1,k))*inv_dy/rhoy(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do nowait
        
    !$omp do
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		pz(i,j,k) = -(prdaux(i,j,k)-prdaux(i,j,k-1))/dzz(k)/rhoz(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel

END SUBROUTINE exp_press


!##################################################################################################################

SUBROUTINE restart_ini()

	!Só deve rodar caso seja restart
	use disc
	use restart, only : it_rest
	use velpre
	use obst
	use tempo_hora_data
	use ls_param
	use mms_m
	use les
	use vartempo

	implicit none
	integer :: ios

	! Garante que as pastas de saída existam antes de qualquer leitura/escrita
	call system('mkdir -p restart')
	call system('mkdir -p arquivos')
	call system('mkdir -p dados')

	open (unit=6662, action='read', file='restart//arquivorestart', status='old', form='formatted', iostat=ios)
	if (ios /= 0) then
		write(*,*) 'ERRO: arquivo de restart "restart//arquivorestart" nao encontrado.'
		write(*,*) 'Gere-o primeiro com irest = 0 e interv_rest > 0 no arquivo 00_suli3d.prm.'
		stop 1
	endif

	! Estado temporal e contadores
	read(6662,*,iostat=ios) it, t, dt, cont, contp

	! Parâmetros do level set e coeficientes de turbulência calibrados
	read(6662,*,iostat=ios) ls_m, d_max, d_min, vol_ini, vol_ins, dx1, dt1
	read(6662,*,iostat=ios) csmag, cmu, cdes

	! Campos de velocidade e pressão
	read(6662,*,iostat=ios) u
	read(6662,*,iostat=ios) v
	read(6662,*,iostat=ios) w
	read(6662,*,iostat=ios) prd0
	read(6662,*,iostat=ios) prd1
	read(6662,*,iostat=ios) prd

	! Campos do level set
	read(6662,*,iostat=ios) rho
	read(6662,*,iostat=ios) ls_mu
	read(6662,*,iostat=ios) ls
	read(6662,*,iostat=ios) hs
	read(6662,*,iostat=ios) hsx
	read(6662,*,iostat=ios) hsy
	read(6662,*,iostat=ios) hsz
	read(6662,*,iostat=ios) kurv

	! Condições de contorno
	read(6662,*,iostat=ios) bxx0, bxx1, bxy0, bxz0, bxxf, bxxf1, bxyf, bxzf
	read(6662,*,iostat=ios) byx0, byy0, byy1, byz0, byxf, byyf, byyf1, byzf
	read(6662,*,iostat=ios) bzx0, bzy0, bzz0, bzz1, bzxf, bzyf, bzzf, bzzf1

	! Obstáculos (ku/kv/kw são recalculados por obstaculo(), mas ficam salvos)
	read(6662,*,iostat=ios) ub
	read(6662,*,iostat=ios) vb
	read(6662,*,iostat=ios) wb
	read(6662,*,iostat=ios) ku
	read(6662,*,iostat=ios) kv
	read(6662,*,iostat=ios) kw

	! Turbulência
	read(6662,*,iostat=ios) mu
	read(6662,*,iostat=ios) mut
	read(6662,*,iostat=ios) nut
	read(6662,*,iostat=ios) ka
	read(6662,*,iostat=ios) omg
	read(6662,*,iostat=ios) xmu
	read(6662,*,iostat=ios) ymu
	read(6662,*,iostat=ios) zmu

	! Histórico dos esquemas multi-estágio (AB2/AB3)
	read(6662,*,iostat=ios) fu0
	read(6662,*,iostat=ios) fv0
	read(6662,*,iostat=ios) fw0
	read(6662,*,iostat=ios) fu1
	read(6662,*,iostat=ios) fv1
	read(6662,*,iostat=ios) fw1

	! Termos fonte do MMS
	read(6662,*,iostat=ios) tf_p
	read(6662,*,iostat=ios) tf_u
	read(6662,*,iostat=ios) tf_v
	read(6662,*,iostat=ios) tf_w
	read(6662,*,iostat=ios) tf_px
	read(6662,*,iostat=ios) tf_py
	read(6662,*,iostat=ios) tf_pz

	close (unit=6662)

	if (ios /= 0) then
		write(*,*) 'ERRO: falha na leitura do arquivo de restart (formato incompatível ou arquivo corrompido).'
		stop 1
	endif

	! Guarda a iteração salva para o programa principal continuar o laço temporal
	it_rest = it

	write(*,'(A,I10,A,F12.4)') " Reinicialização a partir de it = ", it, " e t = ", t

	if (it >= ts) then
		write(*,*) 'AVISO: a simulação salva já atingiu o tempo final (ts). O laço temporal não executará iterações.'
	endif

	! Recalcula os campos do IBM (funções distância obs_ls* e índices
	! ku/kv/kw), que não são salvos no arquivo de restart.
	CALL obstaculo()
	CALL build_distw()
	if (m_turb == 2 .or. m_turb == 3) CALL contorno_les()
	if (m_turb == 4) CALL contorno_sst()

	!Método da integração temporal (recalcula ntt e a_dt com o dt restaurado)
	call coef_tempo()

	if (mms_t == 0) then
	   tf_p = 0.
	   tf_u = 0.
	   tf_v = 0.
	   tf_w = 0.
	endif

	if (wave_t > 0) then
		write(*,*) "Com entrada de onda"
		CALL waves_coef()
	else
		write(*,*) "Sem entrada de onda"
	endif

	alpha1 = 1.5 		!Número de células que farão parte da espessura, pois é uma função suave
	rho_f1 = 1.204		!kg/m³ ar (ls negativo) 20°C
	mu_f1  = 0.000018253 	!Pa/s  !0.00001516!m²/s ar (ls negativo) 20°C !!!$$$$$$$$$ no incompact3d tá como NI
	rho_f2 = 998.0		!kg/m³ água saturada (ls positivo) 20°C.
	mu_f2  = 0.00100798  	!Pa/s  !0.00000101!m²/s água saturada (ls positivo) 20°C !!!$$$$$$$$$$ no incompact3d tá como NI
	sigma  = 0.0728		!N/m tensão superficial da água 20°C
		

	!Coeficientes de integração RK3 TVD
	adtl(1)=1.
	bdtl(1)=0.
	gdtl(1)=1.

	adtl(2)=3./4.
	bdtl(2)=1./4.
	gdtl(2)=1./4.

	adtl(3)=1./3.
	bdtl(3)=2./3.
	gdtl(3)=2./3.

	! Regenera o visu.xdmf considerando o tempo total da simulação
	CALL visu()

END SUBROUTINE restart_ini

!##################################################################################################################

SUBROUTINE restart_salva()

	use disc
	use velpre
	use obst
	use ls_param
	use les
	use mms_m
	use vartempo
	use tempo_hora_data
	implicit none

	! Garante que a pasta de restart exista
	call system('mkdir -p restart')

	open (unit=6661, action='write', file='restart//arquivorestart', status='unknown', form='formatted')

	! Estado temporal e contadores
	write(6661,*) it, t, dt, cont, contp

	! Parâmetros do level set e coeficientes de turbulência calibrados
	write(6661,*) ls_m, d_max, d_min, vol_ini, vol_ins, dx1, dt1
	write(6661,*) csmag, cmu, cdes

	! Campos de velocidade e pressão
	write(6661,*) u
	write(6661,*) v
	write(6661,*) w
	write(6661,*) prd0
	write(6661,*) prd1
	write(6661,*) prd

	! Campos do level set
	write(6661,*) rho
	write(6661,*) ls_mu
	write(6661,*) ls
	write(6661,*) hs
	write(6661,*) hsx
	write(6661,*) hsy
	write(6661,*) hsz
	write(6661,*) kurv

	! Condições de contorno
	write(6661,*) bxx0, bxx1, bxy0, bxz0, bxxf, bxxf1, bxyf, bxzf
	write(6661,*) byx0, byy0, byy1, byz0, byxf, byyf, byyf1, byzf
	write(6661,*) bzx0, bzy0, bzz0, bzz1, bzxf, bzyf, bzzf, bzzf1

	! Obstáculos (ku/kv/kw são recalculados por obstaculo(), mas ficam salvos)
	write(6661,*) ub
	write(6661,*) vb
	write(6661,*) wb
	write(6661,*) ku
	write(6661,*) kv
	write(6661,*) kw

	! Turbulência
	write(6661,*) mu
	write(6661,*) mut
	write(6661,*) nut
	write(6661,*) ka
	write(6661,*) omg
	write(6661,*) xmu
	write(6661,*) ymu
	write(6661,*) zmu

	! Histórico dos esquemas multi-estágio (AB2/AB3)
	write(6661,*) fu0
	write(6661,*) fv0
	write(6661,*) fw0
	write(6661,*) fu1
	write(6661,*) fv1
	write(6661,*) fw1

	! Termos fonte do MMS
	write(6661,*) tf_p
	write(6661,*) tf_u
	write(6661,*) tf_v
	write(6661,*) tf_w
	write(6661,*) tf_px
	write(6661,*) tf_py
	write(6661,*) tf_pz

	close (unit=6661)
	write(*,'(A,I10,A,F12.4,A)') " Salvou restart (it = ", it, ", t = ", t, ") em restart//arquivorestart"

END SUBROUTINE restart_salva


!=========================================================================
! Verifica se alguma velocidade (u, v, w) contém NaN.
! Se encontrar, imprime a localização e para a simulação.
!=========================================================================
SUBROUTINE verifica_nan()

	USE disc, only : nx1, ny, nz, nx, ny1, nz1, it, t, omp_t
	USE velpre, only : u
	!$ USE omp_lib

	implicit none

	integer :: i, j, k
	logical :: encontrou_nan

	encontrou_nan = .false.

	if (omp_t == 0) then
		! Varredura sequencial: aborta no primeiro NaN encontrado
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx1
					if (isnan(u(i,j,k))) call relata_nan(i,j,k)
				enddo
			enddo
		enddo
	else
		! Varredura paralela com redução lógica; em caso de NaN, repete
		! serialmente apenas para localizar/relatar o primeiro deles.
		!$omp parallel do reduction(.or.:encontrou_nan)
		do k = 1, nz
			do j = 1, ny
				do i = 1, nx1
					encontrou_nan = encontrou_nan .or. isnan(u(i,j,k))
				enddo
			enddo
		enddo
		!$omp end parallel do

		if (encontrou_nan) then
			do k = 1, nz
				do j = 1, ny
					do i = 1, nx1
						if (isnan(u(i,j,k))) call relata_nan(i,j,k)
					enddo
				enddo
			enddo
		endif
	endif

CONTAINS

	subroutine relata_nan(i,j,k)

		integer, intent(in) :: i, j, k

		write(*,'(A,3I6)') ' *** NaN detectado em u(i,j,k) = ', i, j, k
		write(*,'(A)')       ' ================================================'
		write(*,'(A,I8)')    ' ERRO FATAL: NaN detectado no passo de tempo it =', it
		write(*,'(A,ES14.6)')' Tempo de simulacao t = ', t
		write(*,'(A)')       ' Simulacao interrompida.'
		write(*,'(A)')       ' ================================================'
		STOP 'NaN detectado nas velocidades - simulacao abortada'

	end subroutine relata_nan

END SUBROUTINE verifica_nan



