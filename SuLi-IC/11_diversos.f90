!Subrotina para calcular diversos


!Implementação em 03/07/2017
! Leonardo Romero Monteiro

!Modificações
!Leonardo Romero Monteiro em 13/07/2017

SUBROUTINE tempo()

	USE vartempo
	USE velpre

	IMPLICIT NONE

	!Contadores
	integer :: i, j, k
	real(8),save :: aux1

	if ((t_tempo == 0) .or. ((t_tempo == 3) .or. (t_tempo == 4) .and. (it == 1))) then !Euler e primeiro do AB2 e AB3
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		u(i,j,k) = u(i,j,k) + dt*Fu(i,j,k)

		enddo
		enddo	
		enddo

		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		v(i,j,k) = v(i,j,k) + dt*Fv(i,j,k)

		enddo
		enddo
		enddo

		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		w(i,j,k) = w(i,j,k) + dt*Fw(i,j,k)

		enddo
		enddo
		enddo
		
		if (t_tempo == 3 .or. t_tempo == 4) then !Só para AB2
			fu0 = Fu
			fv0 = Fv
			fw0 = Fw
		endif
	elseif (t_tempo == 1) then !RK 2
		if (tt == 1) then
			aux1 = dt !Já está modificado para ser 0.5 do dt original

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u0(i,j,k)  = u(i,j,k)
			u(i,j,k)   = u(i,j,k) + aux1*Fu(i,j,k)
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v0(i,j,k)  = v(i,j,k)
			v(i,j,k)   = v(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w0(i,j,k)  = w(i,j,k)
			w(i,j,k)   = w(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo

		elseif (tt == 2) then
			aux1 = dt

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u(i,j,k) = u0(i,j,k) +aux1 * Fu(i,j,k)	
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v(i,j,k) = v0(i,j,k) +aux1 * Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w(i,j,k) = w0(i,j,k) +aux1 * Fw(i,j,k)
			enddo
			enddo
			enddo
		endif
	elseif (t_tempo == 2) then !RK 3
		if (tt == 1) then
			aux1 = dt !Já está modificado para ser 0.5 do dt original

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u0(i,j,k)  = u(i,j,k)
			fu0(i,j,k) = Fu(i,j,k)		
			u(i,j,k)   = u(i,j,k) + aux1*Fu(i,j,k)
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v0(i,j,k)  = v(i,j,k)
			fv0(i,j,k) = Fv(i,j,k)	
			v(i,j,k)   = v(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w0(i,j,k)  = w(i,j,k)
			fw0(i,j,k) = Fw(i,j,k)	
			w(i,j,k)   = w(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo
		elseif (tt == 2) then
			aux1 = 2.*dt
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			fu1(i,j,k) = Fu(i,j,k)		
			u(i,j,k)   = u0(i,j,k) - dt*fu0(i,j,k) + aux1*Fu(i,j,k)	
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			fv1(i,j,k) = Fv(i,j,k)	
			v(i,j,k)   = v0(i,j,k) - dt*fv0(i,j,k) + aux1*Fv(i,j,k)
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			fw1(i,j,k) = Fw(i,j,k)	
			w(i,j,k)   = w0(i,j,k) - dt*fw0(i,j,k) + aux1*Fw(i,j,k)
			enddo
			enddo
			enddo

		elseif (tt == 3) then
			aux1 = dt/6.

			do k = 1, nz
			do j = 1, ny
			do i = 1, nx1
			u(i,j,k) = u0(i,j,k) +aux1 * (fu0(i,j,k) + 4.*fu1(i,j,k) + Fu(i,j,k)) 	
			enddo
			enddo	
			enddo

			do k = 1, nz
			do j = 1, ny1
			do i = 1, nx
			v(i,j,k) = v0(i,j,k) +aux1 * (fv0(i,j,k) + 4.*fv1(i,j,k) + Fv(i,j,k)) 	
			enddo
			enddo
			enddo

			do k = 1, nz1
			do j = 1, ny
			do i = 1, nx
			w(i,j,k) = w0(i,j,k) +aux1 * (fw0(i,j,k) + 4.*fw1(i,j,k) + Fw(i,j,k)) 	
			enddo
			enddo
			enddo
		endif
	elseif (( (t_tempo == 3) .and. (it .ne. 1)).or. (t_tempo == 4) .and. (it .eq. 2)) then ! AB2 e segundo AB3
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		u(i,j,k) = u(i,j,k) + dt*0.5*(3*Fu(i,j,k)-fu0(i,j,k))
		fu0(i,j,k) = Fu(i,j,k)	
		enddo
		enddo	
		enddo

		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		v(i,j,k) = v(i,j,k) + dt*0.5*(3*Fv(i,j,k)-fv0(i,j,k))
		fv0(i,j,k) = Fv(i,j,k)	
		enddo
		enddo
		enddo

		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		w(i,j,k) = w(i,j,k) + dt*0.5*(3*Fw(i,j,k)-fw0(i,j,k))
		fw0(i,j,k) = Fw(i,j,k)	
		enddo
		enddo
		enddo

		if (t_tempo == 4) then
			fu1 = fu0
			fv1 = fv0
			fw1 = fw0
		endif
elseif ( (t_tempo == 4) .and. (it .ge. 3)) then ! AB3

		aux1 = dt/12.

		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
		u(i,j,k) = u(i,j,k) + aux1*(23.*Fu(i,j,k)-16.*fu0(i,j,k)+5.*fu1(i,j,k))
		fu1(i,j,k) = fu0(i,j,k)
		fu0(i,j,k) = Fu(i,j,k)
		enddo
		enddo	
		enddo
	
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
		v(i,j,k) = v(i,j,k) + aux1*(23.*Fv(i,j,k)-16.*fv0(i,j,k)+5.*fv1(i,j,k))
		fv1(i,j,k) = fv0(i,j,k)	
		fv0(i,j,k) = Fv(i,j,k)	
		enddo
		enddo
		enddo
	
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
		w(i,j,k) = w(i,j,k) + aux1*(23.*Fw(i,j,k)-16.*fw0(i,j,k)+5.*fw1(i,j,k))
		fw1(i,j,k) = fw0(i,j,k)	
		fw0(i,j,k) = Fw(i,j,k)	
		enddo
		enddo
		enddo
endif

END SUBROUTINE tempo


!!!!!!!!!!!!################################################
SUBROUTINE complementos()

	USE velpre
	USE param
	USE ls_param
	USE vartempo
	USE mms_m
	USE obst
	IMPLICIT NONE


	!Contadores
	integer :: i, j, k
	real(8),dimension(nx,ny,nz) :: aux
	real(8),dimension(nx1,ny,nz) :: dhsdx
	real(8),dimension(nx,ny1,nz) :: dhsdy
	real(8),dimension(nx,ny,nz1) :: dhsdz, epis_z
	

	!! tensão superficial, MMS, rugosidade, gravidade e camada esponja

	!reinicializando a subrotina
	epis_z = 0.
	
	
	!Camada esponja.
	if (esp_type > 0) then
		CALL sponge_layer(epis_z)
	endif

	aux = hsx*kurv/rho	
	call interpx_cf(aux,nx,ny,nz,dhsdx) !(nx1,ny,nz)
	
	aux = hsy*kurv/rho	
	call interpy_cf(aux,nx,ny,nz,dhsdy) !(nx,ny1,nz)
	
	aux = hsz*kurv/rho	
	call interpz_cf(aux,nx,ny,nz,dhsdz) !(nx,ny,nz1)


	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		Fu(i,j,k) = Fu(i,j,k) + gx &
			    -sigma*dhsdx(i,j,k) +tf_u(i,j,k) - gz/(chezy*chezy)*sqrt(ub(i,j,k)*ub(i,j,k))*u(i,j,k)
	enddo
	enddo	
	
	do j = 1, ny1
	do i = 1, nx
		Fv(i,j,k) = Fv(i,j,k) &
			    -sigma*dhsdy(i,j,k) +tf_v(i,j,k) - gz/(chezy*chezy)*sqrt(vb(i,j,k)*vb(i,j,k))*v(i,j,k)
	enddo
	enddo
	enddo


	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		Fw(i,j,k) = Fw(i,j,k) - (w(i,j,k)-0.)* epis_z(i,j,k) - gz &
			    -sigma*dhsdz(i,j,k) +tf_w(i,j,k) - gz/(chezy*chezy)*sqrt(wb(i,j,k)*wb(i,j,k))*w(i,j,k)
	enddo
	enddo
	enddo

END SUBROUTINE complementos


!!!!!!!!!!!################################################
SUBROUTINE restart_ini()

	!Só deve rodar caso seja restart
	use velpre
	use obst
	use tempo
	use ls_param
	use mms_m
	use smag
	implicit none
	open (unit=6662, action= 'read', file= 'dados//arquivorestart', status= 'old', form='formatted')
	read(6662,*)  u,v,w,ku,kv,kw,prd0,prd1,prd,rho,ls_nu,ls,mod_ls,bxx0,bxx1,bxy0,bxz0,bxxf,bxxf1,bxyf,bxzf,byx0,byy0,&
	&byy1,byz0,byxf,byyf,byyf1,byzf,bzx0,bzy0,bzz0,bzz1,bzxf,bzyf,bzzf,bzzf1,ub,vb,wb,d_max,d_min,ls_m,&
	hsx,hsy,hsz,hs,kurv,xnut,ynut,znut,tf_u,tf_v,tf_w,vol_ini,vol_ins
	 close (unit=6662)

	!Inicialização do código
	it = 0
	t=0.
	cont = 10

	!Método da integração temporal
	write(*,*) "Esquema temporal"
	if ((t_tempo == 0)  .or. (t_tempo == 3)) then
		if (t_tempo == 0) write(*,*) "Euler"
		if (t_tempo == 3) write(*,*) "AB2"
			ntt = 1
			a_dt = 1.0*dt
	elseif (t_tempo == 1) then
		write(*,*) "RK2"
		ntt = 2
		a_dt(1) = 0.5*dt
		a_dt(2) = 1.0*dt
		a_dt(3) = 1.0*dt
	elseif (t_tempo == 2) then
		write(*,*) "RK3"
		ntt = 3
		a_dt(1) = 0.5*dt
		a_dt(2) = 1.0*dt
		a_dt(3) = 1.0*dt
	endif
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
	mi_f1 = 0.000018253 !Pa/s  !0.00001516!m²/s ar (ls negativo) 20°C !!!$$$$$$$$$ no incompact3d tá como NI
	rho_f2 =998.0		!kg/m³ água saturada (ls positivo) 20°C.
	mi_f2 = 0.00100798  !Pa/s  !0.00000101!m²/s água saturada (ls positivo) 20°C !!!$$$$$$$$$$ no incompact3d tá como NI
	sigma = 0.0728		!N/m tensão superficial da água 20°C
		
	rho_m = abs(rho_f2-rho_f1)*0.5

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

END SUBROUTINE restart_ini

SUBROUTINE restart_salva()

	use velpre
	use obst
	use ls_param
	use smag
	use mms_m
	implicit none
	!u,v,w
	!ku,kv,kw
	!prd0,prd1,prd,rho,ls_nu,ls,mod_ls
	open (unit=6661, action= 'write', file= 'dados//arquivorestart', status= 'unknown', form='formatted')
	write(6661,*) u,v,w,ku,kv,kw,prd0,prd1,prd,rho,ls_nu,ls,mod_ls,bxx0,bxx1,bxy0,bxz0,bxxf,bxxf1,bxyf,bxzf,byx0,byy0,&
	&byy1,byz0,byxf,byyf,byyf1,byzf,bzx0,bzy0,bzz0,bzz1,bzxf,bzyf,bzzf,bzzf1,ub,vb,wb,d_max,d_min,ls_m,&
	hsx,hsy,hsz,hs,kurv,xnut,ynut,znut,tf_u,tf_v,tf_w,vol_ini,vol_ins
	close (unit=6661)
	print *, 'Salvou restart'

END SUBROUTINE restart_salva



