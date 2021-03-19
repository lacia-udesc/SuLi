	!Subrotina para calcular o Fu, utilizando o ponto de vista lagrangiano
	!Referencia: Casulli (1990, 1992)

	!!! Implementação 15/04/2014
	! Leonardo Romero Monteiro

	!!! Modificações
	! Leonardo Romero Monteiro - 17/08/2020

	!!! Observações
	! no lagrangiano as velocidades no tempo n = t serão calculadas para o tempo n = t-1

	! são calculadas as Fu, Fv e Fw
	SUBROUTINE convdiff()

	USE dzs
	USE parametros
	USE smag
	USE mms_m
	USE velpre

	IMPLICIT NONE

	!===================================================================================================================
	! auxiliares de velocidades: velocidades lagrangianas
	real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: uint
	real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: vint
	real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: wint

	!
	real(8), dimension(nx1,ny,nz) :: Fu, ch_x
	real(8), dimension(nx,ny1,nz) :: Fv, ch_y
	real(8), dimension(nx,ny,nz1) :: Fw, epis_z

	real(8), dimension(0:nx1,ny,nz) :: ama
	real(8), dimension(nx1,ny1,nz) :: amb
	real(8), dimension(nx1,ny,nz1) :: amd
	real(8), dimension(nx1,ny1,nz) :: bma
	real(8), dimension(nx,0:ny1,nz) :: bmb
	real(8), dimension(nx,ny1,nz1) :: bmd

	!contadores
	integer :: i, j, k

	!===================================================================================================================
	
	!reinicializando a rotina
	eta0 = eta1
	prd0 = prd1

	!reinicializando a subrotina
	epis_z = 0.
	ch_x = 0.
	ch_y = 0.
	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================


	!%%%!-- Advectivo --!%%%!
	if (adv_type == 0) then
		call lagr(uint,vint,wint)
	elseif (adv_type == 1) then
		call classico(uint,vint,wint)
	elseif (adv_type == 2) then
		call rotacional(uint,vint,wint)
	elseif (adv_type == 3) then
		call antissim(uint,vint,wint)
	else
		write(*,*) "sem termo advectivo"
		STOP
	endif


	!Camada esponja
	if (esp_type .ne. 0) CALL sponge_layer(epis_z)
	
	!Chezy
	if (chezy_t .ne. 0) CALL chezy(ch_x,ch_y)	
				    

	!Cálculo dos Fs
 	!call interpx_fc(xnut,nx1,ny,nz,ama(1:nx,1:ny,1:nz)) !(nx,ny,nz)
	call interpy_cf(xnut*dt/dy,nx1,ny,nz,bma) !(nx1,ny1,nz)
	call interpx_cf(ynut*dt/dx,nx,ny1,nz,amb) !(nx1,ny1,nz)
	
	!call interpy_fc(ynut,nx,ny1,nz,bmb(1:nx,1:ny,1:nz)) !(nx,ny,nz)
	call interpx_cf(znut*dt/dx,nx,ny,nz1,amd) !(nx1,ny,nz1)
	call interpy_cf(znut*dt/dy,nx,ny,nz1,bmd) !(nx,ny1,nz1)

      ama(1:nx,1:ny,1:nz)=nut(1:nx,1:ny,1:nz)*dt/dx
      bmb(1:nx,1:ny,1:nz)=nut(1:nx,1:ny,1:nz)*dt/dy

      ama(0,:,:)=ama(1,:,:)
      bmb(:,0,:)=bmb(:,1,:)

      ama(nx1,:,:)=ama(nx,:,:)
      bmb(:,ny1,:)=bmb(:,ny,:)

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		Fu(i,j,k) = uint(i,j,k) + (ama(i,j,k)*(u(i+1,j,k)-u(i,j,k))-&
                            ama(i-1,j,k)*(u(i,j,k)-u(i-1,j,k)))/dx + &
                            (bma(i,j+1,k)*(u(i,j+1,k)-u(i,j,k))-&
                            bma(i,j,k)*(u(i,j,k)-u(i,j-1,k)))/dy + gx * dt  &
			     + (-ch_x(i,j,k) +tf_u(i,j,k))*dt
	enddo
	enddo	

	do j = 1, ny1
	do i = 1, nx
		Fv(i,j,k) = vint(i,j,k) + (amb(i+1,j,k)*(v(i+1,j,k)-v(i,j,k))-&
                            amb(i,j,k)*(v(i,j,k)-v(i-1,j,k)))/dx + &
                            (bmb(i,j,k)*(v(i,j+1,k)-v(i,j,k))-&
                            bmb(i,j-1,k)*(v(i,j,k)-v(i,j-1,k)))/dy  &
			     + (-ch_y(i,j,k) + tf_v(i,j,k))*dt

	enddo
	enddo
	enddo

	! a gravidade está implícita em z por causa da pressão hidrostática
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		Fw(i,j,k) = wint(i,j,k) + (amd(i+1,j,k)*(w(i+1,j,k)-w(i,j,k))-&
                            amd(i,j,k)*(w(i,j,k)-w(i-1,j,k)))/dx + &
                            (bmd(i,j+1,k)*(w(i,j+1,k)-w(i,j,k))-&
                            bmd(i,j,k)*(w(i,j,k)-w(i,j-1,k)))/dy &
                            +(-(w(i,j,k)-0.)*epis_z(i,j,k) +tf_w(i,j,k))*dt
	enddo
	enddo
	enddo

	CALL matrizg(Fu, Fv, Fw)

!==================================================================================================================
END SUBROUTINE convdiff


!####################################################################################################################
!####################################################################################################################


SUBROUTINE matrizg(Fu, Fv, Fw)

	USE dzs
	USE velpre
	USE parametros

	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO TAMBÉM NO PROGRAMA

	real(8), dimension(nx1,ny,nz), intent(in) :: Fu
	real(8), dimension(nx,ny1,nz), intent(in) :: Fv
	real(8), dimension(nx,ny,nz1), intent(in) :: Fw


	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8), dimension(nx1,ny,nz) :: matrizgi, parte1i
	real(8), dimension(nx,ny1,nz) :: matrizgj, parte3j
	real(8), dimension(nx,ny,nz1) :: matrizgk

	real(8), dimension(nx1,ny,nz,nz) :: matrizainvi
	real(8), dimension(nx,ny1,nz,nz) :: matrizainvj	
	real(8), dimension(nx,ny,nz1,nz1) :: matrizainvk

	real(8), dimension(nx1,ny) :: parte1, parte2
	real(8), dimension(nx,ny1) :: parte3, parte4

	!contadores
	integer :: i, j, k, kk

	!auxiliares
	real(8) :: aux1, aux2

	!===================================================================================================================
	
	!%%%!-- Matrizes G --!%%%!

	aux1 = (1-tetah)*dt/dx
	aux2 = (1-tetah)*dt/dy

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		matrizgi(i,j,k)= dzx(i,j,k)*(Fu(i,j,k) -aux1*gz*(eta0(i,j)-eta0(i-1,j)) &
			        - aux1*(prd0(i,j,k)-prd0(i-1,j,k)))

	enddo
    	enddo

	do j = 1, ny1
	do i = 1, nx
		matrizgj(i,j,k)= dzy(i,j,k)*(Fv(i,j,k) -aux2*gz*(eta0(i,j)-eta0(i,j-1)) &
		    	        - aux2*(prd0(i,j,k)-prd0(i,j-1,k)))
	enddo
    	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		matrizgk(i,j,k) = dzz(i,j,k) *(Fw(i,j,k) -(1-tetah)*dt/dzz(i,j,k) * (prd0(i,j,k)-prd0(i,j,k-1)))
	enddo
    	enddo
	enddo
	!!! condição de contorno de topo (desativado) !!!
	!matrizgi(i,j,nz) = matrizgi(i,j,nz) ! + (cwind*uwind*dt)
	!matrizgj(i,j,nz) = matrizgj(i,j,nz) ! + (cwind*vwind*dt)

	!%%%!-- Matrizes A**-1 --!%%%!
	CALL matrizainv(matrizainvi,matrizainvj,matrizainvk)

	!%%%!--Cálculo das Partes (início do desnível)--!%%%!
	!reiniciando
	parte1i = 0.
	parte1  = 0.
	parte2  = 0.
	parte3j = 0.
	parte3  = 0.
	parte4  = 0.
	
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
	do kk = 1, nz
		parte1i(i,j,k) = parte1i(i,j,k) + dzx(i,j,kk) * matrizainvi(i,j,kk,k)
	enddo
		parte1(i,j) = parte1(i,j) + parte1i(i,j,k) * dzx(i,j,k)
		parte2(i,j) = parte2(i,j) + parte1i(i,j,k) * matrizgi(i,j,k)
	enddo
	enddo

	do j = 1, ny1
	do i = 1, nx
	do kk = 1, nz
		parte3j(i,j,k) = parte3j(i,j,k) + dzy(i,j,kk) * matrizainvj(i,j,kk,k)
	enddo
		parte3(i,j) = parte3(i,j) + parte3j(i,j,k) * dzy(i,j,k)
		parte4(i,j) = parte4(i,j) + parte3j(i,j,k) * matrizgj(i,j,k)
    	enddo
	enddo
	enddo

	! Cálculo de eta1 pelo gradiente conjudago (pressão hidrostática)
	CALL gradeta(parte1, parte2, parte3, parte4)

	! Cálculo das velocidades referentes a pressão hidrostática em "n+1"
	CALL poseta(matrizgi,matrizgj,matrizgk,matrizainvi,matrizainvj,matrizainvk)

!==================================================================================================================
END SUBROUTINE matrizg




