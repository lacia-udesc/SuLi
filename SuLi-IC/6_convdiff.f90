	!Subrotina para calcular o Fu, utilizando o ponto de vista lagrangiano
	!Referencia: Casulli (1990, 1992)

	!!! Implementação 15/04/2014
	! Leonardo Romero Monteiro

	!!! Modificações
	! Leonardo Romero Monteiro - 13/01/2015

	!!! Observações
	! no lagrangiano as velocidades no tempo n = t serão calculadas para o tempo n = t-1

	! são calculadas as Fu, Fv e Fw
	SUBROUTINE convdiff()

	USE velpre
	USE parametros
	USE smag
	USE ls_param
	USE vartempo
	USE mms_m
	USE obst
	IMPLICIT NONE

	!===================================================================================================================
	! auxiliares de velocidades: velocidades lagrangianas
	real(8), dimension(nx1,ny,nz) :: uint
	real(8), dimension(nx,ny1,nz) :: vint
	real(8), dimension(nx,ny,nz1) :: wint

	real(8), dimension(0:nx1,ny,nz) :: ama
	real(8), dimension(nx,0:ny1,nz) :: bmb
	real(8), dimension(nx,ny,0:nz1) :: dmd
	real(8), dimension(nx1,ny1,nz)  :: amb, bma
	real(8), dimension(nx1,ny,nz1)  :: amd, dma
	real(8), dimension(nx,ny1,nz1)  :: bmd, dmb

	!real(8), dimension(0:nx1,ny1,nz)  :: vx
	!real(8), dimension(0:nx1,ny,nz1)  :: wx

	!real(8), dimension(nx1,0:ny1,nz)  :: uy
	!real(8), dimension(nx,0:ny1,nz1)  :: wy

	!real(8), dimension(nx1,ny,0:nz1)  :: uz
	!real(8), dimension(nx,ny1,0:nz1)  :: vz

	real(8), dimension(nx1,ny,nz) :: rhox, dhsdx
	real(8), dimension(nx,ny1,nz) :: rhoy, dhsdy
	real(8), dimension(nx,ny,nz1) :: rhoz, dhsdz, epis_z

	!contadores
	integer :: i, j, k
	integer :: ntal
	real(8)    :: tal

	!auxiliares
	real(8) :: aux1, aux2

	!===================================================================================================================
	
	!reinicializando a rotina
	prd0 = prd1

	!reinicializando a subrotina
	epis_z = 0.

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

	!Camada esponja. Ativar somente quando necessário.
	!CALL sponge_layer(epis_z)

	!Cálculo dos Fs

	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,rhoz) !(nx,ny,nz1)

	call interpx_cf(hsx*kurv/rho,nx,ny,nz,dhsdx) !(nx1,ny,nz)
	call interpy_cf(hsy*kurv/rho,nx,ny,nz,dhsdy) !(nx,ny1,nz)
	call interpz_cf(hsz*kurv/rho,nx,ny,nz,dhsdz) !(nx,ny,nz1)

!	call interpy_cf(u(1:nx1,1:ny,1:nz),nx1,ny,nz,uy(1:nx1,1:ny1,1:nz)) !(nx1,ny1,nz)
!	call interpz_cf(u(1:nx1,1:ny,1:nz),nx1,ny,nz,uz(1:nx1,1:ny,1:nz1)) !(nx1,ny,nz1)


!	call interpx_cf(v(1:nx,1:ny1,1:nz),nx,ny1,nz,vx(1:nx1,1:ny1,1:nz)) !(nx1,ny1,nz)
!	call interpz_cf(v(1:nx,1:ny1,1:nz),nx,ny1,nz,vz(1:nx,1:ny1,1:nz1)) !(nx,ny1,nz1)

!	call interpx_cf(w(1:nx,1:ny,1:nz1),nx,ny,nz1,wx(1:nx1,1:ny,1:nz1)) !(nx1,ny,nz1)
!	call interpy_cf(w(1:nx,1:ny,1:nz1),nx,ny,nz1,wy(1:nx,1:ny1,1:nz1)) !(nx,ny1,nz1)

!      vx(0,:,:)=vx(1,:,:)
!      wx(0,:,:)=wx(1,:,:)

!      uy(:,0,:)=uy(:,1,:)
!      wy(:,0,:)=wy(:,1,:)

!      uz(:,:,0)=uz(:,:,1)
!      vz(:,:,0)=vz(:,:,1)


	call interpy_cf(xnut,nx1,ny,nz,bma) !(nx1,ny1,nz)
	call interpz_cf(xnut,nx1,ny,nz,dma) !(nx1,ny,nz1)

	call interpx_cf(ynut,nx,ny1,nz,amb) !(nx1,ny1,nz)
	call interpz_cf(ynut,nx,ny1,nz,dmb) !(nx,ny1,nz1)

	call interpx_cf(znut,nx,ny,nz1,amd) !(nx1,ny,nz1)
	call interpy_cf(znut,nx,ny,nz1,bmd) !(nx,ny1,nz1)

      ama(1:nx,1:ny,1:nz)=nut(1:nx,1:ny,1:nz)/dx
      bma=bma/dy
      dma=dma/dz

      amb=amb/dx
      bmb(1:nx,1:ny,1:nz)=nut(1:nx,1:ny,1:nz)/dy
      dmb=dmb/dz

      amd=amd/dx
      bmd=bmd/dy
      dmd(1:nx,1:ny,1:nz)=nut(1:nx,1:ny,1:nz)/dz

      ama(0,:,:)=ama(1,:,:)
      bmb(:,0,:)=bmb(:,1,:)
      dmd(:,:,0)=dmd(:,:,1)

      ama(nx1,:,:)=ama(nx,:,:)
      bmb(:,ny1,:)=bmb(:,ny,:)
      dmd(:,:,nz1)=dmd(:,:,nz)

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		Fu(i,j,k) = -uint(i,j,k) + ((ama(i,j,k)*2.*(u(i+1,j,k)-u(i,j,k))-&
                            ama(i-1,j,k)*2.*(u(i,j,k)-u(i-1,j,k)))/dx + &
                            bma(i,j+1,k)*((u(i,j+1,k)-u(i,j,k))/dy+(v(i,j+1,k)-v(i-1,j+1,k))/dx) - &
                            bma(i,j,k)  *((u(i,j,k)-u(i,j-1,k))/dy+(v(i,j,k)  -v(i-1,j,k))/dx  ) + &
                            dma(i,j,k+1)*((u(i,j,k+1)-u(i,j,k))/dz+(w(i,j,k+1)-w(i-1,j,k+1))/dx) - &
                            dma(i,j,k)*  ((u(i,j,k)-u(i,j,k-1))/dz+(w(i,j,k)  -w(i-1,j,k))/dx  ))/rhox(i,j,k) + gx &
			    -sigma*dhsdx(i,j,k) +tf_u(i,j,k) - gz/(chezy*chezy)*sqrt(ub(i,j,k)*ub(i,j,k))*u(i,j,k)
	enddo
	enddo	

	do j = 1, ny1
	do i = 1, nx
		Fv(i,j,k) = -vint(i,j,k) + (amb(i+1,j,k)*((v(i+1,j,k)-v(i,j,k))/dx+(u(i+1,j,k)-u(i+1,j-1,k))/dy) - &
                            		    amb(i,j,k)  *((v(i,j,k)-v(i-1,j,k))/dx+(u(i,j,k)  -u(i,j-1,k))/dy  ) + &
                            (bmb(i,j,k)*2.*(v(i,j+1,k)-v(i,j,k)) - &
                            bmb(i,j-1,k)*2.*(v(i,j,k)-v(i,j-1,k)))/dy + &
                            dmb(i,j,k+1)*((v(i,j,k+1)-v(i,j,k))/dz+(w(i,j,k+1)-w(i,j-1,k+1))/dy ) - &
                            dmb(i,j,k)  *((v(i,j,k)-v(i,j,k-1))/dz+(w(i,j,k)  -w(i,j-1,k))/dy   ))/rhoy(i,j,k) &
			    -sigma*dhsdy(i,j,k) +tf_v(i,j,k) - gz/(chezy*chezy)*sqrt(vb(i,j,k)*vb(i,j,k))*v(i,j,k)
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		Fw(i,j,k) = -wint(i,j,k) + (amd(i+1,j,k)*((w(i+1,j,k)-w(i,j,k))/dx+(u(i+1,j,k)-u(i+1,j,k-1))/dz) - &
                            		    amd(i,j,k)  *((w(i,j,k)-w(i-1,j,k))/dx+(u(i,j,k)  -u(i,j,k-1))/dz  ) + &
                            		    bmd(i,j+1,k)*((w(i,j+1,k)-w(i,j,k))/dy+(v(i,j+1,k)-v(i,j+1,k-1))/dz) - &
                            		    bmd(i,j,k)  *((w(i,j,k)-w(i,j-1,k))/dy+(v(i,j,k)  -v(i,j,k-1))/dz  ) + &
                            (dmd(i,j,k)*2.*(w(i,j,k+1)-w(i,j,k))-&
                            dmd(i,j,k-1)*2.*(w(i,j,k)-w(i,j,k-1)))/dz)/rhoz(i,j,k) &
			    - (w(i,j,k)-0.)* epis_z(i,j,k) - gz &
			    -sigma*dhsdz(i,j,k) +tf_w(i,j,k) - gz/(chezy*chezy)*sqrt(wb(i,j,k)*wb(i,j,k))*w(i,j,k)
	enddo
	enddo
	enddo

!==================================================================================================================
END SUBROUTINE convdiff





