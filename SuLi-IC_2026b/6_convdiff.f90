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
	USE param
	USE les
	USE vartempo
	USE cond, only : ccx0, ccy0
	!$ USE omp_lib

	IMPLICIT NONE

	!===================================================================================================================
	! auxiliares de velocidades: velocidades lagrangianas
	real(8),dimension(nx1,ny,nz) :: uint
	real(8),dimension(nx,ny1,nz) :: vint
	real(8),dimension(nx,ny,nz1) :: wint

	real(8),dimension(0:nx1,ny,nz) :: ama
	real(8),dimension(nx,0:ny1,nz) :: bmb
	real(8),dimension(nx,ny,0:nz1) :: dmd
	real(8),dimension(nx1,ny1,nz)  :: amb, bma
	real(8),dimension(nx1,ny,nz1)  :: amd, dma
	real(8),dimension(nx,ny1,nz1)  :: bmd, dmb

	!real(8),dimension(0:nx1,ny1,nz)  :: vx
	!real(8),dimension(0:nx1,ny,nz1)  :: wx

	!real(8),dimension(nx1,0:ny1,nz)  :: uy
	!real(8),dimension(nx,0:ny1,nz1)  :: wy

	!real(8),dimension(nx1,ny,0:nz1)  :: uz
	!real(8),dimension(nx,ny1,0:nz1)  :: vz

	real(8),dimension(nx1,ny,nz) :: rhox
	real(8),dimension(nx,ny1,nz) :: rhoy
	real(8),dimension(nx,ny,nz1) :: rhoz

	!contadores
	integer :: i, j, k

	!auxiliares
	real(8) :: aux1, aux2

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================

	!%%%!-- Advectivo --!%%%!
	call classico(uint,vint,wint)

	!Cálculo dos Fs

	call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
	call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
	call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)


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


	call interpy_cf(xmu,nx1,ny,nz,bma) !(nx1,ny1,nz)
	call interpz_cf(xmu,nx1,ny,nz,dz,dma) !(nx1,ny,nz1)

	call interpx_cf(ymu,nx,ny1,nz,amb) !(nx1,ny1,nz)
	call interpz_cf(ymu,nx,ny1,nz,dz,dmb) !(nx,ny1,nz1)

	call interpx_cf(zmu,nx,ny,nz1,amd) !(nx1,ny,nz1)
	call interpy_cf(zmu,nx,ny,nz1,bmd) !(nx,ny1,nz1)


        !$omp parallel
        !$omp do
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx
                ama(i,j,k) = mu(i,j,k) *inv_dx
                bmb(i,j,k) = mu(i,j,k) *inv_dy
                dmd(i,j,k) = mu(i,j,k) *inv_dzz(k)
                dma(i,j,k) = dma(i,j,k)*inv_dz(k)
                dmb(i,j,k) = dmb(i,j,k)*inv_dz(k)
        enddo
        enddo
        enddo
        !$omp end do

        !$omp do
        do k = 1, nz
        do j = 1, ny
                if (ccx0 == 0) then
                        ! wrap periodico: celula 0 e imagem da celula nx
                        ! e a celula nx1 e imagem da celula 1
                        ama(0,j,k)   = ama(nx,j,k)
                        ama(nx1,j,k) = ama(1,j,k)
                else
                        ama(0,j,k)   = ama(1,j,k)
                        ama(nx1,j,k) = ama(nx,j,k)
                endif
        enddo
        enddo
        !$omp end do nowait
        
        !$omp do
        do k = 1, nz
        do i = 1, nx
                if (ccy0 == 0) then
                	bmb(i,0,k) = bmb(i,ny,k)
                	bmb(i,ny1,k) = bmb(i,1,k)
                else
                	bmb(i,0,k) = bmb(i,1,k)
                	bmb(i,ny1,k) = bmb(i,ny,k)
                endif
        enddo
        enddo
        !$omp end do
	!$omp end parallel



	!$omp parallel do
	do j = 1, ny
	do i = 1, nx
		dmd(i,j,0) = dmd(i,j,1)
		dmd(i,j,nz1) = dmd(i,j,nz)
	enddo
	enddo
	!$omp end parallel do

	! Escalas (evita atribuições de array fora do OpenMP)
	!$omp parallel
	!$omp do
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx1
		amb(i,j,k) = amb(i,j,k) * inv_dx
		bma(i,j,k) = bma(i,j,k) * inv_dy
	enddo
	enddo
	enddo
	!$omp end do nowait

	!$omp do
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx1
		amd(i,j,k) = amd(i,j,k) * inv_dx
	enddo
	enddo
	enddo
	!$omp end do nowait

	!$omp do
	do k = 1, nz1
	do j = 1, ny1
	do i = 1, nx
		bmd(i,j,k) = bmd(i,j,k) * inv_dy
	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel


		!$omp parallel
		!$omp do collapse(3)
        do k = 1, nz
        do j = 1, ny
        do i = 1, nx1

                Fu(i,j,k) = -uint(i,j,k) + ((ama(i,j,k)*2.*(u(i+1,j,k)-u(i,j,k)) - &
                                    ama(i-1,j,k)*2.*(u(i,j,k)-u(i-1,j,k)))*inv_dx + &
                                    bma(i,j+1,k)*((u(i,j+1,k)-u(i,j,k))*inv_dy+(v(i,j+1,k)-v(i-1,j+1,k))*inv_dx) - &
                                    bma(i,j,k)  *((u(i,j,k)-u(i,j-1,k))*inv_dy+(v(i,j,k)  -v(i-1,j,k))*inv_dx  ) + &
                                    dma(i,j,k+1)*((u(i,j,k+1)-u(i,j,k))*inv_dz(k)+(w(i,j,k+1)-w(i-1,j,k+1))*inv_dx) - &
                                    dma(i,j,k)*  ((u(i,j,k)-u(i,j,k-1))*inv_dz(k)+(w(i,j,k)  -w(i-1,j,k))*inv_dx  ))/rhox(i,j,k)
        enddo
        enddo
        enddo
        !$omp end do nowait

		!$omp do collapse(3)
        do k = 1, nz
        do j = 1, ny1
        do i = 1, nx

                Fv(i,j,k) = -vint(i,j,k) + (amb(i+1,j,k)*((v(i+1,j,k)-v(i,j,k))*inv_dx+(u(i+1,j,k)-u(i+1,j-1,k))*inv_dy) - &
                                            amb(i,j,k)  *((v(i,j,k)-v(i-1,j,k))*inv_dx+(u(i,j,k)  -u(i,j-1,k))*inv_dy  ) + &
                                    (bmb(i,j,k)*2.*(v(i,j+1,k)-v(i,j,k)) - &
                                    bmb(i,j-1,k)*2.*(v(i,j,k)-v(i,j-1,k)))*inv_dy + &
                                    dmb(i,j,k+1)*((v(i,j,k+1)-v(i,j,k))*inv_dz(k)+(w(i,j,k+1)-w(i,j-1,k+1))*inv_dy ) - &
                                    dmb(i,j,k)  *((v(i,j,k)-v(i,j,k-1))*inv_dz(k)+(w(i,j,k)  -w(i,j-1,k))*inv_dy   ))/rhoy(i,j,k)
        enddo
        enddo
        enddo
        !$omp end do nowait

		!$omp do collapse(3)
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx


			Fw(i,j,k) = -wint(i,j,k) + (amd(i+1,j,k)*((w(i+1,j,k)-w(i,j,k))*inv_dx+(u(i+1,j,k)-u(i+1,j,k-1))*inv_dz(k)) - &
											amd(i,j,k)  *((w(i,j,k)-w(i-1,j,k))*inv_dx+(u(i,j,k)  -u(i,j,k-1))*inv_dz(k)  ) + &
											bmd(i,j+1,k)*((w(i,j+1,k)-w(i,j,k))*inv_dy+(v(i,j+1,k)-v(i,j+1,k-1))*inv_dz(k)) - &
											bmd(i,j,k)  *((w(i,j,k)-w(i,j-1,k))*inv_dy+(v(i,j,k)  -v(i,j,k-1))*inv_dz(k)  ) + &
								(dmd(i,j,k)*2.*(w(i,j,k+1)-w(i,j,k))-&
								dmd(i,j,k-1)*2.*(w(i,j,k)-w(i,j,k-1)))*inv_dzz(k))/rhoz(i,j,k)
		enddo
		enddo
		enddo
		!$omp end do
		!$omp end parallel

!==================================================================================================================
END SUBROUTINE convdiff

!##################################################################################################################

SUBROUTINE convdiff_les(Fka) ! Termos convectivos e difusivos da energia cinética

	USE velpre
	USE param
	USE les
	USE vartempo
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k
	real(8), dimension(nx,ny,nz) :: muaux, Fka, adv

	real(8),dimension(nx1,ny,nz) :: am
	real(8),dimension(nx,ny1,nz) :: bm
	real(8),dimension(nx,ny,nz1) :: dm

	call adv_esc(ka,nx,ny,nz,adv)  ! Termo advectivo
	

!	!$omp parallel do
!	do k = 1, nz
!	do j = 1, ny
!	do i = 1, nx
!	        muaux(i,j,k) = ls_mu(i,j,k) + mut(i,j,k)!*rho(i,j,k)
!	enddo
!	enddo
!	enddo
!	!$omp end parallel do


	call interpx_cf(mu,nx,ny,nz,am)
	call interpy_cf(mu,nx,ny,nz,bm)
	call interpz_cf(mu,nx,ny,nz,dz,dm)

	!$omp parallel do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		Fka(i,j,k) = - adv(i,j,k) + &
			      ((am(i+1,j,k)*(ka(i+1,j,k)-ka(i,j,k))-am(i,j,k)*(ka(i,j,k)-ka(i-1,j,k)))*inv_dx2 + &
			      (bm(i,j+1,k)*(ka(i,j+1,k)-ka(i,j,k))-bm(i,j,k)*(ka(i,j,k)-ka(i,j-1,k)))*inv_dy2  + &
			      (dm(i,j,k+1)*(ka(i,j,k+1)-ka(i,j,k))-dm(i,j,k)*(ka(i,j,k)-ka(i,j,k-1)))*inv_dz2(k)) /rho(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

END SUBROUTINE convdiff_les


!##################################################################################################################

SUBROUTINE convdiff_sst_k(Fk, sigma_k)
	! Termos convectivos e difusivos para k (SST): mu_eff = mu_lam + sigma_k*mu_t

	USE velpre
	USE param
	USE les
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k
	real(8), intent(in), dimension(nx,ny,nz) :: sigma_k
	real(8), dimension(nx,ny,nz) :: muaux, Fk, adv

	real(8),dimension(nx1,ny,nz) :: am
	real(8),dimension(nx,ny1,nz) :: bm
	real(8),dimension(nx,ny,nz1) :: dm

	call adv_esc(ka,nx,ny,nz,adv)

	!$omp parallel do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		muaux(i,j,k) = ls_mu(i,j,k) + sigma_k(i,j,k)*mut(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

	call interpx_cf(muaux,nx,ny,nz,am)
	call interpy_cf(muaux,nx,ny,nz,bm)
	call interpz_cf(muaux,nx,ny,nz,dz,dm)

	!$omp parallel do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		Fk(i,j,k) = - adv(i,j,k) + &
		      ((am(i+1,j,k)*(ka(i+1,j,k)-ka(i,j,k))-am(i,j,k)*(ka(i,j,k)-ka(i-1,j,k)))*inv_dx2 + &
		      (bm(i,j+1,k)*(ka(i,j+1,k)-ka(i,j,k))-bm(i,j,k)*(ka(i,j,k)-ka(i,j-1,k)))*inv_dy2  + &
		      (dm(i,j,k+1)*(ka(i,j,k+1)-ka(i,j,k))-dm(i,j,k)*(ka(i,j,k)-ka(i,j,k-1)))*inv_dz2(k)) /rho(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

END SUBROUTINE convdiff_sst_k


!##################################################################################################################

SUBROUTINE convdiff_sst_omega(Fom, sigma_om)
	! Termos convectivos e difusivos para omega (SST): mu_eff = mu_lam + sigma_om*mu_t

	USE velpre
	USE param
	USE les
	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k
	real(8), intent(in), dimension(nx,ny,nz) :: sigma_om
	real(8), dimension(nx,ny,nz) :: muaux, Fom, adv

	real(8),dimension(nx1,ny,nz) :: am
	real(8),dimension(nx,ny1,nz) :: bm
	real(8),dimension(nx,ny,nz1) :: dm

	call adv_esc(omg,nx,ny,nz,adv)

	!$omp parallel do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		muaux(i,j,k) = ls_mu(i,j,k) + sigma_om(i,j,k)*mut(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

	call interpx_cf(muaux,nx,ny,nz,am)
	call interpy_cf(muaux,nx,ny,nz,bm)
	call interpz_cf(muaux,nx,ny,nz,dz,dm)

	!$omp parallel do collapse(3)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		Fom(i,j,k) = - adv(i,j,k) + &
		      ((am(i+1,j,k)*(omg(i+1,j,k)-omg(i,j,k))-am(i,j,k)*(omg(i,j,k)-omg(i-1,j,k)))*inv_dx2 + &
		      (bm(i,j+1,k)*(omg(i,j+1,k)-omg(i,j,k))-bm(i,j,k)*(omg(i,j,k)-omg(i,j-1,k)))*inv_dy2  + &
		      (dm(i,j,k+1)*(omg(i,j,k+1)-omg(i,j,k))-dm(i,j,k)*(omg(i,j,k)-omg(i,j,k-1)))*inv_dz2(k)) /rho(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

END SUBROUTINE convdiff_sst_omega



