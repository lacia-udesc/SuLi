!Subrotina para calcular a pressão dinâmica

!Versão: CG
!Referencia: Casulli (1992 e 1999)
!!! Implementação original 15/04/2014

!Versão: GMRES(3) - Generalized Minimal Residual com 3 direções acumuladas
!Referencia: Saad & Schultz (1986)
!!! Implementação GMRES com múltiplas direções - 28/01/2026

!Leonardo Romero Monteiro

module fft_zcoef

	! Coeficientes do operador em z da solução de Poisson por FFT/tridiagonal.
	! Dependem apenas de dz/dzz, então são calculados uma única vez.

	USE disc, only: nz, dz, dzz
	implicit none
	real(8), allocatable :: azL(:), azD(:), azU(:)

contains

	subroutine init_fft_zcoef()
		implicit none
		integer :: k
		allocate(azL(nz), azD(nz), azU(nz))
		azL = 0.0d0
		azD = 0.0d0
		azU = 0.0d0
		! k = 1
		azU(1) = 1.0d0/(dz(1)*dzz(2))
		azD(1) = -azU(1)
		! 1 < k < nz
		do k = 2, nz-1
			azL(k) = 1.0d0/(dz(k)*dzz(k))
			azU(k) = 1.0d0/(dz(k)*dzz(k+1))
			azD(k) = -(azL(k) + azU(k))
		enddo
		! k = nz
		azL(nz) = 1.0d0/(dz(nz)*dzz(nz))
		azD(nz) = -azL(nz)
	end subroutine init_fft_zcoef

end module fft_zcoef


SUBROUTINE graddin()

	USE velpre
	USE obst
	USE ls_param, only: rho_f1
	!$ USE omp_lib
	
	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8) :: alfa, alfamu, alfadi, beta, betamu
	real(8),dimension(nx1,ny,nz) :: matsi
	real(8),dimension(nx,ny1,nz) :: matsj
	real(8),dimension(nx,ny,nz1) :: matsk
	real(8),dimension(nx,ny,nz) :: matq
	real(8),dimension(0:nx1,0:ny1,0:nz1) :: matd, mate, erro, errop

	!contadores
	integer :: i, j, k, cont, conv_press

	!auxiliares
	real(8) :: aux1, aux2

	real(8),dimension(nx1,ny,nz) :: rhox
	real(8),dimension(nx,ny1,nz) :: rhoy
	real(8),dimension(nx,ny,nz1) :: rhoz

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA
	!===================================================================================================================
	CALL divergente(matq)

	!%%%!-- Método do Gradiente Conjugado - Para Pressão Dinâmica --!%%%!

	if  (t_press .eq. 1) then !Projection Method
	
		call interpx_cf(rho,nx,ny,nz,rhox) !(nx1,ny,nz)
		call interpy_cf(rho,nx,ny,nz,rhoy) !(nx,ny1,nz)
		call interpz_cf(rho,nx,ny,nz,dz,rhoz) !(nx,ny,nz1)

	else !Implicit Pressure Algorithm (IPA) 

	! -----------------------------------------------------------------------------------------------------------------
	! FFT Poisson (apenas para IPA / rho constante): resolve (dt/rho) * Laplaciano(prd1) = divergente
	! Implementação via FFTW3 com extensão par (even extension) para Neumann.
	! Requer FFTW3 (fftw3.f03 + libfftw3)
	! -----------------------------------------------------------------------------------------------------------------

		! IBM: anula a fonte de divergência no interior do sólido (obs_lss > 0),
		! para que o Poisson não gere pressões espúrias dentro do obstáculo.
		if (obst_t .ne. 0 .and. (ibm_t == 1 .or. ibm_t == 2)) then
			!$omp parallel do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				if (obs_lss(i,j,k) > 0.) matq(i,j,k) = 0.
			enddo
			enddo
			enddo
			!$omp end parallel do
		endif

		CALL poisson_fft3d_const_rho(prd1, matq)

		! IBM: efetivação do obstáculo como contorno sólido. O solucionador
		! espectral não enxerga a geometria do obstáculo, então a condição de
		! parede é aplicada a posteriori: pressão nula no interior do sólido e
		! extensão de Neumann (média dos vizinhos fluidos) na 1ª camada de
		! células sólidas, anulando o gradiente de pressão através da interface
		! (dP/dn = 0) e evitando correção de velocidade espúria em posdin().
		if (obst_t .ne. 0 .and. (ibm_t == 1 .or. ibm_t == 2)) then
			CALL ibm_press_contorno(prd1)
		endif
		return

	endif

	aux1 = dt / (dx*dx)
	aux2 = dt / (dy*dy)

    !$omp parallel
    !$omp do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		matsi(i,j,k) = aux1/rhox(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		matsj(i,j,k) = aux2/rhoy(i,j,k)
	enddo
	enddo
	enddo
	!$omp end do nowait

    !$omp do
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		matsk(i,j,k) = dt /(dz(k)*dz(k)*rhoz(i,j,k))
	enddo
	enddo
	enddo
	!$omp end do
	!$omp end parallel

	! --------------------------------------------------------------------------
	! IBM: efetivação do obstáculo como contorno sólido na equação de Poisson.
	! - Coeficientes das faces que cruzam/tocam o sólido são anulados, impondo
	!   gradiente de pressão normal nulo (dP/dn = 0) na interface imersa;
	! - Células sólidas têm a fonte (matq) e a pressão (prd1) zeradas; em norm()
	!   a diagonal correspondente é ajustada para 1, removendo o nó do sistema.
	! --------------------------------------------------------------------------
	if (obst_t .ne. 0 .and. (ibm_t == 1 .or. ibm_t == 2)) then
		CALL ibm_mask_poisson(matsi, matsj, matsk, matq, prd1)
	endif

	! aplicação da normalização (pré-condicionamento global)
	call norm(prd1, mate, matd, matsi, matsj, matsk,matq,1)


    !$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		erro(i,j,k) = - mate(i,j,k) + matsi(i+1,j,k) * mate(i+1,j,k) + matsi(i,j,k) * mate(i-1,j,k) &
			+ matsj(i,j+1,k) * mate(i,j+1,k) + matsj(i,j,k) * mate(i,j-1,k) &
			+ matsk(i,j,k+1) * mate(i,j,k+1) + matsk(i,j,k) * mate(i,j,k-1) + matq(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

	CALL cc_geral(erro,nx,ny,nz ,1)	
	
	!! para quando fazer o precondicionamento interno	
	CALL soma_vec(erro,erro,nx,ny,nz,alfamu)

    !$omp parallel do
	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, nx+1
		errop(i,j,k) = erro(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do	

	cont = 0
	aux1 = dabs(alfamu)
	aux2 = (0.01/(nx*ny*nz))

	CALL cg(mate, matsi, matsj, matsk, erro, errop, cont, aux1, aux2, alfamu)

	! se pular pressão ele vai avisar
	! write(*,*) "o número de iterações foi; ", cont
	if (cont == 10000) write(*,*) "pulou pressão; ", "erro =", aux1
	contp = contp + cont
		
	!%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	! Desnormalização de matriz e para a pressão dinâmica do precondicionamento global
	call norm(prd1, mate, matd, matsi, matsj, matsk,matq,2)

END SUBROUTINE graddin

!===============================================================================================================
! IBM: efetivação do obstáculo como contorno sólido na equação de Poisson (CG / Projection Method)
!
! A interface imersa passa a atuar como parede:
! 1) os coeficientes das faces que cruzam o sólido (ou tocam uma célula sólida)
!    são anulados => condição de Neumann (dP/dn = 0) na parede;
! 2) as células sólidas têm a fonte de divergência (matq) e a pressão (prd1)
!    zeradas; em norm() a diagonal correspondente é ajustada para 1, de modo
!    que o nó sólido se torna identidade e fica desacoplado do sistema linear.
!===============================================================================================================

SUBROUTINE ibm_mask_poisson(matsi, matsj, matsk, matq, prd1)

	USE disc, only: nx, ny, nz, nx1, ny1, nz1
	USE obst, only: obs_lsx, obs_lsy, obs_lsz, obs_lss
	!$ USE omp_lib

	IMPLICIT NONE

	real(8),intent(inout),dimension(nx1,ny,nz) :: matsi
	real(8),intent(inout),dimension(nx,ny1,nz) :: matsj
	real(8),intent(inout),dimension(nx,ny,nz1) :: matsk
	real(8),intent(inout),dimension(nx,ny,nz)  :: matq
	real(8),intent(inout),dimension(0:nx+1,0:ny+1,0:nz+1) :: prd1

	integer :: i, j, k
	logical :: lesq, ldir

	! faces x (entre as células i-1 e i)
	!$omp parallel do private(lesq,ldir)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		lesq = (i > 1  .and. obs_lss(i-1,j,k) > 0.d0)
		ldir = (i <= nx .and. obs_lss(i  ,j,k) > 0.d0)
		if (obs_lsx(i,j,k) >= 0.d0 .or. lesq .or. ldir) matsi(i,j,k) = 0.d0
	enddo
	enddo
	enddo
	!$omp end parallel do

	! faces y (entre as células j-1 e j)
	!$omp parallel do private(lesq,ldir)
	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		lesq = (j > 1  .and. obs_lss(i,j-1,k) > 0.d0)
		ldir = (j <= ny .and. obs_lss(i,j  ,k) > 0.d0)
		if (obs_lsy(i,j,k) >= 0.d0 .or. lesq .or. ldir) matsj(i,j,k) = 0.d0
	enddo
	enddo
	enddo
	!$omp end parallel do

	! faces z (entre as células k-1 e k)
	!$omp parallel do private(lesq,ldir)
	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		lesq = (k > 1  .and. obs_lss(i,j,k-1) > 0.d0)
		ldir = (k <= nz .and. obs_lss(i,j,k  ) > 0.d0)
		if (obs_lsz(i,j,k) >= 0.d0 .or. lesq .or. ldir) matsk(i,j,k) = 0.d0
	enddo
	enddo
	enddo
	!$omp end parallel do

	! células sólidas: fonte de divergência e pressão zeradas
	!$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		if (obs_lss(i,j,k) > 0.d0) then
			matq(i,j,k) = 0.d0
			prd1(i,j,k) = 0.d0
		endif
	enddo
	enddo
	enddo
	!$omp end parallel do

END SUBROUTINE ibm_mask_poisson

!===============================================================================================================
! IBM: efetivação do obstáculo como contorno sólido no Poisson FFT (IPA)
!
! O solucionador espectral não enxerga a geometria do obstáculo; a condição de
! parede é aplicada a posteriori:
! 1) pressão zerada no interior do sólido;
! 2) na 1ª camada de células sólidas (adjacentes à interface) a pressão recebe
!    a média dos vizinhos fluidos => gradiente de pressão normal nulo através
!    da parede (dP/dn = 0), evitando correção de velocidade espúria em posdin().
!===============================================================================================================

SUBROUTINE ibm_press_contorno(prd1)

	USE disc, only: nx, ny, nz
	USE obst, only: obs_lss
	!$ USE omp_lib

	IMPLICIT NONE

	real(8),intent(inout),dimension(0:nx+1,0:ny+1,0:nz+1) :: prd1

	integer :: i, j, k, nviz
	real(8) :: soma
	real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: p_ref

	! cópia do campo resolvido (a extrapolação não depende da ordem de varredura)
	!$omp parallel do
	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, nx+1
		p_ref(i,j,k) = prd1(i,j,k)
	enddo
	enddo
	enddo
	!$omp end parallel do

	!$omp parallel do private(soma,nviz)
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		if (obs_lss(i,j,k) > 0.d0) then
			soma = 0.d0
			nviz = 0
			if (i > 1  .and. obs_lss(i-1,j,k) <= 0.d0) then; soma = soma + p_ref(i-1,j,k); nviz = nviz + 1; endif
			if (i < nx .and. obs_lss(i+1,j,k) <= 0.d0) then; soma = soma + p_ref(i+1,j,k); nviz = nviz + 1; endif
			if (j > 1  .and. obs_lss(i,j-1,k) <= 0.d0) then; soma = soma + p_ref(i,j-1,k); nviz = nviz + 1; endif
			if (j < ny .and. obs_lss(i,j+1,k) <= 0.d0) then; soma = soma + p_ref(i,j+1,k); nviz = nviz + 1; endif
			if (k > 1  .and. obs_lss(i,j,k-1) <= 0.d0) then; soma = soma + p_ref(i,j,k-1); nviz = nviz + 1; endif
			if (k < nz .and. obs_lss(i,j,k+1) <= 0.d0) then; soma = soma + p_ref(i,j,k+1); nviz = nviz + 1; endif
			if (nviz > 0) then
				prd1(i,j,k) = soma / dble(nviz)
			else
				prd1(i,j,k) = 0.d0
			endif
		endif
	enddo
	enddo
	enddo
	!$omp end parallel do

	CALL cc_geral(prd1,nx,ny,nz,1)

END SUBROUTINE ibm_press_contorno

!===============================================================================================================


SUBROUTINE norm(prda, mate, matd, matsi, matsj, matsk,matq,nlock)

	USE disc
	USE cond
	USE obst, only: obs_lss
		
	!$ USE omp_lib
	
	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8),intent(inout),dimension(nx1,ny,nz) :: matsi
	real(8),intent(inout),dimension(nx,ny1,nz) :: matsj
	real(8),intent(inout),dimension(nx,ny,nz1) :: matsk

	real(8),intent(inout),dimension(nx,ny,nz) :: matq
	
	real(8),intent(inout),dimension(0:nx1,0:ny1,0:nz1) :: mate, prda, matd

	!real(16), dimension(0:nx1,0:ny1,0:nz1) :: inv_sqrt_d !!este só funciona se estiver em real(16)

	!contadores
	integer :: i, j, k, nlock

	if (nlock == 1) then

		!! início - normalização

		!$omp parallel do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			matd(i,j,k) = matsi(i+1,j,k) + matsi(i,j,k) + matsj(i,j+1,k) + matsj(i,j,k) + matsk(i,j,k+1) + matsk(i,j,k)
		enddo
		enddo
		enddo
		!$omp end parallel do

		! IBM: células sólidas têm a diagonal ajustada para 1 (identidade).
		! Como as faces que tocam o sólido foram anuladas em ibm_mask_poisson(),
		! a soma dos coeficientes seria nula e a normalização dividiria por zero.
		! Com diagonal unitária o nó sólido fica desacoplado do sistema e a
		! pressão (zerada) permanece nula durante as iterações do CG.
		if (obst_t .ne. 0 .and. (ibm_t == 1 .or. ibm_t == 2)) then
			!$omp parallel do
			do k = 1, nz
			do j = 1, ny
			do i = 1, nx
				if (obs_lss(i,j,k) > 0.d0) matd(i,j,k) = 1.d0
			enddo
			enddo
			enddo
			!$omp end parallel do
		endif

		
		! Condições de contorno, von Neumann
		CALL cc_geral(matd,nx,ny,nz ,1)	
			
		! Normalização da pressão dinâmica
		! Altera dados

		!$omp parallel
		!$omp do
		do k = 0, nz1
		do j = 0, ny1
		do i = 0, nx1
			mate(i,j,k) = prda(i,j,k) * dsqrt(matd(i,j,k))
			!inv_sqrt_d(i,j,k) = 1.0d0 / dsqrt(matd(i,j,k))
		enddo
		enddo
		enddo
		!$omp end do nowait

		!Normalização das matrizes s e cálculo do primeiro erro

		!$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx1
			matsi(i,j,k) = matsi(i,j,k)/dsqrt(matd(i,j,k)*matd(i-1,j,k))
			!matsi(i,j,k) = matsi(i,j,k) * inv_sqrt_d(i,j,k) * inv_sqrt_d(i-1,j,k)
		enddo
		enddo
		enddo
		!$omp end do nowait

		!$omp do
		do k = 1, nz
		do j = 1, ny1
		do i = 1, nx
			matsj(i,j,k) = matsj(i,j,k)/dsqrt(matd(i,j,k)*matd(i,j-1,k))
			!matsj(i,j,k) = matsj(i,j,k) * inv_sqrt_d(i,j,k) * inv_sqrt_d(i,j-1,k)
		enddo
		enddo
		enddo
		!$omp end do nowait

		!$omp do
		do k = 1, nz1
		do j = 1, ny
		do i = 1, nx
			matsk(i,j,k) = matsk(i,j,k)/dsqrt(matd(i,j,k)*matd(i,j,k-1))
			!matsk(i,j,k) = matsk(i,j,k) * inv_sqrt_d(i,j,k) * inv_sqrt_d(i,j,k-1)
		enddo
		enddo
		enddo
		!$omp end do nowait

		!$omp do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			matq(i,j,k) = matq(i,j,k)/dsqrt(matd(i,j,k))
			!matq(i,j,k) = matq(i,j,k) * inv_sqrt_d(i,j,k)
		enddo
		enddo
		enddo
		!$omp end do
		!$omp end parallel

		!! fim - précondicionamento

	else

		! Condições de contorno
      	CALL cc_geral(mate,nx,ny,nz ,1)	
	
		!! volta do précondicioamento

		!$omp parallel do
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, nx+1
			prda(i,j,k) = mate(i,j,k) / dsqrt(matd(i,j,k))
		enddo
		enddo
		enddo
		!$omp end parallel do


	endif


END SUBROUTINE norm

!===================================================================================================================



SUBROUTINE soma_vec(a,b,dimx,dimy,dimz,soma)

	USE disc, only: omp_t
	!$ USE omp_lib
	
	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8),intent(in),dimension(0:dimx+1,0:dimy+1,0:dimz+1) :: a,b

	real(8),intent(out) :: soma
	integer,intent(in) :: dimx, dimy, dimz
	integer :: 	i, j, k

	soma = 0.d0

	if (omp_t == 0) then
		do k = 1, dimz
		do j = 1, dimy
		do i = 1, dimx
			soma = soma + a(i,j,k) * b(i,j,k)
		enddo
		enddo
		enddo
	else
		!$omp parallel do reduction(+:soma)
		do k = 1, dimz
		do j = 1, dimy
		do i = 1, dimx
			soma = soma + a(i,j,k) * b(i,j,k)
		enddo
		enddo
		enddo
        !$omp end parallel do 
	endif
	
END SUBROUTINE soma_vec



!!!!#########!!!!!!!!!##################!!!!!!!!!!#########!!!!
!!!!#########!!!!!!!!!		CG - INÍCIO	!!!!!!!!!!#########!!!!
!!!!#########!!!!!!!!!##################!!!!!!!!!!#########!!!!

SUBROUTINE cg(mate, matsi, matsj, matsk, erro, errop, cont, aux1, aux2, alfamu)

	USE disc, only: nx, ny, nz, nx1, ny1, nz1, pcg_prec_t
	!$ USE omp_lib
	
	IMPLICIT NONE

	!===================================================================================================================
	!DECLARADO SOMENTE NA SUBROTINA (ou não precisam de entrada)

	real(8) :: alfa, alfamu, alfadi, beta, betamu
	real(8) :: rz_old, rz_new
	real(8),dimension(nx1,ny,nz) :: matsi
	real(8),dimension(nx,ny1,nz) :: matsj
	real(8),dimension(nx,ny,nz1) :: matsk
	real(8),dimension(0:nx1,0:ny1,0:nz1) :: mate, erro, errop, mp, zpc

	!contadores
	integer :: i, j, k, cont

	!auxiliares
	real(8) :: aux1, aux2

	!===================================================================================================================
	!RESOLUÇÃO DO PROBLEMA (PCG com precondicionador selecionável)
	!============
	integer :: pcg_prec
	pcg_prec = pcg_prec_t
	if (pcg_prec < 0 .or. pcg_prec > 3) pcg_prec = 0

	! residual norm para critério de parada
	CALL soma_vec(erro,erro,nx,ny,nz,alfamu)
	aux1 = dabs(alfamu)

	! z = M^{-1} r
	CALL precond_apply(pcg_prec, erro, zpc, matsi, matsj, matsk)
	CALL soma_vec(erro,zpc,nx,ny,nz,rz_old)

	! p = z
    !$omp parallel do
	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, nx+1
		errop(i,j,k) = zpc(i,j,k)
	enddo
	enddo
	enddo
    !$omp end parallel do
	CALL cc_geral(errop,nx,ny,nz ,1)

	do while ((aux1 > aux2) .and. (cont < 10000) )
		cont = cont + 1

		! Ap = A p  (A = I - offdiag)
        !$omp parallel do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			mp(i,j,k) = errop(i,j,k) - errop(i+1,j,k)*matsi(i+1,j,k) - errop(i-1,j,k)*matsi(i,j,k) &
			- errop(i,j+1,k)*matsj(i,j+1,k) - errop(i,j-1,k)*matsj(i,j,k) &
			- errop(i,j,k+1)*matsk(i,j,k+1) - errop(i,j,k-1)*matsk(i,j,k)
		enddo
		enddo
		enddo
        !$omp end parallel do

		CALL soma_vec(errop,mp,nx,ny,nz,alfadi)

		alfa = rz_old / alfadi

        !$omp parallel do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			mate(i,j,k) = mate(i,j,k) + alfa * errop(i,j,k)
			erro(i,j,k) = erro(i,j,k) - alfa * mp(i,j,k)
		enddo
		enddo
		enddo
        !$omp end parallel do

		CALL soma_vec(erro,erro,nx,ny,nz,alfamu)
		aux1 = dabs(alfamu)

		! z = M^{-1} r
		CALL precond_apply(pcg_prec, erro, zpc, matsi, matsj, matsk)
		CALL soma_vec(erro,zpc,nx,ny,nz,rz_new)

		beta = rz_new / rz_old
		rz_old = rz_new

        !$omp parallel do
		do k = 1, nz
		do j = 1, ny
		do i = 1, nx
			errop(i,j,k) = zpc(i,j,k) + beta * errop(i,j,k)
		enddo
		enddo
		enddo
        !$omp end parallel do
		CALL cc_geral(errop,nx,ny,nz ,1)

	enddo

	return

END SUBROUTINE cg


SUBROUTINE precond_apply(prec_t, r, zpc, matsi, matsj, matsk)

	USE disc, only: nx, ny, nz, nx1, ny1, nz1

	IMPLICIT NONE

	integer,intent(in) :: prec_t
	integer :: i, j, k
	real(8),intent(in),dimension(nx1,ny,nz) :: matsi
	real(8),intent(in),dimension(nx,ny1,nz) :: matsj
	real(8),intent(in),dimension(nx,ny,nz1) :: matsk
	real(8),intent(in),dimension(0:nx1,0:ny1,0:nz1) :: r
	real(8),intent(out),dimension(0:nx1,0:ny1,0:nz1) :: zpc


	if (prec_t == 0) then
		! cópia explícita
		!$omp parallel do
		do k = 0, nz+1
		do j = 0, ny+1
		do i = 0, nx+1
			zpc(i,j,k) = r(i,j,k)
		enddo
		enddo
		enddo
		!$omp end parallel do
	elseif (prec_t == 1) then
		CALL precond_sgs(r, zpc, matsi, matsj, matsk)
	else
		write(*,*) 'ERRO: precondicionador PCG inválido.'
		stop
	endif

END SUBROUTINE precond_apply

!######################################################################################


SUBROUTINE precond_sgs(r, zpc, matsi, matsj, matsk)

	USE disc, only: nx, ny, nz, nx1, ny1, nz1

	IMPLICIT NONE

	real(8),intent(in),dimension(nx1,ny,nz) :: matsi
	real(8),intent(in),dimension(nx,ny1,nz) :: matsj
	real(8),intent(in),dimension(nx,ny,nz1) :: matsk
	real(8),intent(in),dimension(0:nx1,0:ny1,0:nz1) :: r
	real(8),intent(out),dimension(0:nx1,0:ny1,0:nz1) :: zpc

	real(8),dimension(0:nx1,0:ny1,0:nz1) :: y_sgs
	integer :: i, j, k

	! SGS red-black (2 cores) consistente com a ordenação por blocos (red depois black).
	! Para stencil 7-pontos, não há acoplamento entre nós da mesma cor; então o bloco-diagonal é identidade.
	!
	! (D+L) y = r  com L apenas no bloco (black <- red):
	!   y_red   = r_red
	!   y_black = r_black + Sum(s * y_red(vizinhos))
	!
	! (D+U) z = y  com U apenas no bloco (red <- black):
	!   z_black = y_black
	!   z_red   = y_red + Sum(s * z_black(vizinhos))

	! inicialização explícita (evita atribuição de array fora do OpenMP)
	!$omp parallel do
	do k = 0, nz+1
	do j = 0, ny+1
	do i = 0, nx+1
		y_sgs(i,j,k) = 0.0d0
		zpc(i,j,k)   = 0.0d0
	enddo
	enddo
	enddo
	!$omp end parallel do

	! y_red = r_red
    !$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		if (iand(i+j+k,1) == 0) then
			y_sgs(i,j,k) = r(i,j,k)
		endif
	enddo
	enddo
	enddo
    !$omp end parallel do
	CALL cc_geral(y_sgs,nx,ny,nz ,1)

	! y_black = r_black + sum(s * y_red(neigh))
    !$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		if (iand(i+j+k,1) == 1) then
			y_sgs(i,j,k) = r(i,j,k) + matsi(i+1,j,k) * y_sgs(i+1,j,k) + matsi(i,j,k) * y_sgs(i-1,j,k) &
								  + matsj(i,j+1,k) * y_sgs(i,j+1,k) + matsj(i,j,k) * y_sgs(i,j-1,k) &
								  + matsk(i,j,k+1) * y_sgs(i,j,k+1) + matsk(i,j,k) * y_sgs(i,j,k-1)
		endif
	enddo
	enddo
	enddo
    !$omp end parallel do
	CALL cc_geral(y_sgs,nx,ny,nz ,1)

	! z_black = y_black
    !$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		if (iand(i+j+k,1) == 1) then
			zpc(i,j,k) = y_sgs(i,j,k)
		endif
	enddo
	enddo
	enddo
    !$omp end parallel do
	CALL cc_geral(zpc,nx,ny,nz ,1)

	! z_red = y_red + sum(s * z_black(neigh))
    !$omp parallel do
	do k = 1, nz
	do j = 1, ny
	do i = 1, nx
		if (iand(i+j+k,1) == 0) then
			zpc(i,j,k) = y_sgs(i,j,k) + matsi(i+1,j,k) * zpc(i+1,j,k) + matsi(i,j,k) * zpc(i-1,j,k) &
								  + matsj(i,j+1,k) * zpc(i,j+1,k) + matsj(i,j,k) * zpc(i,j-1,k) &
								  + matsk(i,j,k+1) * zpc(i,j,k+1) + matsk(i,j,k) * zpc(i,j,k-1)
		endif
	enddo
	enddo
	enddo
    !$omp end parallel do
	CALL cc_geral(zpc,nx,ny,nz ,1)

END SUBROUTINE precond_sgs

!!!!#########!!!!!!!!!##################!!!!!!!!!!#########!!!!
!!!!#########!!!!!!!!!		CG - FIM	!!!!!!!!!!#########!!!!
!!!!#########!!!!!!!!!##################!!!!!!!!!!#########!!!!




!===================================================================================================================
! FFT-based Poisson solver (constant rho, uniform dx/dy/dz)
!===================================================================================================================

SUBROUTINE poisson_fft3d_const_rho(prd1, rhs)

	USE disc, only: nx, ny, nz, dx, dy, dz, dt, dzz, inv_dx2, inv_dy2, twopi
	USE cond, only: ccx0, ccy0
	USE ls_param, only: rho_f1
	use, intrinsic :: iso_c_binding

	IMPLICIT NONE
	include 'fftw3.f03'

	real(8),intent(inout),dimension(0:nx+1,0:ny+1,0:nz+1) :: prd1
	real(8),intent(in),dimension(nx,ny,nz) :: rhs

	integer :: i, j, k
	integer :: nxE, nyE, nzE
	integer :: ii, jj, kk
	integer :: nTot
	real(c_double) :: rho0, dz0
	real(c_double) :: kx, ky, kz
	real(c_double) :: lam, denom
	real(c_double) :: inv_dz2
	logical, save :: printed_once = .false.

	complex(c_double_complex), allocatable :: dataC(:,:,:)
	type(C_PTR) :: plan_f, plan_b

	! Tamanhos estendidos: Neumann -> extensão par para virar periódico; Periódico -> mantém tamanho.
	nxE = nx
	if (ccx0 .ne. 0) nxE = 2*nx
	nyE = ny
	if (ccy0 .ne. 0) nyE = 2*ny
	! z é sempre Neumann no código (cc_geral): estende
	nzE = 2*nz

	! Caso dz não seja uniforme, use FFT em x/y + solução tridiagonal em z
	! (mantém espectral horizontal e permite dz(k) variável).
	dz0 = dzz(1)
	if (abs(dzz(nz)-dz0) > 1.0d-12*max(1.0d0,abs(dz0))) then
		if (.not. printed_once) then
			write(*,*) 'Poisson (t_press=2): FFTW em x/y + tridiagonal em z (dz variável).'
			printed_once = .true.
		endif
		CALL poisson_fftxy_tridiag_z_const_rho(prd1, rhs)
		return
	endif
	if (.not. printed_once) then
		write(*,*) 'Poisson (t_press=2): FFTW 3D (dz uniforme).'
		printed_once = .true.
	endif
	rho0 = rho_f1
	inv_dz2 = 1.0d0/(dz0*dz0)

	allocate(dataC(nxE, nyE, nzE))

	! Monta RHS estendido (even extension para Neumann)
	!$omp parallel do private(ii,jj,kk)
	do k = 1, nzE
		do j = 1, nyE
			do i = 1, nxE
				if (nxE == nx) then
					ii = i
				else
					if (i <= nx) then
						ii = i
					else
						ii = 2*nx - i + 1
					endif
				endif

				if (nyE == ny) then
					jj = j
				else
					if (j <= ny) then
						jj = j
					else
						jj = 2*ny - j + 1
					endif
				endif

				if (k <= nz) then
					kk = k
				else
					kk = 2*nz - k + 1
				endif

				! Consistência com o operador usado no CG: (dt/rho) * Laplaciano(p) = -rhs
				dataC(i,j,k) = cmplx(-rhs(ii,jj,kk), 0.0d0, kind=c_double)
			enddo
		enddo
	enddo
	!$omp end parallel do

	! Planos FFT (C-order): para array Fortran dataC(nxE,nyE,nzE), use (nzE,nyE,nxE)
	plan_f = fftw_plan_dft_3d(nzE, nyE, nxE, dataC, dataC, FFTW_FORWARD, FFTW_ESTIMATE)
	plan_b = fftw_plan_dft_3d(nzE, nyE, nxE, dataC, dataC, FFTW_BACKWARD, FFTW_ESTIMATE)

	call fftw_execute_dft(plan_f, dataC, dataC)

	! Resolve no espaço espectral: (dt/rho0)*lambda(k) * p_hat = rhs_hat
	! lambda(k) = (2cos(kx)-2)/dx^2 + (2cos(ky)-2)/dy^2 + (2cos(kz)-2)/dz^2
	!$omp parallel do private(kx,ky,kz,lam,denom)
	do k = 1, nzE
		kz = twopi * dble(k-1) / dble(nzE)
		do j = 1, nyE
			ky = twopi * dble(j-1) / dble(nyE)
			do i = 1, nxE
				kx = twopi * dble(i-1) / dble(nxE)

				lam = (2.0d0*cos(kx) - 2.0d0)*inv_dx2 + (2.0d0*cos(ky) - 2.0d0)*inv_dy2 + (2.0d0*cos(kz) - 2.0d0)*inv_dz2
				denom = (dt/rho0) * lam

				if (dabs(denom) < 1.0d-30) then
					dataC(i,j,k) = (0.0d0, 0.0d0)
				else
					dataC(i,j,k) = dataC(i,j,k) / denom
				endif
			enddo
		enddo
	enddo
	!$omp end parallel do

	call fftw_execute_dft(plan_b, dataC, dataC)

	! Normalização (FFTW não normaliza a inversa)
	nTot = nxE * nyE * nzE
	dataC = dataC / dble(nTot)

	! Extrai solução no domínio original e aplica CCs (ghosts)
	!$omp parallel do
	do k = 1, nz
		do j = 1, ny
			do i = 1, nx
				prd1(i,j,k) = dble(dataC(i,j,k))
			enddo
		enddo
	enddo
	!$omp end parallel do

	CALL cc_geral(prd1, nx, ny, nz ,1)

	call fftw_destroy_plan(plan_f)
	call fftw_destroy_plan(plan_b)
	call fftw_cleanup()

	deallocate(dataC)

	return

END SUBROUTINE poisson_fft3d_const_rho


!===================================================================================================================
! Poisson solver: FFT em x/y (periodic ou Neumann via extensão par) + FD em z (malha não-uniforme)
! Resolve (dt/rho0) * (Lx + Ly + Lz) p = -rhs, com Lz discretizado em forma conservativa em z.
! z possui Neumann (derivada normal zero) coerente com cc_geral.
!===================================================================================================================

SUBROUTINE poisson_fftxy_tridiag_z_const_rho(prd1, rhs)

	USE disc, only: nx, ny, nz, dx, dy, dz, dt, dzz, inv_dx2, inv_dy2, twopi
	USE cond, only: ccx0, ccy0
	USE ls_param, only: rho_f1
	USE fft_zcoef
	use, intrinsic :: iso_c_binding

	IMPLICIT NONE
	include 'fftw3.f03'

	real(8),intent(inout),dimension(0:nx+1,0:ny+1,0:nz+1) :: prd1
	real(8),intent(in),dimension(nx,ny,nz) :: rhs

	integer :: i, j, k
	integer :: nxE, nyE
	integer :: ii, jj
	integer(c_int) :: rank, howmany
	integer(c_int) :: n(2), inembed(2), onembed(2)
	integer(c_int) :: istride, ostride
	integer(c_int) :: idist, odist
	integer(c_int) :: flags
	integer :: nTotXY
	real(c_double) :: rho0
	real(c_double) :: kx, ky
	real(c_double) :: lamxy
	real(c_double) :: scale
	complex(c_double_complex), allocatable :: dataC(:,:,:)
	complex(c_double_complex), allocatable :: a(:), b(:), c(:), d(:), cp(:), dp(:)
	type(C_PTR) :: plan_f, plan_b
	complex(c_double_complex) :: denom

	! Extensão par (even extension) em x/y quando BC não for periódica.
	nxE = nx
	if (ccx0 .ne. 0) nxE = 2*nx
	nyE = ny
	if (ccy0 .ne. 0) nyE = 2*ny

	rho0 = rho_f1
	nTotXY = nxE * nyE
	scale = 1.0d0 / dble(nTotXY)

	allocate(dataC(nxE, nyE, nz))
	! Coeficientes do operador em z: calculados uma única vez (dependem só de dz/dzz)
	if (.not. allocated(azL)) call init_fft_zcoef()

	! Monta RHS estendido em x/y (mantém z original)
	!$omp parallel do
	do k = 1, nz
		do j = 1, nyE
			do i = 1, nxE
				if (nxE == nx) then
					ii = i
				else
					if (i <= nx) then
						ii = i
					else
						ii = 2*nx - i + 1
					endif
				endif

				if (nyE == ny) then
					jj = j
				else
					if (j <= ny) then
						jj = j
					else
						jj = 2*ny - j + 1
					endif
				endif

				dataC(i,j,k) = cmplx(-rhs(ii,jj,k), 0.0d0, kind=c_double)
			enddo
		enddo
	enddo
	!$omp end parallel do

	! FFT 2D (x/y) para cada plano k (batch = nz)
	rank = 2
	n = (/ int(nyE, c_int), int(nxE, c_int) /)
	howmany = int(nz, c_int)
	inembed = n
	onembed = n
	istride = 1_c_int
	ostride = 1_c_int
	idist = int(nTotXY, c_int)
	odist = int(nTotXY, c_int)
	flags = FFTW_ESTIMATE

	plan_f = fftw_plan_many_dft(rank, n, howmany, dataC, inembed, istride, idist, dataC, onembed, ostride, odist, FFTW_FORWARD, flags)
	plan_b = fftw_plan_many_dft(rank, n, howmany, dataC, inembed, istride, idist, dataC, onembed, ostride, odist, FFTW_BACKWARD, flags)

	call fftw_execute_dft(plan_f, dataC, dataC)

	! Resolve, para cada modo (i,j), um sistema tridiagonal em z:
	! (dt/rho0) * (lamxy*I + Lz) p_hat(:,i,j) = rhs_hat(:,i,j)

	! Importante (OpenMP): a/b/c/d/cp/dp precisam ser privados POR THREAD e alocados dentro
	! da região paralela; caso contrário, "private" em allocatable deixa o status indefinido.
	!$omp parallel private(i,j,k,kx,ky,lamxy,a,b,c,d,cp,dp,denom)
	allocate(a(nz), b(nz), c(nz), d(nz), cp(nz), dp(nz))
	!$omp do
	do j = 1, nyE
		ky = twopi * dble(j-1) / dble(nyE)
		do i = 1, nxE
			kx = twopi * dble(i-1) / dble(nxE)
			lamxy = (2.0d0*cos(kx) - 2.0d0)*inv_dx2 + (2.0d0*cos(ky) - 2.0d0)*inv_dy2

			! Monta tridiagonal (complex por simplicidade; coeficientes são reais)
			if (dabs(lamxy) < 1.0d-30) then
				! (kx,ky)=(0,0): resolve somente em z (Lz). O sistema com Neumann em z tem
				! um modo nulo (constante), então fixamos um gauge (p(1)=0) em vez de zerar
				! toda a coluna (o que eliminaria variações verticais do modo médio horizontal).
				!
				! Além disso, removemos a média (ponderada por dz) do RHS para satisfazer a
				! condição de compatibilidade do Neumann (integral(rhs) = 0) neste modo.
				denom = (0.0d0, 0.0d0)
				do k = 1, nz
					denom = denom + dataC(i,j,k) * cmplx(dz(k), 0.0d0, kind=c_double)
				enddo
				denom = denom / cmplx(sum(dz(1:nz)), 0.0d0, kind=c_double)
				do k = 1, nz
					d(k) = dataC(i,j,k) - denom
				enddo

				! Gauge: p(1)=0
				a(1) = (0.0d0, 0.0d0)
				b(1) = (1.0d0, 0.0d0)
				c(1) = (0.0d0, 0.0d0)
				d(1) = (0.0d0, 0.0d0)

				do k = 2, nz-1
					a(k) = cmplx((dt/rho0) * azL(k), 0.0d0, kind=c_double)
					b(k) = cmplx((dt/rho0) * azD(k), 0.0d0, kind=c_double)
					c(k) = cmplx((dt/rho0) * azU(k), 0.0d0, kind=c_double)
					! d(k) já foi setado acima
				enddo

				a(nz) = cmplx((dt/rho0) * azL(nz), 0.0d0, kind=c_double)
				b(nz) = cmplx((dt/rho0) * azD(nz), 0.0d0, kind=c_double)
				c(nz) = (0.0d0, 0.0d0)
				! d(nz) já foi setado acima
			else
				a(1) = (0.0d0, 0.0d0)
				b(1) = cmplx((dt/rho0) * (azD(1) + lamxy), 0.0d0, kind=c_double)
				c(1) = cmplx((dt/rho0) * azU(1), 0.0d0, kind=c_double)
				d(1) = dataC(i,j,1)

				do k = 2, nz-1
					a(k) = cmplx((dt/rho0) * azL(k), 0.0d0, kind=c_double)
					b(k) = cmplx((dt/rho0) * (azD(k) + lamxy), 0.0d0, kind=c_double)
					c(k) = cmplx((dt/rho0) * azU(k), 0.0d0, kind=c_double)
					d(k) = dataC(i,j,k)
				enddo

				a(nz) = cmplx((dt/rho0) * azL(nz), 0.0d0, kind=c_double)
				b(nz) = cmplx((dt/rho0) * (azD(nz) + lamxy), 0.0d0, kind=c_double)
				c(nz) = (0.0d0, 0.0d0)
				d(nz) = dataC(i,j,nz)
			endif

			! Thomas (eliminação)
			cp(1) = c(1) / b(1)
			dp(1) = d(1) / b(1)
			do k = 2, nz-1
				denom = b(k) - a(k) * cp(k-1)
				cp(k) = c(k) / denom
				dp(k) = (d(k) - a(k) * dp(k-1)) / denom
			enddo
			denom = b(nz) - a(nz) * cp(nz-1)
			dp(nz) = (d(nz) - a(nz) * dp(nz-1)) / denom

			! Back-substitution
			dataC(i,j,nz) = dp(nz)
			do k = nz-1, 1, -1
				dataC(i,j,k) = dp(k) - cp(k) * dataC(i,j,k+1)
			enddo
		enddo
	enddo
	!$omp end do
	deallocate(a, b, c, d, cp, dp)
	!$omp end parallel

	call fftw_execute_dft(plan_b, dataC, dataC)

	! Normalização da inversa (FFTW não normaliza)
	dataC = dataC * scale

	! Extrai solução no domínio original e aplica CCs
	!$omp parallel do
	do k = 1, nz
		do j = 1, ny
			do i = 1, nx
				prd1(i,j,k) = dble(dataC(i,j,k))
			enddo
		enddo
	enddo
	!$omp end parallel do

	CALL cc_geral(prd1, nx, ny, nz ,1)

	call fftw_destroy_plan(plan_f)
	call fftw_destroy_plan(plan_b)
	call fftw_cleanup()

	deallocate(dataC)

	return

END SUBROUTINE poisson_fftxy_tridiag_z_const_rho

!######################################################################################
SUBROUTINE pressh()

	USE velpre
	USE param
!	!$ USE omp_lib

	IMPLICIT NONE

	integer :: i, j, k
	real(8),save :: intz

	!WORK Z-PENCILS (integração em k tem dependência; paraleliza em (i,j))

    !$omp parallel do
	   do k = nz,1 ,-1
	   do j = 1, ny
	   do i = 1, nx
	      if (k == nz) then !top boundary conditions
		 prd1(i,j,k)   = dz(k)*gz*rho(i,j,k)

	      elseif (k == nz-1) then !top boundary conditions

		 prd1(i,j,k)   = gz*(rho(i,j,k)*dz(k)+ rho(i,j,k+1)*dz(k+1))*0.5 + prd1(i,j,k+1)

	      elseif (k < nz-1) then

		 prd1(i,j,k)   = gz * (rho(i,j,k)*dz(k) + 4.*rho(i,j,k+1)*dz(k+1) + rho(i,j,k+2)*dz(k+2)) /3. + prd1(i,j,k+2)

	      endif

	   enddo
	   enddo
	   enddo
    !$omp end parallel do

END SUBROUTINE pressh

