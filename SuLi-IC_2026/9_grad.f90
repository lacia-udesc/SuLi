!Subrotina para calcular a pressão dinâmica

!Versão: CG
!Referencia: Casulli (1992 e 1999)
!!! Implementação original 15/04/2014

!Versão: GMRES(3) - Generalized Minimal Residual com 3 direções acumuladas
!Referencia: Saad & Schultz (1986)
!!! Implementação GMRES com múltiplas direções - 28/01/2026

!Leonardo Romero Monteiro

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
		CALL poisson_fft3d_const_rho(prd1, matq)
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

	CALL cc_geral(erro,nx,ny,nz)	
	
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


SUBROUTINE norm(prda, mate, matd, matsi, matsj, matsk,matq,nlock)

	USE disc
	USE cond
		
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

		
		! Condições de contorno, von Neumann
		CALL cc_geral(matd,nx,ny,nz)	
			
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
      	CALL cc_geral(mate,nx,ny,nz)	
	
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

	soma = 0.

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
	if (pcg_prec < 0 .or. pcg_prec > 3) pcg_prec = 2

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
	CALL cc_geral(errop,nx,ny,nz)

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
		CALL cc_geral(errop,nx,ny,nz)

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


	if (prec_t == 0 .or. prec_t == 1) then
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
	elseif (prec_t == 2) then
		CALL precond_zline(r, zpc, matsk)
	elseif (prec_t == 3) then
		CALL precond_sgs(r, zpc, matsi, matsj, matsk)
	else
		write(*,*) 'ERRO: precondicionador PCG inválido.'
		stop
	endif

END SUBROUTINE precond_apply


SUBROUTINE precond_zline(r, zpc, matsk)

	USE disc, only: nx, ny, nz, nx1, ny1, nz1
	IMPLICIT NONE

	real(8),intent(in),dimension(nx,ny,nz1) :: matsk
	real(8),intent(in),dimension(0:nx1,0:ny1,0:nz1) :: r
	real(8),intent(out),dimension(0:nx1,0:ny1,0:nz1) :: zpc

	integer :: i, j, k
	real(8) :: den, a_k, c_k
	real(8) :: cp(nz), dp(nz)

	! inicializa
	zpc = 0.0d0

	! Resolve (por coluna i,j): x_k - a_k x_{k-1} - c_k x_{k+1} = r_k
	! com a_k = matsk(i,j,k) e c_k = matsk(i,j,k+1). Diagonal = 1.
    !$omp parallel do
	do j = 1, ny
	do i = 1, nx
		! k = 1 (sem subdiagonal)
		c_k = -matsk(i,j,2)
		cp(1) = c_k
		dp(1) = r(i,j,1)

		do k = 2, nz-1
			a_k = -matsk(i,j,k)
			c_k = -matsk(i,j,k+1)
			den = 1.0d0 - a_k * cp(k-1)
			cp(k) = c_k / den
			dp(k) = (r(i,j,k) - a_k * dp(k-1)) / den
		enddo

		! k = nz (sem superdiagonal)
		a_k = -matsk(i,j,nz)
		den = 1.0d0 - a_k * cp(nz-1)
		cp(nz) = 0.0d0
		dp(nz) = (r(i,j,nz) - a_k * dp(nz-1)) / den

		! back substitution
		zpc(i,j,nz) = dp(nz)
		do k = nz-1, 1, -1
			zpc(i,j,k) = dp(k) - cp(k) * zpc(i,j,k+1)
		enddo
	enddo
	enddo
    !$omp end parallel do

	CALL cc_geral(zpc,nx,ny,nz)

END SUBROUTINE precond_zline


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
	CALL cc_geral(y_sgs,nx,ny,nz)

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
	CALL cc_geral(y_sgs,nx,ny,nz)

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
	CALL cc_geral(zpc,nx,ny,nz)

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
	CALL cc_geral(zpc,nx,ny,nz)

END SUBROUTINE precond_sgs

!!!!#########!!!!!!!!!##################!!!!!!!!!!#########!!!!
!!!!#########!!!!!!!!!		CG - FIM	!!!!!!!!!!#########!!!!
!!!!#########!!!!!!!!!##################!!!!!!!!!!#########!!!!




!===================================================================================================================
! FFT-based Poisson solver (constant rho, uniform dx/dy/dz)
!===================================================================================================================

SUBROUTINE poisson_fft3d_const_rho(prd1, rhs)

	USE disc, only: nx, ny, nz, dx, dy, dt, dzz
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

	! Para este caso, assume dz uniforme
	dz0 = dzz(1)
	if (abs(dzz(nz)-dz0) > 1.0d-12*max(1.0d0,abs(dz0))) then
		write(*,*) 'ERRO: t_press=2 com FFT requer dz uniforme (dzz(k) constante).'
		stop
	endif
	rho0 = rho_f1

	allocate(dataC(nxE, nyE, nzE))

	! Monta RHS estendido (even extension para Neumann)
	!$omp parallel do
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
	!$omp parallel do
	do k = 1, nzE
		kz = 2.0d0*acos(-1.0d0) * dble(k-1) / dble(nzE)
		do j = 1, nyE
			ky = 2.0d0*acos(-1.0d0) * dble(j-1) / dble(nyE)
			do i = 1, nxE
				kx = 2.0d0*acos(-1.0d0) * dble(i-1) / dble(nxE)

				lam = (2.0d0*cos(kx) - 2.0d0)/(dx*dx) + (2.0d0*cos(ky) - 2.0d0)/(dy*dy) + (2.0d0*cos(kz) - 2.0d0)/(dz0*dz0)
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

	CALL cc_geral(prd1, nx, ny, nz)

	call fftw_destroy_plan(plan_f)
	call fftw_destroy_plan(plan_b)
	call fftw_cleanup()

	deallocate(dataC)

	return

END SUBROUTINE poisson_fft3d_const_rho

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

