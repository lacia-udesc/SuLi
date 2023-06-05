!Características básicas do domínio, escoamento e método numérico SuLi-IC
module disc
	
	integer :: numb_threads
	real(8),parameter ::  pi = acos(-1.) 

	!Discretizações espaciais em x e y (metros), discretização temporal (segundos)
	real(8) :: dx, dy, dz, lx, ly, lz
	integer :: nx, ny, nz
	integer :: nx1, ny1, nz1
	integer :: ts
	
	real(8) :: t, dt, t_i, t_a, t_s, dt0 
	!Número de células para x, y e z (-); número de pontos para x, y e z (-); tempo de simulação (segundos)

	!Número de tempo por arquivo plotado
	real(8) :: dt_frame 

	!Para fazer dz variável no espaço inicialmente criar uma função ...
	real(8) :: uinicial

	integer :: t_plot ! 0 = modo simples (velocidade, Level Set e IBM), 1 = modo completo (pressão, vorticidade, viscosidade)

	integer :: t_tempo ! 0 = Euler Explícito, 1 = RK 2, 2 = RK 3, 3 = AB2, 4 = AB3
	integer :: t_tempo_var  ! 0 = dt constante, 1 = dt adaptativo

	integer :: der  ! 1 = upwind, 2 = centrado, 3 = upwind 2nd order (centrado só para advectivo clássico)
	integer :: adv_type  ! 1 = advectivo clássico, 2 = rotacional, 3 = antissimétrico
		
	integer :: obst_t  ! 0 = sem obst, 1 = dunas, 2 = dunas2, 3 = gaussiano3D, 4 = beji, 5 = delft degrau, 6 = delft 1_2, 7 = SBRH calombos e buracos, 8 = fennema1990, 9 = aureli2008, bd_koshizuka1995eKleefsman2005

	integer :: m_turb ! 0 = sem modelo, 1 = LES Smagorinsky-Lilly Clássico, 2 = LES Smagorinsky-Lilly Direcional

	integer :: esp_type  ! 0 = sem camada esponja, 1 = leva em consideração a profundidade, 2 = não leva em consideração a profundidade, 3 = Método da Tangente Hiperbólica

	integer :: wave_t ! 0 = sem onda, 1 = Stokes I, 2 = Stokes II, 5 = Stokes V

	integer :: t_press  ! 0 = aproximacao hidrostatica de pressao, 1 = aproximacao nao-hidrostatica (1 mais indicado)

	integer :: mms_t   ! 0 = sem MMS, 1 = MMS permanente, 2 = MMS não permanente

	integer :: it, tt ,ntt

	real(8),dimension(3) :: a_dt
	real(8) :: ampl, lambdax, lambday, prof, m

end module disc

module restart
	 real(8) :: irest !0=essa simulacao nao é um restart de outra 1 = o arquivo de restart vai ser lido e usado pra continuar a simulacao
	 real(8) :: interv_rest !de quantas em quantas iteracoes salva o restart
end module restart

!Condicoes de contorno
module cond 

	integer :: ccx0 !condicao de contorno parede x=0 --> 0 é periodico, 1 é free-slip, 2 é no-slip, 3 é prescrita, 4 é fluxo validacao
	integer :: ccxf !condicao de contorno parede x=xf --> 0 é periodico, 1 é free-slip, 2 é no-slip, 3 é prescrita, 4 é saida livre

	!Só pode usar condição periódica no final quando usar no começo e vice-versa

	integer :: ccy0 !condicao de contorno parede y=0 --> 0 é periodico, 1 é free-slip e 2 é no-slip, 3 é prescrita
	integer :: ccyf !condicao de contorno parede y=yf --> 0 é periodico, 1 é free-slip e 2 é no-slip, 3 é prescrita
	integer :: ccz0 !condicao de contorno parede z=0 --> 1 é free-slip, 2 é no-slip, 3 é prescrita
	integer :: cczf !condicao de contorno parede z=zf --> 1 é free-slip, 3 é prescrita


end module cond

module param

	!Parâmetros
	!Viscosidade cinemática (m²/s), coeficiente de chezy (m(1/2)/s), aceleração da gravidade (m/s²) e implicitness parameter $Patnaik et al. 1987$ (-) 
	real(8) :: chezy, decliv
	real(8) :: gx, gz
	
	!real(8), parameter :: gx = 0. , gz = 0.

	!wind stress coefficient (m/s) e velocidade do vento para x e y (m/s)

	!real(8) :: cwind, uwind, vwind

end module param

module tempo

	integer :: cont

	! hora e data
	integer :: agora(8), agora1(8)
	real(8) :: ciclo, prev


end module tempo

module obst

	!Velocidade de fundo (m/s)
	!Velocidade de fundo (m/s)
	real(8), allocatable, dimension(:,:,:)  :: ub, vb, wb

	!Obstáculo: indicam até que altura as velocidades tem que ser zeradas (até qual índice k)
	integer, allocatable, dimension(:,:)  :: ku, kv, kw

end module obst

subroutine init_variables1

USE disc !allocatable
USE obst 

implicit none

	allocate(ub(0:nx1+1,0:ny+1,0:nz+1))
	allocate(vb(0:nx+1,0:ny1+1,0:nz+1))
	allocate(wb(0:nx+1,0:ny+1,0:nz1+1))
	allocate(ku(0:nx1+1,0:ny+1))
	allocate(kv(0:nx+1,0:ny1+1))
	allocate(kw(0:nx+1,0:ny+1))
	
end subroutine init_variables1


module velpre

	USE disc

	!Velocidades para x e z (m/s)
	real(8), allocatable, dimension(:,:,:) :: u,v,w

	real(8), allocatable, dimension(:,:) :: bxx0, bxx1, blx1, bxy0, bxz0, bxxf, bxxf1, bxyf, bxzf

	real(8), allocatable, dimension(:,:) :: byx0, byy0, byy1, byz0, byxf, byyf, byyf1, byzf

	real(8), allocatable, dimension(:,:) :: bzx0, bzy0, bzz0, bzz1, bzxf, bzyf, bzzf, bzzf1

	!Pressão não-hidrostática (m²/s²)
	real(8), allocatable, dimension(:,:,:) :: prd1, prd0, prd, rho, ls_nu

	real(8) :: d_max, d_min, b_eta0, b_eta1

end module velpre

subroutine init_variables2

	USE velpre
	
	implicit none
	
	allocate(u(0:nx1+1,0:ny+1,0:nz+1))
	allocate(v(0:nx+1,0:ny1+1,0:nz+1))
	allocate(w(0:nx+1,0:ny+1,0:nz1+1))

	allocate(bxx0(0:ny+1,0:nz+1), bxx1(0:ny+1,0:nz+1), blx1(0:ny+1,0:nz+1), bxy0(0:ny1+1,0:nz+1))
	allocate(bxz0(0:ny+1,0:nz1+1), bxxf(0:ny+1,0:nz+1), bxxf1(0:ny+1,0:nz+1))
	allocate(bxyf(0:ny1+1,0:nz+1), bxzf(0:ny+1,0:nz1+1), byx0(0:nx1+1,0:nz+1))
	
	allocate(byy0(0:nx+1,0:nz+1), byy1(0:nx+1,0:nz+1)) 
	allocate(byz0(0:nx+1,0:nz1+1))
	allocate(byyf(0:nx+1,0:nz+1), byyf1(0:nx+1,0:nz+1), byzf(0:nx+1,0:nz1+1))
	
	allocate(bzx0(0:nx1+1,0:ny+1), bzy0(0:nx+1,0:ny1+1))
	allocate(bzz0(0:nx+1,0:ny+1), bzz1(0:nx+1,0:ny+1))
	allocate(bzxf(0:nx1+1,0:ny+1), bzyf(0:nx+1,0:ny1+1), bzzf(0:nx+1,0:ny+1), bzzf1(0:nx+1,0:ny+1))
	
	allocate(prd1(0:nx+1,0:ny+1,0:nz+1), prd0(0:nx+1,0:ny+1,0:nz+1), prd(0:nx+1,0:ny+1,0:nz+1))
	allocate(rho(nx,ny,nz), ls_nu(nx,ny,nz))

end subroutine init_variables2



module vartempo

	USE disc
	
	real(8), dimension(:,:,:), allocatable :: Fu, u0, fu0, fu1
	real(8), dimension(:,:,:), allocatable :: Fv, v0, fv0, fv1
	real(8), dimension(:,:,:), allocatable :: Fw, w0, fw0, fw1

end module vartempo

subroutine init_variables3

	USE vartempo
	
	implicit none
	
	allocate(Fu(nx1,ny,nz))
	allocate(u0(nx1,ny,nz))
	allocate(fu0(nx1,ny,nz))
	allocate(fu1(nx1,ny,nz))

	allocate(Fv(nx,ny1,nz))
	allocate(v0(nx,ny1,nz))
	allocate(fv0(nx,ny1,nz))
	allocate(fv1(nx,ny1,nz))

	allocate(Fw(nx,ny,nz1))
	allocate(w0(nx,ny,nz1))
	allocate(fw0(nx,ny,nz1))
	allocate(fw1(nx,ny,nz1))

end subroutine init_variables3



module wave_c
	USE disc

	real(8) :: p_w, n_w, c_w, l_w, a_w, f_w, h0_f, l0_w
	real(8) :: avel1, avel2, avel3, avel4, avel5
	real(8) :: aeta1, aeta2, aeta3, aeta4, aeta5
	real(8), allocatable :: kp(:)

end module wave_c

subroutine init_variables4

	USE wave_c
	
	implicit none
	
	allocate(kp(0:nz+1))	

end subroutine init_variables4




module smag

	USE disc

	real(8) :: csmag
	real(8), allocatable :: nut(:,:,:)
	real(8), allocatable :: xnut(:,:,:)
	real(8), allocatable :: ynut(:,:,:)
	real(8), allocatable :: znut(:,:,:) 

end module smag

subroutine init_variables5

	USE smag
	
	implicit none
	
	allocate(nut(nx,ny,nz))
	allocate(xnut(nx1,ny,nz))
	allocate(ynut(nx,ny1,nz))
	allocate(znut(nx,ny,nz1))

end subroutine init_variables5




!Variáveis para o LES
module ls_param

	USE disc

	real(8), allocatable :: ls(:,:,:), mod_ls(:,:,:), kurv(:,:,:), hs(:,:,:), ddlsdx(:,:,:), ddlsdy(:,:,:), ddlsdz(:,:,:)
	real(8), allocatable :: hsx(:,:,:), hsy(:,:,:), hsz(:,:,:)
	real(8) :: dtau, alpha1, mi_f1, mi_f2, rho_f1, rho_f2 , vol_ini, vol_ins, ls_m, rho_m, sigma
	real(8), dimension(3) :: adtl, bdtl, gdtl
	real(8) :: dx1
	real(8) :: dt1
	integer, parameter :: t_hs = 0
	integer :: tipo  !1 = Onda; 2 = Onda 2; 3 = Barragem; 4 = Gota; 5 = fennema1990 ou aureli2008; 6 = MMS; 7 = Koshizuka et al., 1995
	
end module ls_param

subroutine init_variables6

	USE ls_param
	
	implicit none
	
	allocate(ls(nx,ny,nz), mod_ls(nx,ny,nz), kurv(nx,ny,nz), hs(nx,ny,nz), ddlsdx(nx,ny,nz), ddlsdy(nx,ny,nz), ddlsdz(nx,ny,nz))
	allocate(hsx(nx,ny,nz), hsy(nx,ny,nz), hsz(nx,ny,nz))

end subroutine init_variables6



module mms_m

	USE disc

	real(8), parameter :: a  = 0.!5
	real(8), parameter :: h0 = 2.

	! dt / rho?
	real(8), parameter :: coef = 0. ! 1.
	real(8) :: erro_t ! erro médio acumulado
	real(8), dimension(:,:,:), allocatable :: tf_p, erro_p
	real(8), dimension(:,:,:), allocatable :: tf_u, erro_u, tf_px
	real(8), dimension(:,:,:), allocatable :: tf_v, erro_v, tf_py
	real(8), dimension(:,:,:), allocatable :: tf_w, erro_w, tf_pz

end module mms_m

subroutine init_variables7

	USE mms_m
	
	implicit none
	
	! Allocate the arrays with dimensions nx, ny, nz
	allocate(tf_p(nx, ny, nz), erro_p(nx, ny, nz))
	allocate(tf_u(nx1, ny, nz), erro_u(nx1, ny, nz), tf_px(nx1, ny, nz))
	allocate(tf_v(nx, ny1, nz), erro_v(nx, ny1, nz), tf_py(nx, ny1, nz))
	allocate(tf_w(nx, ny, nz1), erro_w(nx, ny, nz1), tf_pz(nx, ny, nz1))

	! Use the arrays here

	! Deallocate the arrays when they are no longer needed
	!deallocate(tf_p, erro_p)
	!deallocate(tf_u, erro_u, tf_px)
	!deallocate(tf_v, erro_v, tf_py)
	!deallocate(tf_w, erro_w, tf_pz)

end subroutine init_variables7


