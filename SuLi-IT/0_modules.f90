! características básicas do domínio, escoamento e método numérico SuLi-IT
module disc

! definição do número pi
real(8),parameter :: pi = acos(-1.)   

! discretizações espacial em x, y e z (metros); discretização temporal (segundos)
real(8),parameter :: dx = 2.*pi/50., dy = 2.*pi/50., dz1 = pi/25., dztopo = pi/25., dt = 0.001*pi

! número de tempo por arquivo plotado
real(8),parameter :: dt_frame = 0.1*pi


! número de célular para x, y e z (-)
integer,parameter :: nx=ceiling(2.*pi/dx) , ny=ceiling(2.*pi/dy), nz=ceiling(pi/dz1)!int((0.055-dztopo+dz1)/dz1) !6 para baixo  + 1.5 para cima

!nz=int(10./dz1-0.1+0.5) porque a última célula é maior (0.5)

! tempo de simulação (segundos)
integer,parameter :: ts = ceiling(0.1*pi/dt)


real(8),parameter :: uinicial = 0.

integer,parameter :: der = 1 ! 1 = upwind, 2 = centrado (centrado só para advectivo clássico)
integer,parameter :: adv_type = 1 ! 0 = lagrangiano, 1 = advectivo clássico, 2 = rotacional, 3 = antissimétrico

integer,parameter :: obst_t = 0 ! 0 = sem obst, 1 = dunas, 2 = dunas2, 3 = gaussiano3D, 4 = beji, 5 = delft degrau, 6 = delft 1_2, 7 = SBRH calombos e buracos, 8 = Rigotti

integer,parameter :: m_turb = 0 ! 0 = sem modelo, 1 = LES Smagorinsky-Lilly Clássico, 2 = LES Smagorinsky-Lilly Direcional


integer,parameter :: esp_type = 0 ! 0 = sem camada esponja, 1 = leva em consideração a profundidade, 2 = não leva em consideração a profundidade, 3 = Método da Tangente Hiperbólica

integer,parameter :: wave_t = 0 ! 0 = sem onda, 1 = Stokes I, 2 = Stokes II, 5 = Stokes V

integer,parameter :: mms_t = 2 ! 0 = sem MMS, 1 = MMS permanente, 2 = MMS não permanente


integer,parameter :: chezy_t = 0 ! 0 = sem atrito no fundo, 1 = com atrito no fundo


! número de pontos para x, y e z (-)
integer,parameter :: nx1=nx+1, ny1=ny+1, nz1=nz+1
! para fazer dz variável no espaço inicialmente criar uma função ...
real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: dz = dz1
! índice de tempo e tempo instantâneo
integer :: it
real(8) :: t

end module disc

!###################################################################################
!###################################################################################

!condiçõs de contorno
module cond 

integer,parameter :: ccx0=2 !condicao de contorno parede x=0 --> 0 é periodico, 1 é free-slip, 2 é no-slip, 3 é prescrita, 4 é fluxo validacao, 5 é retorno
integer,parameter :: ccxf=2 !condicao de contorno parede x=xf --> 0 é periodico, 1 é free-slip, 2 é no-slip, 3 é prescrita, 4 é saida livre, 5 é retorno
!só pode usar condição periódica no final quando usar no começo e vice-versa
integer,parameter :: ccy0=2 !condicao de contorno parede y=0 --> 0 é periodico, 1 é free-slip e 2 é no-slip, 3 é prescrita
integer,parameter :: ccyf=2 !condicao de contorno parede y=yf --> 0 é periodico, 1 é free-slip e 2 é no-slip, 3 é prescrita
integer,parameter :: ccz0=2 !condicao de contorno parede z=0 --> 1 é free-slip, 2 é no-slip, 3 é prescrita
integer,parameter :: cczf=3 !condicao de contorno parede z=zf --> 1 é free-slip, 3 é prescrita


end module cond

!###################################################################################
!###################################################################################

!define obstáculos no fundo do canal
module obst

USE disc

! velocidade de fundo (m/s)
real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: ub
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: vb
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: wb

!obstáculo
integer, dimension(0:nx1+1,0:ny+1) :: ku    !indicam até que altura as velocidades tem que ser zeradas (até qual índice k)
integer, dimension(0:nx+1,0:ny1+1) :: kv
integer, dimension(0:nx+1,0:ny+1)  :: kw
integer, parameter :: elev = 5. *dz1      !comprimento a adicionar abaixo

real(8), parameter :: amp = 0.25, comp = 1., fase = -pi/2. !amplitude, comprimento e fase da onda

end module obst

!###################################################################################
!###################################################################################

! define componenstes de velocidade e pressão
module velpre

USE disc

! velocidades para x e z (m/s)
real(8),dimension(0:nx1+1,0:ny+1,0:nz+1) :: u
real(8),dimension(0:nx+1,0:ny1+1,0:nz+1) :: v
real(8),dimension(0:nx+1,0:ny+1,0:nz1+1) :: w

! auxiliares de contorno
real(8),dimension(0:ny+1,0:nz+1)  :: bxx0, bxx1
real(8),dimension(0:ny1+1,0:nz+1) :: bxy0
real(8),dimension(0:ny+1,0:nz1+1) :: bxz0

real(8),dimension(0:ny+1,0:nz+1)  :: bxxf, bxxf1
real(8),dimension(0:ny1+1,0:nz+1) :: bxyf
real(8),dimension(0:ny+1,0:nz1+1) :: bxzf

real(8),dimension(0:nx1+1,0:nz+1)  :: byx0
real(8),dimension(0:nx+1,0:nz+1) :: byy0, byy1
real(8),dimension(0:nx+1,0:nz1+1) :: byz0

real(8),dimension(0:nx1+1,0:nz+1)  :: byxf
real(8),dimension(0:nx+1,0:nz+1) :: byyf, byyf1
real(8),dimension(0:nx+1,0:nz1+1) :: byzf

real(8),dimension(0:nx1+1,0:ny+1)  :: bzx0
real(8),dimension(0:nx+1,0:ny1+1) :: bzy0
real(8),dimension(0:nx+1,0:ny+1) :: bzz0, bzz1

real(8),dimension(0:nx1+1,0:ny+1)  :: bzxf
real(8),dimension(0:nx+1,0:ny1+1) :: bzyf
real(8),dimension(0:nx+1,0:ny+1) :: bzzf, bzzf1

!pressão não-hidrostática (m²/s²)
real(8),dimension(0:nx+1,0:ny+1,0:nz+1) :: prd1, prd0, prd

!pressão hidrostática (m)
real(8),dimension(0:nx+1,0:ny+1) :: eta1, eta0

real(8),dimension(nx,ny) :: deltah
real(8) :: d_max, d_min, b_eta0, b_eta1

end module velpre

!###################################################################################
!###################################################################################

! define parâmetros do escoamento
module parametros

!!! Parâmetros !!!
! implicitness parameter $Patnaik et al. 1987$ (-) , coeficiente de chezy (m**(1/2)/s), declividade na direção x
real(8), parameter :: tetah = 0.6, chez = 99999999999999999999999., decliv = 0.

! componentes gravitacionais
real(8) :: gx = 9.80665 * sin(atan(decliv)) , gz = 9.80665 * cos(atan(decliv)), decliv1

! wind stress coefficient (m/s) e velocidade do vento para x e y (m/s)
!real(8) :: cwind, uwind, vwind

end module parametros


!###################################################################################
!###################################################################################

! defini variações da malha na direção z
module dzs

USE disc

real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: dzx
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: dzy
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: dzz


end module dzs

!###################################################################################
!###################################################################################

! define características para contar o tempo de simulação
module tempo

integer(8) :: cont

! hora e data
integer :: agora(8), agora1(8)
real(8) :: ciclo, prev


end module tempo

!###################################################################################
!###################################################################################

! defini variáveis para contorno com ondas
module wave_c
USE disc

real(8) :: p_w, n_w, c_w, l_w, a_w, f_w, h0_f, l0_w
real(8) :: avel1, avel2, avel3, avel4, avel5
real(8) :: aeta1, aeta2, aeta3, aeta4, aeta5

real(8), dimension(0:nz+1) :: kp


end module wave_c

!###################################################################################
!###################################################################################

!Adicionado por Luísa Lucchese --> 07/16
! variáveis para o LES
module smag

use disc

!valor original de viscosidade cinemática
real(8), parameter :: nuor=0.00000101
  
real(8), parameter :: csmag=0.5 !0.17
real(8), dimension(nx,ny,nz)  :: nut 
real(8), dimension(nx1,ny,nz) :: xnut 
real(8), dimension(nx,ny1,nz) :: ynut 
real(8), dimension(nx,ny,nz1) :: znut 
end module smag

!###################################################################################
!###################################################################################

! Varáveis para a aplicação de soluções manufaturadas
module mms_m

use disc

real(8), parameter :: a  = 0.1
real(8), parameter :: h0 = pi

real(8) :: erro_t ! erro médio acumulado

real(8), dimension(nx,ny,nz)  :: tf_press
real(8), dimension(nx1,ny,nz) :: tf_u, erro_u
real(8), dimension(nx,ny1,nz) :: tf_v, erro_v
real(8), dimension(nx,ny,nz1) :: tf_w, erro_w

end module mms_m


!###################################################################################
!###################################################################################

!Análise final do escoamento
module analis

real(8) :: ut_m, ut_m0

end module analis

!###################################################################################



