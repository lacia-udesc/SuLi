!Subrotina de verificação do código utilizando soluções manufaturadas

!!! Implementação 15/05/2017
!!! Modificação   08/08/2020
! Leonardo Romero Monteiro


! obs.:1 ao utilizar esta subrotina, deve-se adiconar os termos fontes para as equações de NS.
! obs.:2 os casos permanentes só são rodados uma vez, no final do código.

! dados para caso 1:




!##! dados para caso 2:
!x = 2.*pi/50., dy = 2.*pi/50., dz1 = pi/25., dztopo = pi/25., dt = 0.001*pi, dt_frame = 0.1*pi
!nx=ceiling(2.*pi/dx) , ny=ceiling(2.*pi/dy), nz=ceiling(pi/dz1) , ts = ceiling(0.1*pi/dt)
!mms_t = 2, ccx0=3, ccxf=3, ccy0=3, ccyf=3, ccz0=3, cczf=3


! cálculo das condições iniciais
SUBROUTINE mms_i()

USE velpre
use mms_m

IMPLICIT NONE
!==================================================================================================================
integer :: i, j, k
real(8) :: c,x,y,h,hpi,dpi
real(8), dimension(nx,ny)  :: h_m
real(8), dimension(nx1,ny) :: h_mx
real(8), dimension(nx,ny1) :: h_my
real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: zx
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: zy
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: zz

hpi=0.5*pi
dpi=2.0*pi

t = 0

! para calcular as alturas para diferentes posições da malha
CALL zs(zx,zy,zz)
		
! cálculo das condições iniciais para verificar com as soluções manufaturadas
if (mms_t == 1) then ! Weng, 2009 - Caso permanente

	do j = 0, ny+1
	y = (j-0.5)*dy
	do i = 0, nx+1
	x = (i-0.5)*dx
	 eta0(i,j) = a*sin(x)*sin(y)
	enddo
	enddo
	
	do j = 1, ny
	y = (j-0.5)*dy
	do i = 1, nx1
	x = (i-1. )*dx
	 h_mx(i,j) = a*sin(x)*sin(y) + h0
	enddo
	enddo

	do j = 1, ny1
	y = (j-1. )*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 h_my(i,j) = a*sin(x)*sin(y) + h0
	enddo
	enddo

	do j = 1, ny
	y = (j-0.5 )*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 h_m(i,j) = eta0(i,j) + h0
	enddo
	enddo
	
	do k = 1, nz
	do j = 1, ny
	y = (j-0.5)*dy
	do i = 1, nx1
	x = (i-1.0)*dx
	 u(i,j,k) =  sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/h_mx(i,j)) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/h_mx(i,j))-1.)
	enddo
	enddo
	enddo


	do k = 1, nz
	do j = 1, ny1
	y = (j-1.0)*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 v(i,j,k) = -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/h_my(i,j)) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/h_my(i,j))-1.)
	enddo
	enddo
	enddo


	do k = 1, nz1
	do j = 1, ny
	y = (j-0.5)*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 w(i,j,k) = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/h_m(i,j)) &
		     -2.*pi*zz(i,j,k)/h_m(i,j)*cos(2.*pi*zz(i,j,k)/h_m(i,j))+2.*pi)/(2.*pi)
	enddo
	enddo
	enddo


elseif (mms_t == 2) then ! Weng, 2009 - Caso não permanente

	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi
	 eta0(i,j) = a*cos(x*0.5)*cos(y*0.5)*cos(t)
	 h_m(i,j) = eta0(i,j) + h0
	enddo
	enddo

	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx1
	x = (i-1. )*dx-pi
	 h_mx(i,j) = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	enddo
	enddo

	do j = 1, ny1
	y = (j-1. )*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi
	 h_my(i,j) = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	enddo
	enddo


	do k = 1, nz
	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx1
	x = (i-1.0)*dx-pi
	 u(i,j,k) = sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/h_mx(i,j)))*sin(t)
	enddo
	enddo
	enddo

	do k = 1, nz
	do j = 1, ny1
	y = (j-1.0)*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi

	v(i,j,k) = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/h_my(i,j)))*sin(t)
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi

	w(i,j,k) = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/h_m(i,j)))*sin(t)
	enddo
	enddo
	enddo
	
endif

open (unit=2031989, action= 'write', file= 'dados//erros_mms.txt', status= 'unknown')
write(2031989,*) "dx", "erro_u1", "dy", "erro_v1", "dz1", "erro_w1"
! cálculo das condições iniciais


END SUBROUTINE mms_i


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!@@@!!!!!!!!!!!!!!!!!!!!!!@@@!!!!!!!!!!!!!!!!!!

! verificação do erro com a solução analítica
SUBROUTINE mms()

USE velpre
USE mms_m

IMPLICIT NONE
!==================================================================================================================
integer :: i, j, k

real(8), dimension(nx1,ny,nz) :: u_m
real(8), dimension(nx,ny1,nz) :: v_m
real(8), dimension(nx,ny,nz1) :: w_m
real(8), dimension(nx,ny,nz)  :: p_m
real(8), dimension(nx,ny)     :: h_m
real(8), dimension(nx1,ny)    :: h_mx
real(8), dimension(nx,ny1)    :: h_my
real(8), dimension(0:nx1+1,0:ny+1,0:nz+1)   :: zx
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1)   :: zy
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1)   :: zz

real(8) :: x,y,hpi,dpi

real(8) :: erro_u1, erro_v1, erro_w1

hpi=0.5*pi
dpi=2.0*pi

erro_t = 0.

! para calcular as alturas para diferentes posições da malha
CALL zs(zx,zy,zz)


! cálculo das soluções manufaturadas
if (mms_t == 1) then ! Weng, 2009 - Caso permanente

	do j = 1, ny
	y = (j-0.5)*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 h_m(i,j) = a*sin(x)*sin(y) + h0
	enddo
	enddo

	do j = 1, ny
	y = (j-0.5)*dy
	do i = 1, nx1
	x = (i-1. )*dx
	 h_mx(i,j) = a*sin(x)*sin(y) + h0
	enddo
	enddo

	do j = 1, ny1
	y = (j-1. )*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 h_my(i,j) = a*sin(x)*sin(y) + h0
	enddo
	enddo

	do k = 1, nz
	do j = 1, ny
	y = (j-0.5)*dy
	do i = 1, nx1
	x = (i-1.0)*dx
	 u_m(i,j,k) =  sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/h_mx(i,j)) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/h_mx(i,j))-1.)
	enddo
	enddo
	enddo

	do k = 1, nz
	do j = 1, ny1
	y = (j-1.0)*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 v_m(i,j,k) = -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/h_my(i,j)) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/h_my(i,j))-1.)
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	y = (j-0.5)*dy
	do i = 1, nx
	x = (i-0.5)*dx
	 w_m(i,j,k) = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/h_m(i,j)) &
		     -2.*pi*zz(i,j,k)/h_m(i,j)*cos(2.*pi*zz(i,j,k)/h_m(i,j))+2.*pi)/(2.*pi)
	enddo
	enddo
	enddo


elseif (mms_t == 2) then ! Weng, 2009 - Caso não permanente

	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi
	 h_m(i,j) = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	enddo
	enddo

	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx1
	x = (i-1. )*dx-pi
	 h_mx(i,j) = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	enddo
	enddo

	do j = 1, ny1
	y = (j-1. )*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi
	 h_my(i,j) = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	enddo
	enddo



	do k = 1, nz
	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx1
	x = (i-1.0)*dx-pi
	 u_m(i,j,k) = sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/h_mx(i,j)))*sin(t)
	enddo
	enddo
	enddo


	do k = 1, nz
	do j = 1, ny1
	y = (j-1.0)*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi
	 v_m(i,j,k) = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/h_my(i,j)))*sin(t)
	enddo
	enddo
	enddo
	

	do k = 1, nz1
	do j = 1, ny
	y = (j-0.5)*dy-pi
	do i = 1, nx
	x = (i-0.5)*dx-pi
	 w_m(i,j,k) = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/h_m(i,j)))*sin(t)
	enddo
	enddo
	enddo



!	do k = 1, nz
!	do j = 1, ny
!	y = (j-0.5)*dy
!	do i = 1, nx
!	x = (i-0.5)*dx
!	 p_m(i,j,k) = c*cos(x)*cos(y)*cos(pi*0.5*zz(i,j,k)/h_m(i,j))
!	enddo
!	enddo
!	enddo

endif


! cálculo dos erros
! a primeira opção é para que seja visualizado o campo exato no paraview
erro_u = 0. 
erro_v = 0.
erro_w = 0.

	do k = 1, nz
	do j = 1, ny
	do i = 1, nx1
		!erro_u(i,j,k)= u_m(i,j,k)
		erro_u(i,j,k)= (u_m(i,j,k)-u(i,j,k))**2
	enddo
	enddo
	enddo

	do k = 1, nz
	do j = 1, ny1
	do i = 1, nx
		!erro_v(i,j,k) = v_m(i,j,k)
		erro_v(i,j,k) = (v_m(i,j,k)-v(i,j,k))**2
	enddo
	enddo
	enddo

	do k = 1, nz1
	do j = 1, ny
	do i = 1, nx
		!erro_w(i,j,k) = w_m(i,j,k)
		erro_w(i,j,k) = (w_m(i,j,k)-w(i,j,k))**2
	enddo
	enddo
	enddo

	! calculando o erro do domínio inteiro
	erro_u1 = sqrt(sum(erro_u)/(nx1*ny*nz) ) 
	erro_v1 = sqrt(sum(erro_v)/(nx*ny1*nz) ) 
	erro_w1 = sqrt(sum(erro_w)/(nx*ny*nz1) ) 


	!!! preparar para plotar (tirar se a primeira opção tiver sido utilizado)
	erro_u = sqrt(abs(erro_u))
	erro_v = sqrt(abs(erro_v))
	erro_w = sqrt(abs(erro_w))


! plotagem dos erros
if(mod(it, ceiling(dt_frame/dt)).eq.0) then
	write(2031989,*) dx, erro_u1, dy, erro_v1, dz1, erro_w1
endif

! erro acumualado para avaliação temporal
erro_t = (erro_u1 + erro_v1 + erro_w1)/3. + erro_t


! vai aparecer no final da simulação no temrinal
if (it == ts) write(*,*) "dt = ", dt, " ",  "erro total acumulado=", erro_t/ts

END SUBROUTINE mms


!!##!!!##!!!!!!!!!!!!###!!##!!!!#!!!#!#!!#!#!#!#!!#!#!#!!#!!!!!!!!!!####!!!!!!#!!#!#!!!!!
!!##!!!##!!!!!!!!!!!!###!!##!!!!#!!!#!#!!#!#!#!#!!#!#!#!!#!!!!!!!!!!####!!!!!!#!!#!#!!!!!

!! definição das condições de contorno prescritas. 
SUBROUTINE mms_bc()

USE dzs
USE velpre
USE parametros
USE cond
use mms_m


IMPLICIT NONE


real(8), dimension(0:nx1+1,0:ny+1,0:nz+1) :: zx
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1) :: zy
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1) :: zz

real(8) :: x,y,ha,hpi,dpi

integer :: i,j,k

!==================================================================================================================


hpi=0.5*pi
dpi=2.0*pi

! para calcular as alturas para diferentes posições da malha
CALL zs(zx,zy,zz)

if ((mms_t == 1) .and. (t == 0.)) then ! Weng, 2009 - Caso permanente

	!! c.c. bx
	do k = 0, nz+1
	do j = 0, ny+1
	y = (j-0.5)*dy

	i = 0
	x = (i-1.0)*dx
	ha = a*sin(x)*sin(y) + h0
	bxx0(j,k)   = sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)

	i = 1
	x = (i-1.0)*dx
	ha = a*sin(x)*sin(y) + h0
	bxx1(j,k)   = sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)

	i = nx1
	x = (i-1.0)*dx
	ha = a*sin(x)*sin(y) + h0
	bxxf(j,k)   = sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)

	i = nx1+1
	x = (i-1.0)*dx
	ha = a*sin(x)*sin(y) + h0
	bxxf1(j,k)  = sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)


	i = 0
	x = (i-.5)*dx
	y = (j-1.)*dy
	ha = a*sin(x)*sin(y) + h0
	bxy0(j,k)  = -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)
	y = (j-.5)*dy
	ha = a*sin(x)*sin(y) + h0
	bxz0(j,k)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		     -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)
	i = nx+1
	x = (i-.5)*dx
	y = (j-1.)*dy
	ha = a*sin(x)*sin(y) + h0
	bxyf(j,k)  = -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)
	y = (j-.5)*dy
	ha = a*sin(x)*sin(y) + h0
	bxzf(j,k)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		     -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)

	enddo
	enddo

	!! c.c. by
	do k = 0, nz+1
	do i = 0, nx+1
	x = (i-0.5)*dx

	j = 0
	y = (j-1.0)*dy
	ha = a*sin(x)*sin(y) + h0
	byy0(i,k)  =  -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)
	j = 1
	y = (j-1.0)*dy
	ha = a*sin(x)*sin(y) + h0
	byy1(i,k)  =  -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)
	j = ny1
	y = (j-1.0)*dy
	ha = a*sin(x)*sin(y) + h0
	byyf(i,k)  =  -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)
	j = ny1+1
	y = (j-1.0)*dy
	ha = a*sin(x)*sin(y) + h0
	byyf1(i,k)  = -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)


	j = 0
	y = (j-.5)*dy
	x = (i-1.)*dx
	ha = a*sin(x)*sin(y) + h0
	byx0(i,k)  = sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)
	x = (i-.5)*dx
	ha = a*sin(x)*sin(y) + h0
	byz0(i,k)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		     -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)
	j = ny+1
	y = (j-.5)*dy
	x = (i-1.)*dx
	ha = a*sin(x)*sin(y) + h0
	byxf(i,k)  =  sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)
	x = (i-.5)*dx
	ha = a*sin(x)*sin(y) + h0
	byzf(i,k)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		     -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)
	enddo
	enddo



	!! c.c. bz
	do j = 0, ny+1
	do i = 0, nx+1
	y = (j-0.5)*dy
	x = (i-0.5)*dx
	ha = a*sin(x)*sin(y) + h0

	k = 0
	bzz0(i,j)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		     -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)
	k = 1
	bzz1(i,j)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		     -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)
	k = nz1
	bzzf(i,j)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		     -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)
	k = nz1+1
	bzzf1(i,j)  = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(2.*pi*zz(i,j,k)/ha) &
		      -2.*pi*zz(i,j,k)/ha*cos(2.*pi*zz(i,j,k)/ha)+2.*pi)/(2.*pi)


	k= 0
	x = (i-1.)*dx
	y = (j-.5)*dy
	ha = a*sin(x)*sin(y) + h0
	bzx0(i,j)  =  sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)
	x = (i-.5)*dx
	y = (j-1.)*dy
	ha = a*sin(x)*sin(y) + h0
	bzy0(i,j)  = -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)

	k = nz+1
	x = (i-1.)*dx
	y = (j-.5)*dy
	ha = a*sin(x)*sin(y) + h0
	bzxf(i,j)  =  sin(x)*cos(y)*sin(pi*0.5*zx(i,j,k)/ha) - cos(x)*sin(y)*(cos(pi*2.*zx(i,j,k)/ha)-1.)
	x = (i-.5)*dx
	y = (j-1.)*dy
	ha = a*sin(x)*sin(y) + h0
	bzyf(i,j)  = -cos(x)*sin(y)*sin(pi*0.5*zy(i,j,k)/ha) + sin(x)*cos(y)*(cos(pi*2.*zy(i,j,k)/ha)-1.)

	enddo
	enddo



elseif (mms_t == 2) then ! Weng, 2009 - Caso não permanente

	!! c.c. bx
	do k = 0, nz+1
	do j = 0, ny+1
	y = (j-0.5)*dy-pi

	i = 0
	x = (i-1.0)*dx-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxx0(j,k)   = sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)

	i = 1
	x = (i-1.0)*dx-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxx1(j,k)   = sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)

	i = nx1
	x = (i-1.0)*dx-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxxf(j,k)   = sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)

	i = nx1+1
	x = (i-1.0)*dx-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxxf1(j,k)  = sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)


	i = 0
	x = (i-.5)*dx-pi
	y = (j-1.)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxy0(j,k)  = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)
	x = (i-.5)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxz0(j,k)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)

	i = nx+1
	x = (i-.5)*dx-pi
	y = (j-1.)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxyf(j,k)  = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)
	x = (i-.5)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bxzf(j,k)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)

	enddo
	enddo


	!! c.c. by
	do k = 0, nz+1
	do i = 0, nx+1
	x = (i-0.5)*dx-pi

	j = 0
	y = (j-1.0)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byy0(i,k)  =  -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)
	j = 1
	y = (j-1.0)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byy1(i,k)  =  -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)
	j = ny1
	y = (j-1.0)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byyf(i,k)  =  -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)
	j = ny1+1
	y = (j-1.0)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byyf1(i,k) =  -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)


	j = 0
	x = (i-1.)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byx0(i,k)  = sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)
	x = (i-.5)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byz0(i,k)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)

	j = ny+1
	x = (i-1.)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byxf(i,k)  =  sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)
	x = (i-.5)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	byzf(i,k)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)

	enddo
	enddo


	!! c.c. bz
	do j = 0, ny+1
	do i = 0, nx+1
	y = (j-0.5)*dy-pi
	x = (i-0.5)*dx-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0

	k = 0
	bzz0(i,j)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)

	k = 1
	bzz1(i,j)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)

	k = nz1
	!bzzf(i,j)  = -a*cos(x*0.5)*cos(y*0.5)*sin(t)
	bzzf(i,j)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)
	
	k = nz1+1
	!bzzf1(i,j)  = -a*cos(x*0.5)*cos(y*0.5)*sin(t)
	bzzf1(i,j)  = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zz(i,j,k)/ha))*sin(t)


	k= 0
	x = (i-1.)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bzx0(i,j)  =  sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)
	x = (i-.5)*dx-pi
	y = (j-1.)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	bzy0(i,j)  = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)

	k = nz+1
	x = (i-1.)*dx-pi
	y = (j-.5)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	!bzxf(i,j)  =  sin(y)*cos(x*0.5)*cos(x*0.5)*sin(t)
	bzxf(i,j)  =  sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zx(i,j,k)/ha))*sin(t)
	x = (i-.5)*dx-pi
	y = (j-1.)*dy-pi
	ha = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	!bzyf(i,j)  = -sin(x)*cos(y*0.5)*cos(y*0.5)*sin(t)
	bzyf(i,j)  = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zy(i,j,k)/ha))*sin(t)
	
	enddo
	enddo
endif


END SUBROUTINE mms_bc


!!##!!!##!!!!!!!!!!!!###!!##!!!!#!!!#!#!!#!#!#!#!!#!#!#!!#!!!!!!!!!!####!!!!!!#!!#!#!!!!!
!!##!!!##!!!!!!!!!!!!###!!##!!!!#!!!#!#!!#!#!#!#!!#!#!#!!#!!!!!!!!!!####!!!!!!#!!#!#!!!!!


! termos-fonte para o caso 1
SUBROUTINE termo_fonte1()

use mms_m
use dzs

IMPLICIT NONE

integer :: i,j,k

real(8) :: x,y,hpi,dpi,h1,h2,h3,zh1,h4,h5,dhdx,dhdy,d2hdx2,d2hdy2,d2hdxy,z1,visc

real(8), dimension(0:nx1+1,0:ny+1,0:nz+1)   :: zx
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1)   :: zy
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1)   :: zz


real(8) :: ux, vx, wx
real(8) :: uy, vy, wy
real(8) :: uz, vz, wz

real(8), dimension(nx,ny)     :: h
real(8), dimension(nx1,ny)    :: hx
real(8), dimension(nx,ny1)    :: hy

real(8) :: dudx,dudy,dudz,d2udx2,d2udy2,d2udz2,d2udxy,d2udxz,d2udyz
real(8) :: dvdx,dvdy,dvdz,d2vdx2,d2vdy2,d2vdz2,d2vdxy,d2vdxz,d2vdyz 
real(8) :: dwdx,dwdy,dwdz,d2wdx2,d2wdy2,d2wdz2,d2wdxy,d2wdxz,d2wdyz

hpi=0.5*pi
dpi=2.0*pi

visc = 0.00100798/998.


CALL zs(zx,zy,zz)

!-------- derivadas para o termo fonte Fu
do k = 1, nz
do j = 1, ny
y = (j-0.5)*dy
do i = 1, nx1
x = (i-1.0)*dx
z1 = zx(i,j,k)

	dhdx=   a*cos(x)*sin(y)
	dhdy=   a*sin(x)*cos(y)
	d2hdx2=-a*sin(x)*sin(y)	
	d2hdy2=-a*sin(x)*sin(y)
	d2hdxy= a*cos(x)*cos(y)

	h1 = a*sin(x)*sin(y) + h0
	h2 = h1*h1
	h3 = h2*h1
	h4 = h3*h1
	h5 = h4*h1

	zh1 = z1/h1

	ux =  sin(x)*cos(y)*sin(hpi*zh1) - cos(x)*sin(y)*(cos(dpi*zh1)-1.)
	vx = -cos(x)*sin(y)*sin(hpi*zh1) + sin(x)*cos(y)*(cos(dpi*zh1)-1.)
	wx = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(dpi*zh1) &
		     -dpi*zh1*cos(dpi*zh1)+dpi)/(dpi)


	dudx= cos(x)*cos(y)*sin(hpi*zh1)+sin(x)*sin(y)*(cos(dpi*zh1)-1.) &
	     -sin(x)*cos(y)*(hpi*z1/h2)*cos(hpi*zh1)*dhdx &
	     -cos(x)*sin(y)*(dpi*z1/h2)*sin(dpi*zh1)*dhdx
        
	dudy=-sin(x)*sin(y)*sin(hpi*zh1)-cos(x)*cos(y)*(cos(dpi*zh1)-1.) &
	     +sin(x)*cos(y)*(-hpi*z1/h2)*cos(hpi*zh1)*dhdy &
	     +cos(x)*sin(y)*(-dpi*z1/h2)*sin(dpi*zh1)*dhdy

	dudz=(hpi/h1)*sin(x)*cos(y)*cos(hpi*zh1)+(dpi/h1)*cos(x)*sin(y)*sin(dpi*zh1)
  
	d2udx2=-ux &
	       	+2.*cos(x)*cos(y)*(-hpi*z1/h2)*cos(hpi*zh1)*dhdx &
		+2.*sin(x)*sin(y)*( dpi*z1/h2)*sin(dpi*zh1)*dhdx &
		+sin(x)*cos(y)*(-(hpi*z1/h2)**2*sin(hpi*zh1)+ (pi*z1/h3)*cos(hpi*zh1))*dhdx**2 &
		+cos(x)*sin(y)*( (dpi*z1/h2)**2*cos(dpi*zh1)+4.*pi*z1/h3*sin(dpi*zh1))*dhdx**2 &
	       -(sin(x)*cos(y)*hpi*z1/h2*cos(hpi*zh1)+cos(x)*sin(y)*(dpi*z1/h2)*sin(dpi*zh1))*d2hdx2


	d2udy2=-ux &
		-2.*sin(x)*sin(y)*(-hpi*z1/h2)*cos(hpi*zh1)*dhdy &
		-2.*cos(x)*cos(y)*( dpi*z1/h2)*sin(dpi*zh1)*dhdy &
		+sin(x)*cos(y)*(-(hpi*z1/h2)**2*sin(hpi*zh1)+ (pi*z1/h3)*cos(hpi*zh1))*dhdy**2 &
		+cos(x)*sin(y)*( (dpi*z1/h2)**2*cos(dpi*zh1)+4.*pi*z1/h3*sin(dpi*zh1))*dhdy**2 &
		-(sin(x)*cos(y)*hpi*z1/h2*cos(hpi*zh1)-cos(x)*sin(y)*(dpi*z1/h2)*sin(dpi*zh1))*d2hdy2

       
	d2udz2=-(hpi/h1)**2.*sin(x)*cos(y)*sin(hpi*zh1) &
		+(dpi/h1)**2.*cos(x)*sin(y)*cos(dpi*zh1)
        

	!d2vdxy=sin(x)*cos(y)*sin(hpi*zh1)-cos(x)*sin(y)*(cos(dpi*zh1)-1.0) &
	!	+sin(x)*sin(y)*(-hpi*z1/h2)*cos(hpi*zh1)*dhdy &
	!	+cos(x)*cos(y)*( dpi*z1/h2)*sin(dpi*zh1)*dhdy &
	!	+cos(x)*cos(y)*( hpi*z1/h2)*cos(hpi*zh1)*dhdx &
	!	-sin(x)*sin(y)*( dpi*z1/h2)*sin(dpi*zh1)*dhdx &
	!	+cos(x)*sin(y)*( (hpi*z1/h2)**2*sin(hpi*zh1)- (pi*z1/h3)*cos(hpi*zh1))*dhdx*dhdy &
	!	+sin(x)*cos(y)*(-(dpi*z1/h2)**2*cos(dpi*zh1)-(4.*pi*z1/h3)*sin(dpi*zh1))*dhdx*dhdy &
	!	+(cos(x)*sin(y)*(hpi*z1/h2)*cos(hpi*zh1)+sin(x)*cos(y)*dpi*z1/h2*sin(dpi*zh1))*d2hdxy

	!d2wdxz=-a/dpi*sin(2.*x)*(dpi*dpi*z1/h2)*sin(dpi*zh1) &
	!	+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2) &
	!	*((8.*pi*pi*z1/h3)*sin(dpi*zh1)+(8.*pi*pi*pi*z1*z1/h4)*cos(dpi*zh1))*dhdx

	tf_u(i,j,k) = ux*dudx + vx*dudy + wx*dudz &
		     -visc*(d2udx2+d2udy2+d2udz2)
		     !-visc*(d2udx2+d2vdxy+d2wdxz)

enddo
enddo
enddo


!-------- derivadas para o termo fonte Fv
do k = 1, nz
do j = 1, ny1
y = (j-1.0)*dy
do i = 1, nx
x = (i-0.5)*dx
z1 = zy(i,j,k)

	dhdx=   a*cos(x)*sin(y)
	dhdy=   a*sin(x)*cos(y)
	d2hdx2=-a*sin(x)*sin(y)	
	d2hdy2=-a*sin(x)*sin(y)
	d2hdxy= a*cos(x)*cos(y)

	h1 = a*sin(x)*sin(y) + h0
	h2 = h1*h1
	h3 = h2*h1
	h4 = h3*h1
	h5 = h4*h1

	zh1 = z1/h1

	uy =  sin(x)*cos(y)*sin(hpi*zh1) - cos(x)*sin(y)*(cos(dpi*zh1)-1.)
	vy = -cos(x)*sin(y)*sin(hpi*zh1) + sin(x)*cos(y)*(cos(dpi*zh1)-1.)
	wy = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(dpi*zh1) &
		     -dpi*zh1*cos(dpi*zh1)+dpi)/(dpi)

	   dvdx= sin(x)*sin(y)*sin(hpi*zh1)+cos(x)*cos(y)*(cos(dpi*zh1)-1.0) &
		+cos(x)*sin(y)*(hpi*z1/h2)*cos(hpi*zh1)*dhdx &
		+sin(x)*cos(y)*(dpi*z1/h2)*sin(dpi*zh1)*dhdx

	   dvdy=-cos(x)*cos(y)*sin(hpi*zh1)-sin(x)*sin(y)*(cos(dpi*zh1)-1.0) &
		+cos(x)*sin(y)*(hpi*z1/h2)*cos(hpi*zh1)*dhdy &
		+sin(x)*cos(y)*(dpi*z1/h2)*sin(dpi*zh1)*dhdy

	   dvdz=-(hpi/h1)*cos(x)*sin(y)*cos(hpi*zh1)-sin(x)*cos(y)*(dpi/h1)*sin(dpi*zh1)

	   d2vdx2=-vy &
		+2.*sin(x)*sin(y)*(-hpi*z1/h2)*cos(hpi*zh1)*dhdx &
		+2.*cos(x)*cos(y)*( dpi*z1/h2)*sin(dpi*zh1)*dhdx &
		-cos(x)*sin(y)*(-(hpi*z1/h2)**2*sin(hpi*zh1)+   pi*z1/h3*cos(hpi*zh1))*dhdx**2 &
		-sin(x)*cos(y)*( (dpi*z1/h2)**2*cos(dpi*zh1)+4.*pi*z1/h3*sin(dpi*zh1))*dhdx**2 &
		+(cos(x)*sin(y)*hpi*z1/h2*cos(hpi*zh1)+sin(x)*cos(y)*dpi*z1/h2*sin(dpi*zh1))*d2hdx2

	   d2vdy2=-vy &
		-2.*cos(x)*cos(y)*(-hpi*z1/h2)*cos(hpi*zh1)*dhdy &
		-2.*sin(x)*sin(y)*( dpi*z1/h2)*sin(dpi*zh1)*dhdy &
		+(cos(x)*sin(y)*hpi*z1/h2*cos(hpi*zh1)+sin(x)*cos(y)*(dpi*z1/h2)*sin(dpi*zh1))*d2hdy2 &  
		+ cos(x)*sin(y)*((hpi*z1/h2)**2*sin(hpi*zh1)-   pi*z1/h3*cos(hpi*zh1))*dhdy*dhdy &
		- sin(x)*cos(y)*((dpi*z1/h2)**2*cos(dpi*zh1)+4.*pi*z1/h3*sin(dpi*zh1))*dhdy*dhdy
		
	   d2vdz2=(hpi/h1)**2*cos(x)*sin(y)*sin(hpi*zh1)-sin(x)*cos(y)*(dpi/h1)**2*cos(dpi*zh1)

	 !  d2udxy=-cos(x)*sin(y)*sin(hpi*zh1)+sin(x)*cos(y)*(cos(dpi*zh1)-1.0) &
	!	+(sin(x)*sin(y)*dhdx-cos(x)*cos(y)*dhdy)*(hpi*z1/h2)*cos(hpi*zh1) &
	!	+(sin(x)*sin(y)*dhdy-cos(x)*cos(y)*dhdx)*(dpi*z1/h2)*sin(dpi*zh1) &
	!	-(sin(x)*cos(y)*(hpi*z1/h2)*cos(hpi*zh1)+cos(x)*sin(y)*(dpi*z1/h2)*sin(dpi*zh1))*d2hdxy &
	!	+ sin(x)*cos(y)*(-(hpi*z1/h2)**2*sin(hpi*zh1)+   pi*z1/h3*cos(hpi*zh1))*dhdx*dhdy &
	!	+ cos(x)*sin(y)*( (dpi*z1/h2)**2*cos(dpi*zh1)+4.*pi*z1/h3*sin(dpi*zh1))*dhdx*dhdy

	  ! d2wdyz= a/dpi*sin(2.*y)*(dpi*dpi*z1/h2)*sin(dpi*zh1) &
	!	+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2) &
	!	*((8.*pi*pi*z1/h3)*sin(dpi*zh1)+(8.*pi*pi*pi*z1*z1/h4)*cos(dpi*zh1))*dhdy	

	tf_v(i,j,k) = uy*dvdx + vy*dvdy + wy*dvdz &
		     -visc*(d2vdx2+d2vdy2+d2vdz2)
		     !-visc*(d2udxy+d2vdy2+d2wdyz) 

enddo
enddo
enddo

!-------- derivatives of w
do k = 1, nz1
y = (j-0.5)*dy
do j = 1, ny
do i = 1, nx
x = (i-0.5)*dx
z1 = zz(i,j,k)

	dhdx=   a*cos(x)*sin(y)
	dhdy=   a*sin(x)*cos(y)
	d2hdx2=-a*sin(x)*sin(y)	
	d2hdy2=-a*sin(x)*sin(y)
	d2hdxy= a*cos(x)*cos(y)

	h1 = a*sin(x)*sin(y) + h0
	h2 = h1*h1
	h3 = h2*h1
	h4 = h3*h1
	h5 = h4*h1

	zh1 = z1/h1

	uz =  sin(x)*cos(y)*sin(pi*0.5*zh1) - cos(x)*sin(y)*(cos(pi*2.*zh1)-1.)
	vz = -cos(x)*sin(y)*sin(pi*0.5*zh1) + sin(x)*cos(y)*(cos(pi*2.*zh1)-1.)
	wz = -a*(sin(x)**2.*cos(y)**2.-cos(x)**2.*sin(y)**2.)*(sin(dpi*zh1) &
	     -dpi*zh1*cos(dpi*zh1)+dpi)/(dpi)

	   dwdx=-a/dpi*sin(2.*x)*(sin(dpi*zh1)-(dpi*zh1)*cos(dpi*zh1)+dpi) &
		+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2)*(dpi*dpi*z1*z1/h3)*sin(dpi*zh1)*dhdx

	   dwdy= a/dpi*sin(2.*y)*(sin(dpi*zh1)-(dpi*zh1)*cos(dpi*zh1)+dpi) &
		+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2)*(dpi*dpi*z1*z1/h3)*sin(dpi*zh1)*dhdy

	   dwdz=-a*dpi*z1/h2*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2)*sin(dpi*zh1)

	   d2wdx2=-a/pi*cos(2.*x)*(sin(dpi*zh1)-(dpi*zh1)*cos(dpi*zh1)+dpi) &
		+a/pi*sin(2.*x)*(dpi*dpi*z1*z1/h3)*sin(dpi*zh1)*dhdx &
		+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2)*(-(dpi**3.*z1**3./h5) &
		*cos(dpi*zh1)-(12.*pi*pi*z1*z1/h4)*sin(dpi*zh1))*dhdx*dhdx &
		+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2)*(4.*pi*pi*z1*z1/h3)*sin(dpi*zh1)*d2hdx2

	   d2wdy2= a/pi*cos(2.*y)*(sin(dpi*zh1)-(dpi*zh1)*cos(dpi*zh1)+dpi) &
		-a/pi*sin(2.*y)*(dpi*dpi*z1*z1/h3)*sin(dpi*zh1)*dhdy &
		+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2)*(-(dpi**3.*z1**3./h5) &
		*cos(dpi*zh1)-(12.*pi*pi*z1*z1/h4)*sin(dpi*zh1))*dhdy*dhdy &
		+a/dpi*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2)*(dpi*dpi*z1*z1/h3)*sin(dpi*zh1)*d2hdy2  
		
	   d2wdz2=-a*(sin(x)**2*cos(y)**2-cos(x)**2*sin(y)**2) &
		*(dpi*dpi*z1/h3*cos(dpi*zh1)+(dpi/h2)*sin(dpi*zh1))
		
	 !  d2vdyz=-cos(x)*cos(y)*( hpi/h1)*cos(hpi*zh1)+sin(x)*sin(y)*(dpi/h1)*sin(dpi*zh1) &
	!	+cos(x)*sin(y)*(-pi*pi*z1/4./h3*sin(hpi*zh1)+hpi/h2*cos(hpi*zh1))*dhdy &
	!	+sin(x)*cos(y)*( 4.*pi*pi*z1/h3*cos(dpi*zh1)+dpi/h2*sin(dpi*zh1))*dhdy

	  ! d2udxz= cos(x)*cos(y)*(hpi/h1)*cos(hpi*zh1)-sin(x)*sin(y)*(dpi/h1)*sin(dpi*zh1) &
	!	+sin(x)*cos(y)*(pi*pi*z1/4./h3*sin(hpi*zh1)-hpi/h2*cos(hpi*zh1))*dhdx &
	!	-cos(x)*sin(y)*(4.*pi*pi*z1/h3*cos(dpi*zh1)+dpi/h2*sin(dpi*zh1))*dhdx

	tf_w(i,j,k) = uz*dwdx + vz*dwdy + wz*dwdz &
		     -visc*(d2wdx2+d2wdy2+d2wdz2)
		     !-visc*(d2udxz+d2vdyz+d2wdz2) 

enddo
enddo
enddo


END SUBROUTINE termo_fonte1


!!##!!!##!!!!!!!!!!!!###!!##!!!!#!!!#!#!!#!#!#!#!!#!#!#!!#!!!!!!!!!!####!!!!!!#!!#!#!!!!!
!!##!!!##!!!!!!!!!!!!###!!##!!!!#!!!#!#!!#!#!#!#!!#!#!#!!#!!!!!!!!!!####!!!!!!#!!#!#!!!!!

! termos-fonte para o caso 2
SUBROUTINE termo_fonte2()

use mms_m
use dzs

IMPLICIT NONE

integer :: i,j,k

real(8) :: c,s,x,y,hpi,dpi,h1,h2,h3,zh1,h4,h5,dhdx,dhdy,d2hdx2,d2hdy2,d2hdxy,z,visc

real(8) :: dcdx,dcdy,dcdz,dsdx,dsdy,dsdz,d2cdx2,d2cdy2,d2cdz2,d2cdxy,d2cdxz,d2cdyz

real(8) :: ux, vx, wx, dudt
real(8) :: uy, vy, wy, dvdt
real(8) :: uz, vz, wz, dwdt

real(8) :: dudx,dudy,dudz,d2udx2,d2udy2,d2udz2,d2udxy,d2udxz,d2udyz
real(8) :: dvdx,dvdy,dvdz,d2vdx2,d2vdy2,d2vdz2,d2vdxy,d2vdxz,d2vdyz 
real(8) :: dwdx,dwdy,dwdz,d2wdx2,d2wdy2,d2wdz2,d2wdxy,d2wdxz,d2wdyz

real(8), dimension(0:nx1+1,0:ny+1,0:nz+1)   :: zx
real(8), dimension(0:nx+1,0:ny1+1,0:nz+1)   :: zy
real(8), dimension(0:nx+1,0:ny+1,0:nz1+1)   :: zz


hpi=0.5*pi
dpi=2.0*pi

visc = 0.00100798/998.

CALL zs(zx,zy,zz)

!-------- derivadas para o termo fonte Fu
do k = 1, nz
do j = 1, ny
y = (j-0.5)*dy-pi
do i = 1, nx1
x = (i-1.0)*dx-pi
z = zx(i,j,k)

	dhdx=   -a*0.5*sin(x*0.5)*cos(y*0.5)*cos(t)
	dhdy=   -a*0.5*cos(x*0.5)*sin(y*0.5)*cos(t)
	d2hdx2= -a*0.25*cos(x*0.5)*cos(y*0.5)*cos(t)
	d2hdy2= -a*0.25*cos(x*0.5)*cos(y*0.5)*cos(t)
	d2hdxy=  a*0.25*sin(x*0.5)*sin(y*0.5)*cos(t)

	h1 = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	h2 = h1*h1
	h3 = h2*h1

	zh1 = z/h1


	c = cos(hpi*(1.0-zh1))
	s = sin(hpi*(1.0-zh1)) 

	dcdx = -hpi*z/h2*s*dhdx
	dcdy = -hpi*z/h2*s*dhdy
	dcdz =  hpi/h1*s

	dsdx =  hpi*z/h2*c*dhdx 	
	dsdy =  hpi*z/h2*c*dhdy	
	dsdz = -hpi/h1*c

	d2cdx2 = -hpi*z/h2*dsdx*dhdx +pi*z/h3*s*dhdx*dhdx -hpi*z/h2*s*d2hdx2
	d2cdy2 = -hpi*z/h2*dsdy*dhdy +pi*z/h3*s*dhdy*dhdy -hpi*z/h2*s*d2hdy2
	d2cdz2 = -(hpi/h1)**2.*c
	d2cdxy = -hpi*z/h2*dsdy*dhdx +pi*z/h3*s*dhdx*dhdy -hpi*z/h2*s*d2hdxy
	d2cdxz =  hpi*z/h2*hpi/h1*c*dhdx-hpi/h2*s*dhdx
	d2cdyz =  hpi*z/h2*hpi/h1*c*dhdy-hpi/h2*s*dhdy

	ux =  sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zh1))*sin(t)
	vx = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zh1))*sin(t)
	wx = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zh1))*sin(t)

	dudt = sin(y)*cos(0.5*x)*cos(0.5*x)*c*cos(t)  

	dudx= -0.5*sin(x)*sin(y)*c*sin(t)+sin(y)*cos(0.5*x)*cos(0.5*x)*dcdx*sin(t)
        
	dudy= cos(0.5*x)*cos(0.5*x)*cos(y)*c*sin(t)+sin(y)*cos(0.5*x)*cos(0.5*x)*dcdy*sin(t)

	dudz= cos(0.5*x)*cos(0.5*x)*sin(y)*dcdz*sin(t)
  
	d2udx2= -0.5*cos(x)*sin(y)*c*sin(t)-sin(x)*sin(y)*dcdx*sin(t)+sin(y)*cos(0.5*x)*cos(0.5*x)*d2cdx2*sin(t)

	d2udy2= -cos(0.5*x)*cos(0.5*x)*sin(y)*c*sin(t) +2.0*cos(y)*cos(0.5*x)*cos(0.5*x)*dcdy*sin(t) &
		+sin(y)*cos(0.5*x)*cos(0.5*x)*d2cdy2*sin(t)

	d2udz2= cos(0.5*x)*cos(0.5*x)*sin(y)*d2cdz2*sin(t)
        
	d2vdxy= 0.5*sin(y)*cos(x)*c*sin(t)+0.5*sin(x)*sin(y)*dcdx*sin(t) &
	       -cos(0.5*y)*cos(0.5*y)*cos(x)*dcdy*sin(t)-sin(x)*cos(0.5*y)*cos(0.5*y)*d2cdxy*sin(t)

	d2wdxz= 0.5*a*sin(0.5*x)*cos(0.5*y)*dcdz*sin(t)-a*cos(0.5*x)*cos(0.5*y)*d2cdxz*sin(t)

	tf_u(i,j,k) = dudt + ux*dudx + vx*dudy + wx*dudz &
		      -visc*(d2udx2+d2udy2+d2udz2)
		     !-visc*(d2udx2+d2udy2+d2udz2+d2udx2+d2vdxy+d2wdxz)


!	dpdx=-sin(x)*cos(y)*cos(hpi*zh1)+cos(x)*cos(y)*sin(hpi*zh1)*hpi*zh1/h1*dhdx
!	tf_px(i,j,k) =  dpdx * coef

enddo
enddo
enddo


!-------- derivadas para o termo fonte Fv
do k = 1, nz
do j = 1, ny1
y = (j-1.0)*dy-pi
do i = 1, nx
x = (i-0.5)*dx-pi
z = zy(i,j,k)

	dhdx=   -a*0.5*sin(x*0.5)*cos(y*0.5)*cos(t)
	dhdy=   -a*0.5*cos(x*0.5)*sin(y*0.5)*cos(t)
	d2hdx2= -a*0.25*cos(x*0.5)*cos(y*0.5)*cos(t)
	d2hdy2= -a*0.25*cos(x*0.5)*cos(y*0.5)*cos(t)
	d2hdxy=  a*0.25*sin(x*0.5)*sin(y*0.5)*cos(t)

	h1 = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	h2 = h1*h1
	h3 = h2*h1

	zh1 = z/h1

	c = cos(hpi*(1.0-zh1))
	s = sin(hpi*(1.0-zh1)) 

	dcdx = -hpi*z/h2*s*dhdx
	dcdy = -hpi*z/h2*s*dhdy
	dcdz =  hpi/h1*s

	dsdx =  hpi*z/h2*c*dhdx 	
	dsdy =  hpi*z/h2*c*dhdy	
	dsdz = -hpi/h1*c

	d2cdx2 = -hpi*z/h2*dsdx*dhdx +pi*z/h3*s*dhdx*dhdx -hpi*z/h2*s*d2hdx2
	d2cdy2 = -hpi*z/h2*dsdy*dhdy +pi*z/h3*s*dhdy*dhdy -hpi*z/h2*s*d2hdy2
	d2cdz2 = -(hpi/h1)**2.*c
	d2cdxy = -hpi*z/h2*dsdy*dhdx +pi*z/h3*s*dhdx*dhdy -hpi*z/h2*s*d2hdxy
	d2cdxz =  hpi*z/h2*hpi/h1*c*dhdx-hpi/h2*s*dhdx
	d2cdyz =  hpi*z/h2*hpi/h1*c*dhdy-hpi/h2*s*dhdy

	uy =  sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zh1))*sin(t)
	vy = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zh1))*sin(t)
	wy = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zh1))*sin(t)

          dvdt=-sin(x)*cos(0.5*y)*cos(0.5*y)*c*cos(t)   
   
	   dvdx= -cos(x)*cos(0.5*y)*cos(0.5*y)*c*sin(t)-sin(x)*cos(0.5*y)*cos(0.5*y)*dcdx*sin(t)

	   dvdy= 0.5*sin(x)*sin(y)*c*sin(t)-sin(x)*cos(0.5*y)*cos(0.5*y)*dcdy*sin(t)

	   dvdz= -cos(0.5*y)*cos(0.5*y)*sin(x)*dcdz*sin(t)

	   d2vdx2= sin(x)*cos(0.5*y)*cos(0.5*y)*c*sin(t) -2.0*cos(x)*cos(0.5*y)*cos(0.5*y)*dcdx*sin(t) &
		  -sin(x)*cos(0.5*y)*cos(0.5*y)*d2cdx2*sin(t)

	   d2vdy2= 0.5*sin(x)*cos(y)*c*sin(t)+sin(x)*sin(y)*dcdy*sin(t)-sin(x)*cos(0.5*y)*cos(0.5*y)*d2cdy2*sin(t)
		
	   d2vdz2= -sin(x)*cos(0.5*y)*cos(0.5*y)*d2cdz2*sin(t)

	   d2udxy= -0.5*sin(x)*cos(y)*c*sin(t)-0.5*sin(x)*sin(y)*dcdy*sin(t) &
		   +cos(0.5*x)*cos(0.5*x)*cos(y)*dcdx*sin(t)+sin(y)*cos(0.5*x)*cos(0.5*x)*d2cdxy*sin(t)

	   d2wdyz= 0.5*a*cos(0.5*x)*sin(0.5*y)*dcdz*sin(t)-a*cos(0.5*x)*cos(0.5*y)*d2cdyz*sin(t)	

	tf_v(i,j,k) = dvdt + uy*dvdx + vy*dvdy + wy*dvdz &
		     -visc*(d2vdx2+d2vdy2+d2vdz2) 
		     !-visc*(d2vdx2+d2vdy2+d2vdz2+d2udxy+d2vdy2+d2wdyz) 

!	dpdy=-cos(x)*sin(y)*cos(hpi*zh1)+cos(x)*cos(y)*sin(hpi*zh1)*hpi*zh1/h1*dhdy
!	tf_py(i,j,k) =  dpdy * coef
enddo
enddo
enddo

!-------- derivatives of w
do k = 1, nz1
do j = 1, ny
y = (j-0.5)*dy-pi
do i = 1, nx
x = (i-0.5)*dx-pi
z = zz(i,j,k)

	dhdx=   -a*0.5*sin(x*0.5)*cos(y*0.5)*cos(t)
	dhdy=   -a*0.5*cos(x*0.5)*sin(y*0.5)*cos(t)
	d2hdx2= -a*0.25*cos(x*0.5)*cos(y*0.5)*cos(t)
	d2hdy2= -a*0.25*cos(x*0.5)*cos(y*0.5)*cos(t)
	d2hdxy=  a*0.25*sin(x*0.5)*sin(y*0.5)*cos(t)

	h1 = a*cos(x*0.5)*cos(y*0.5)*cos(t) + h0
	h2 = h1*h1
	h3 = h2*h1

	zh1 = z/h1

	c = cos(hpi*(1.0-zh1))
	s = sin(hpi*(1.0-zh1)) 

	dcdx = -hpi*z/h2*s*dhdx
	dcdy = -hpi*z/h2*s*dhdy
	dcdz =  hpi/h1*s

	dsdx =  hpi*z/h2*c*dhdx 	
	dsdy =  hpi*z/h2*c*dhdy	
	dsdz = -hpi/h1*c

	d2cdx2 = -hpi*z/h2*dsdx*dhdx +pi*z/h3*s*dhdx*dhdx -hpi*z/h2*s*d2hdx2
	d2cdy2 = -hpi*z/h2*dsdy*dhdy +pi*z/h3*s*dhdy*dhdy -hpi*z/h2*s*d2hdy2
	d2cdz2 = -(hpi/h1)**2.*c
	d2cdxy = -hpi*z/h2*dsdy*dhdx +pi*z/h3*s*dhdx*dhdy -hpi*z/h2*s*d2hdxy
	d2cdxz =  hpi*z/h2*hpi/h1*c*dhdx-hpi/h2*s*dhdx
	d2cdyz =  hpi*z/h2*hpi/h1*c*dhdy-hpi/h2*s*dhdy

	uz =  sin(y)*cos(x*0.5)*cos(x*0.5)*cos(hpi*(1.-zh1))*sin(t)
	vz = -sin(x)*cos(y*0.5)*cos(y*0.5)*cos(hpi*(1.-zh1))*sin(t)
	wz = -a*cos(x*0.5)*cos(y*0.5)*cos(hpi*(1.-zh1))*sin(t)

          dwdt=  -a*cos(0.5*x)*cos(0.5*y)*c*cos(t)

	   dwdx= 0.5*a*sin(0.5*x)*cos(0.5*y)*c*sin(t)-a*cos(0.5*x)*cos(0.5*y)*dcdx*sin(t)

	   dwdy= 0.5*a*cos(0.5*x)*sin(0.5*y)*c*sin(t)-a*cos(0.5*x)*cos(0.5*y)*dcdy*sin(t)

	   dwdz= -a*cos(0.5*x)*cos(0.5*y)*dcdz*sin(t)

	   d2wdx2= 0.25*a*cos(0.5*x)*cos(0.5*y)*c*sin(t)+a*sin(0.5*x)*cos(0.5*y)*dcdx*sin(t) &
		  -a*cos(0.5*x)*cos(0.5*y)*d2cdx2*sin(t)

	   d2wdy2= 0.25*a*cos(0.5*x)*cos(0.5*y)*c*sin(t)+a*cos(0.5*x)*sin(0.5*y)*dcdy*sin(t) &
		  -a*cos(0.5*x)*cos(0.5*y)*d2cdy2*sin(t)
		
	   d2wdz2= -a*cos(0.5*x)*cos(0.5*y)*d2cdz2*sin(t)
		
	   d2vdyz= 0.5*sin(x)*sin(y)*dcdz*sin(t)-sin(x)*cos(0.5*y)*cos(0.5*y)*d2cdyz*sin(t)

	   d2udxz=-0.5*sin(x)*sin(y)*dcdz*sin(t)+sin(y)*cos(0.5*x)*cos(0.5*x)*d2cdxz*sin(t)

	tf_w(i,j,k) = dwdt + uz*dwdx + vz*dwdy + wz*dwdz &
		     -visc*(d2wdx2+d2wdy2+d2wdz2) 
		     !-visc*(d2wdx2+d2wdy2+d2wdz2+d2udxz+d2vdyz+d2wdz2) 

!	dpdz=-cos(x)*cos(y)*sin(hpi*zh1)*hpi/h1
!	tf_pz(i,j,k) =  dpdz * coef

enddo
enddo
enddo


END SUBROUTINE termo_fonte2


