!Subrotina visco --> viscosidade por Smagorinsky-Lilly

!!! ImplementaÃ§Ã£o 19/07/16
! LuÃ­sa Vieira Lucchese

!!! ModificaÃ§Ã£o 18/09/16
! Leonardo Romero Monteiro

!! Modificação 2 31/12/16
! Luísa V. Lucchese

SUBROUTINE visco()

	USE parametros
	USE tempo
        USE smag
        USE velpre

	IMPLICIT NONE
	!===================================================================================================================
	integer :: i, j, k

        real(8) :: p1, p2, p3
                          
        real(8), dimension(nx,ny,nz) :: aux1

        real(8), dimension(nx1,ny,nz) :: dudx_i, dudy_i, dudz_i
        real(8), dimension(nx,ny1,nz) :: dvdx_i, dvdy_i, dvdz_i
        real(8), dimension(nx,ny,nz1) :: dwdx_i, dwdy_i, dwdz_i


        real(8), dimension(nx,ny,nz) :: dudx, dudy, dudz
        real(8), dimension(nx,ny,nz) :: dvdx, dvdy, dvdz
        real(8), dimension(nx,ny,nz) :: dwdx, dwdy, dwdz

        real(8), dimension(nx1,ny,nz) :: dudx_x, dudy_x, dudz_x, dwdx_x, dvdx_x
        real(8), dimension(nx,ny1,nz) :: dudy_y, dvdx_y, dwdy_y, dvdz_y, dvdy_y
        real(8), dimension(nx,ny,nz1) :: dwdy_z, dvdz_z, dwdz_z, dwdx_z, dudz_z
        real(8), dimension(nx1,ny1,nz) :: dvdx_a, dudy_a
        real(8), dimension(nx1,ny,nz1) :: dwdx_a, dudz_a 
        real(8), dimension(nx,ny1,nz1) :: dwdy_a, dvdz_a

	real(8) :: aux2,aux3,aux4,deltag,nuor
        real(8), dimension(nx,ny1,nz) :: ynut_a
        real(8), dimension(nx1,ny,nz) :: xnut_a
        real(8), dimension(nx,ny,nz1) :: znut_a

	real(8), dimension(nx1,ny,nz) :: ls_nux
	real(8), dimension(nx,ny1,nz) :: ls_nuy
	real(8), dimension(nx,ny,nz1) :: ls_nuz



if (m_turb == 0) then ! sem modelo de turbulÃªncia

     xnut_a=0.
     ynut_a=0.
     znut_a=0.

else
    !computando as derivadas

  call derivax(u,nx1,ny,nz,dudx_i); call derivax(v,nx,ny1,nz,dvdx_i); call derivax(w,nx,ny,nz1,dwdx_i)
  call derivay(u,nx1,ny,nz,dudy_i); call derivay(v,nx,ny1,nz,dvdy_i); call derivay(w,nx,ny,nz1,dwdy_i)
  call derivaz(u,nx1,ny,nz,dudz_i); call derivaz(v,nx,ny1,nz,dvdz_i); call derivaz(w,nx,ny,nz1,dwdz_i)


if (m_turb == 1) then ! LES Smagorinsky-Lilly clássico

 !as interpolações abaixo são interpolando tudo para o centro da celula

       ! interpolação para u
       call interpx_fc(dudx_i,nx1,ny,nz,dudx)!dudx(i,j,k)=(dudx_i(i+1,j,k)+dudx_i(i,j,k))*0.5
       call interpx_fc(dudy_i,nx1,ny,nz,dudy)!dudy(i,j,k)=(dudy_i(i+1,j,k)+dudy_i(i,j,k))*0.5
       call interpx_fc(dudz_i,nx1,ny,nz,dudz)!dudz(i,j,k)=(dudz_i(i+1,j,k)+dudz_i(i,j,k))*0.5

       ! interpolação para v
       call interpy_fc(dvdx_i,nx,ny1,nz,dvdx)!dvdx(i,j,k)=(dvdx_i(i,j+1,k)+dvdx_i(i,j,k))*0.5
       call interpy_fc(dvdy_i,nx,ny1,nz,dvdy)!dvdy(i,j,k)=(dvdy_i(i,j+1,k)+dvdy_i(i,j,k))*0.5
       call interpy_fc(dvdz_i,nx,ny1,nz,dvdz)!dvdz(i,j,k)=(dvdz_i(i,j+1,k)+dvdz_i(i,j,k))*0.5

       ! interpolação para w
       call interpz_fc(dwdx_i,nx,ny,nz1,dwdx)!dwdx(i,j,k)=(dwdx_i(i,j,k+1)+dwdx_i(i,j,k))*0.5
       call interpz_fc(dwdy_i,nx,ny,nz1,dwdy)!dwdy(i,j,k)=(dwdy_i(i,j,k+1)+dwdy_i(i,j,k))*0.5
       call interpz_fc(dwdz_i,nx,ny,nz1,dwdz)!dwdz(i,j,k)=(dwdz_i(i,j,k+1)+dwdz_i(i,j,k))*0.5

  do k=1,nz
  do j=1,ny
  do i=1,nx
   p1=dudx(i,j,k)+0.5*(dudy(i,j,k)+dvdx(i,j,k)+dudz(i,j,k)+dwdx(i,j,k))
   p2=dvdy(i,j,k)+0.5*(dvdx(i,j,k)+dudy(i,j,k)+dvdz(i,j,k)+dwdy(i,j,k))
   p3=dwdz(i,j,k)+0.5*(dwdy(i,j,k)+dvdz(i,j,k)+dudz(i,j,k)+dwdx(i,j,k))

   nut(i,j,k)=(dx*dy*dz)**(2./3.)*csmag*csmag*sqrt(2.*(p1+p2+p3)*(p1+p2+p3))*rho(i,j,k)
  enddo
  enddo
  enddo

  call interpx_cf(nut,nx,ny,nz,xnut_a)
  call interpy_cf(nut,nx,ny,nz,ynut_a)
  call interpz_cf(nut,nx,ny,nz,znut_a)


elseif (m_turb == 2) then ! LES Smagorinsky-Lilly direcional ! arrumar para mih!!

    !xnut, ynut e znut interpolados para o mesmo local de cada velocidade u,v e
    !w. 

!xnut
 call interpx_cf(dvdx_i,nx,ny1,nz,dvdx_a) !nx1,ny1,nz
 call interpy_fc(dvdx_a,nx1,ny1,nz,dvdx_x) !nx1,ny,nz

 call interpx_cf(dwdx_i,nx,ny,nz1,dwdx_a)  !nx1,ny,nz1
 call interpz_fc(dwdx_a,nx1,ny,nz1,dwdx_x)  !nx1,ny,nz

!ynut
 call interpy_cf(dudy_i,nx1,ny,nz,dudy_a)
 call interpx_fc(dudy_a,nx1,ny1,nz,dudy_y) !saida nx ny1 nz

 call interpy_cf(dwdy_i,nx,ny,nz1,dwdy_a)
 call interpz_fc(dwdy_a,nx,ny1,nz1,dwdy_y) 

!znut
 call interpz_cf(dudz_i,nx1,ny,nz,dudz_a)
 call interpx_fc(dudz_a,nx1,ny,nz1,dudz_z)

 call interpz_cf(dvdz_i,nx,ny1,nz,dvdz_a)
 call interpy_fc(dvdz_a,nx,ny1,nz1,dvdz_z)

  !aqui, faremos as contabilizacoes para xnut.
  !as p são variaveis auxiliares, que serao sobrescritas 
  !na avaliacao de ynut e znut
  !a notacao _a denomina auxiliar.


  do k=1,nz
  do j=1,ny
  do i=1,nx1
  !contabilizacoes xnut
   p1=dudx_i(i,j,k)+0.5*(dudy_i(i,j,k)+dvdx_x(i,j,k)+dudz_i(i,j,k)+dwdx_x(i,j,k))
   xnut_a(i,j,k)=dx*dx*csmag*csmag*sqrt(2.*p1*p1)
  enddo
  enddo
  enddo

  do k=1,nz
  do j=1,ny1
  do i=1,nx
  !contabilizacoes ynut
   p1=dvdy_i(i,j,k)+0.5*(dvdx_i(i,j,k)+dudy_y(i,j,k)+dvdz_i(i,j,k)+dwdy_y(i,j,k))
   ynut_a(i,j,k)=dy*dy*csmag*csmag*sqrt(2.*p1*p1)
  enddo
  enddo
  enddo

  do k=1,nz1
  do j=1,ny
  do i=1,nx
  !contabilizacoes znut
   p1=dwdz_i(i,j,k)+0.5*(dwdy_i(i,j,k)+dvdz_z(i,j,k)+dudz_z(i,j,k)+dwdx_i(i,j,k))
   znut_a(i,j,k)=dz*dz*csmag*csmag*sqrt(2.*p1*p1)
  enddo
  enddo
  enddo
  !xnut_a(0,:,:)=xnut_a(1,:,:)
  !ynut_a(:,0,:)=ynut_a(:,1,:)
  !znut_a(:,:,0)=znut_a(:,:,1)


 endif

endif

  !valor original de viscosidade
!  xnut_a(:,:,nz-1)  = 0.!znut_a(:,:,nz-1) 
!  xnut_a(:,:,nz)  = 0.

	call interpx_cf(ls_nu,nx,ny,nz,ls_nux) !(nx1,ny,nz)
	call interpy_cf(ls_nu,nx,ny,nz,ls_nuy) !(nx,ny1,nz)
	call interpz_cf(ls_nu,nx,ny,nz,ls_nuz) !(nx,ny,nz1)


  do k=1,nz
  do j=1,ny
  do i=1,nx1
     xnut(i,j,k)=ls_nux(i,j,k)+xnut_a(i,j,k)
  enddo
  enddo
  enddo

!  ynut_a(:,:,nz-1)  = 0.!znut_a(:,:,nz-1) 
!  ynut_a(:,:,nz)  = 0.

  do k=1,nz
  do j=1,ny1
  do i=1,nx
     ynut(i,j,k)=ls_nuy(i,j,k)+ynut_a(i,j,k)
  enddo
  enddo
  enddo

!  znut_a(:,:,nz)  = 0.!znut_a(:,:,nz-1) 
!  znut_a(:,:,nz1) = 0.!znut_a(:,:,nz) 

  do k=1,nz1
  do j=1,ny
  do i=1,nx
     znut(i,j,k)=ls_nuz(i,j,k)+znut_a(i,j,k)
  enddo
  enddo
  enddo







END SUBROUTINE visco
