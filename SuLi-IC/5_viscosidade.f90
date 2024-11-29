!Subrotina visco --> viscosidade por modelo de turbulência

!!! Implementação 19/07/16
! Luísa Vieira Lucchese

!!! Modificação 
! Leonardo Romero Monteiro em 18/09/16
! Luísa V. Lucchese em 31/12/16
! Bruna Fernanda Soares em 15/08/23

SUBROUTINE visco()

	USE param
	USE tempo
        USE les
        USE velpre

	IMPLICIT NONE
	
	integer :: i, j, k

	real(8),save :: p1, p2, p3, lcar
	
	real(8),dimension(nx,ny,nz) :: Fka, diss, des, e0, fe0, fe1

	real(8),dimension(nx1,ny,nz) :: dudx_i, dudy_i, dudz_i
	real(8),dimension(nx,ny1,nz) :: dvdx_i, dvdy_i, dvdz_i
	real(8),dimension(nx,ny,nz1) :: dwdx_i, dwdy_i, dwdz_i

	real(8),dimension(nx,ny,nz) :: dudx, dudy, dudz
    	real(8),dimension(nx,ny,nz) :: dvdx, dvdy, dvdz
    	real(8),dimension(nx,ny,nz) :: dwdx, dwdy, dwdz

	real(8),dimension(nx1,ny,nz) :: dudx_x, dudy_x, dudz_x, dwdx_x, dvdx_x
	real(8),dimension(nx,ny1,nz) :: dudy_y, dvdx_y, dwdy_y, dvdz_y, dvdy_y
	real(8),dimension(nx,ny,nz1) :: dwdy_z, dvdz_z, dwdz_z, dwdx_z, dudz_z
	real(8),dimension(nx1,ny1,nz) :: dvdx_a, dudy_a
	real(8),dimension(nx1,ny,nz1) :: dwdx_a, dudz_a 
   	real(8),dimension(nx,ny1,nz1) :: dwdy_a, dvdz_a

	real(8),save :: aux2,aux3,aux4,deltag,nuor

	!===================================================================================================================
	if (m_turb == 0) then ! 0 = sem modelo de turbulência

		mut = 0.
		ka  = 0.

	else
    
		!computando as derivadas
		call derivax(u,nx1,ny,nz,dudx_i); call derivax(v,nx,ny1,nz,dvdx_i); call derivax(w,nx,ny,nz1,dwdx_i)
		call derivay(u,nx1,ny,nz,dudy_i); call derivay(v,nx,ny1,nz,dvdy_i); call derivay(w,nx,ny,nz1,dwdy_i)
  		call derivaz(u,nx1,ny,nz,dudz_i); call derivaz(v,nx,ny1,nz,dvdz_i); call derivaz(w,nx,ny,nz1,dwdz_i)
  		
 		! as interpolações abaixo são interpolando tudo para o centro da celula

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

		!===================================================================================================================
		if (m_turb == 1) then ! 1 = LES Smagorinsky-Lilly clássico
		
			ka  = 0.
			
			do k=1,nz
			do j=1,ny
			do i=1,nx
	   			p1=dudx(i,j,k)+0.5*(dudy(i,j,k)+dvdx(i,j,k)+dudz(i,j,k)+dwdx(i,j,k)) ! p1+p2+p3 = Sij
	   			p2=dvdy(i,j,k)+0.5*(dvdx(i,j,k)+dudy(i,j,k)+dvdz(i,j,k)+dwdy(i,j,k))
	   			p3=dwdz(i,j,k)+0.5*(dwdy(i,j,k)+dvdz(i,j,k)+dudz(i,j,k)+dwdx(i,j,k))

	   			mut(i,j,k)=(dx*dy*dz)**(2./3.)*csmag*csmag*sqrt(2.*(p1+p2+p3)*(p1+p2+p3))*rho(i,j,k) ! viscosidade dinâmica modelada
			enddo
			enddo
			enddo

		!===================================================================================================================
		elseif (m_turb .ge. 2) then ! 2 = LES - energia cinemática turbulenta, 3 = DES (Heinz, 2020)
			
		! comprimento característico	
		lcar = max(prof,ampl)
			
			do tt = 1, ntt
		
				dt   = a_dt(tt)
						
				nut  = cmu*ka(1:nx,1:ny,1:nz)**(0.5)*delta

				diss = ka(1:nx,1:ny,1:nz)**(1.5)/lcar ! taxa de dissipação > denominador = comprimento caracteristico
					! sem obstáculo > comprimento caracteristico = prof
					! com obstáculo > comprimento caracteristico = 0.2m altura do degrau (backward facing step/T1 - Delft, 1980)

				if (m_turb == 3) then
					des = ka(1:nx,1:ny,1:nz)**(1.5)/(cdes*delta*diss)  
					des = max(1.,des) ! modificação da equação k para DES
				else
					des = 1.
				endif
						
				CALL convdiff_les(Fka)       ! termos convectivos e difusivos da energia cinética
				
				do k = 1, nz                 ! adição dos outros termos da energia cinética
				do j = 1, ny
				do i = 1, nx
					p1=dudx(i,j,k)+0.5*(dudy(i,j,k)+dvdx(i,j,k)+dudz(i,j,k)+dwdx(i,j,k)) ! p1+p2+p3 = Sij
	   				p2=dvdy(i,j,k)+0.5*(dvdx(i,j,k)+dudy(i,j,k)+dvdz(i,j,k)+dwdy(i,j,k))
	   				p3=dwdz(i,j,k)+0.5*(dwdy(i,j,k)+dvdz(i,j,k)+dudz(i,j,k)+dwdx(i,j,k))
 				  												
					Fka(i,j,k) = Fka(i,j,k) + nut(i,j,k)*2.*(p1+p2+p3)*(p1+p2+p3) - des(i,j,k)*diss(i,j,k)
				enddo
				enddo
				enddo
				
				CALL tempo_esc(Fka,nx,ny,nz,ka(1:nx,1:ny,1:nz),e0,fe0,fe1) ! evolução temporal da energia cinética turbulenta (ka)

				CALl contorno_les

				do k = 1, nz
				do j = 1, ny
				do i = 1, nx				
				 
					if (ka(i,j,k) < 0) ka(i,j,k) = 0.
					
				enddo
				enddo
				enddo

			enddo
			
			nut = cmu*ka(1:nx,1:ny,1:nz)**(0.5)*delta ! viscosidade modelada
			mut = nut*rho     ! viscosidade dinâmica modelada							
		endif

	endif
	
	!===================================================================================================================
	
	mu = ls_mu + mut ! viscosidade dinâmica total
	
	call interpx_cf(mu,nx,ny,nz,xmu)
	call interpy_cf(mu,nx,ny,nz,ymu)
	call interpz_cf(mu,nx,ny,nz,zmu)

END SUBROUTINE visco
