  subroutine visu ()

  USE disc

  implicit none
	
  real(8) :: lx, ly, lz
  integer :: nfiles, icrfile, file1, filen, ifile, dig1, dig2, dig3, dig4, dig5
  real(8), allocatable :: y1(:), y2(:), y3(:)
  integer :: i, j, k, num, aig, ii, nfil

  character(100) :: a
  character(5) :: chits

!  write (*,*) 'file'
  a = "campos_"! nome do arquivo de leitura 



  file1 = 10 !valor do primeiro arquivo de leitura ex: campo_000010
  filen = ts / ceiling(dt_frame/real(dt)) + file1 !valor do último arquivo de leitura ex: campo_000800
  icrfile = 1 !incremento entre arquivos de leitura
  nfiles = int((filen-file1)/icrfile) ! número total de arquivos lidos

 lx = real(dx) * real ( nx, 8 )
 ly = real(dy) * real ( ny, 8 )
 lz = real(dz) * real ( nz, 8 )

! chits = "0010"

  nfil=41
  open(nfil,file='arquivos//visu.xdmf')

  write(nfil,'(A22)')'<?xml version="1.0" ?>'
  write(nfil,*)'<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>'
  write(nfil,*)'<Xdmf xmlns:xi="http://www.w3.org/2001/XInclude" Version="2.0">'
  write(nfil,*)'<Domain>'

  write(nfil,'(/)')
  write(nfil,*)'    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
  write(nfil,*)'        <Time TimeType="HyperSlab">'
  write(nfil,*)'            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
  write(nfil,*)'           <!--Start, Stride, Count-->'
  write(nfil,*)'            0.0', dt_frame
  write(nfil,*)'            </DataItem>'
  write(nfil,*)'        </Time>'

  do ifile = file1, filen, icrfile

   dig1 =    ifile/10000 + 48
   dig2 = ( ifile - 10000*( ifile/10000 ) )/1000 + 48
   dig3 = ( ifile - 1000*( ifile/1000 ) )/100 + 48
   dig4 = ( ifile - 100*( ifile/100 ) )/10 + 48
   dig5 = ( ifile - 10*( ifile/10 ) )/1 + 48
   chits(1:5) = char(dig1)//char(dig2)//char(dig3)//char(dig4)//char(dig5)


   ! write(*,*) ifile, trim(a)//chits

    write(nfil,'(/)')
    write(nfil,*)'        <Grid Name="'//trim(a)//chits//'" GridType="Uniform">'
  write(nfil,*)'    <Topology name="topo" TopologyType="3DSMesh"'
  write(nfil,*)'        Dimensions="',nz1,ny1,nx1,'">'
  write(nfil,*)'    </Topology>'

  write(nfil,*)'    <Geometry name="geo" Type="X_Y_Z">'
  write(nfil,*)'    <DataItem Name="X" Format="Binary" Dimensions="',nz1,ny1,nx1,'" DataType="Float" Precision="4" Endian="little"'
  write(nfil,*)'    Seek="',nx1*ny1*nz1*4*0+4,'">'
  write(nfil,*)'    '//trim(a)//chits
  write(nfil,*)'    </DataItem>'

  write(nfil,*)'    <DataItem Name="Y" Format="Binary" Dimensions="',nz1,ny1,nx1,'" DataType="Float" Precision="4" Endian="little"'
  write(nfil,*)'    Seek="',nx1*ny1*nz1*4*1+4,'">'
  write(nfil,*)'    '//trim(a)//chits
  write(nfil,*)'    </DataItem>'

  write(nfil,*)'    <DataItem Name="Z" Format="Binary" Dimensions="',nz1,ny1,nx1,'" DataType="Float" Precision="4" Endian="little"'
  write(nfil,*)'    Seek="',nx1*ny1*nz1*4*2+4,'">'
  write(nfil,*)'    '//trim(a)//chits
  write(nfil,*)'    </DataItem>'
  write(nfil,*)'    </Geometry>'


    write(nfil,*)'            <Attribute Name="V" AttributeType="Vector" Center="Node">'
    write(nfil,*)'               <DataItem Dimensions="',nz1,ny1,nx1,'3" Function="JOIN($0, $1, $2)"'
    write(nfil,*)'                         ItemType="Function">'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*3+4,'" Dimensions="',nz1,ny1,nx1,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*4+4,'" Dimensions="',nz1,ny1,nx1,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*5+4,'" Dimensions="',nz1,ny1,nx1,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'


    write(nfil,*)'            <Attribute Name="level_set" Center="Cell">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="4" Endian="little" Seek="',nx1*ny1*nz1*4*6+4,'"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  '//trim(a)//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'


    write(nfil,*)'            <Attribute Name="obstaculo" Center="Cell">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="4" Endian="little" Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*1+4,'"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  '//trim(a)//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'


    if (t_plot == 1) then

    write(nfil,*)'            <Attribute Name="pressure" Center="Cell">'
    write(nfil,*)'               <DataItem Format="Binary" '
    write(nfil,*)'                DataType="Float" Precision="4" Endian="little" Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*2+4,'"'
    write(nfil,*)'                Dimensions="',nz,ny,nx,'">'
    write(nfil,*)'                  '//trim(a)//chits
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'


    write(nfil,*)'            <Attribute Name="vorticity" AttributeType="Vector" Center="Cell">'
    write(nfil,*)'               <DataItem Dimensions="',nz,ny,nx,'3" Function="JOIN($0, $1, $2)"'
    write(nfil,*)'                         ItemType="Function">'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*3+4,'" Dimensions="',nz,ny,nx,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*4+4,'" Dimensions="',nz,ny,nx,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*5+4,'" Dimensions="',nz,ny,nx,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'


    write(nfil,*)'            <Attribute Name="visc" AttributeType="Vector" Center="Cell">'
    write(nfil,*)'               <DataItem Dimensions="',nz,ny,nx,'3" Function="JOIN($0, $1, $2)"'
    write(nfil,*)'                         ItemType="Function">'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*6+4,'" Dimensions="',nz,ny,nx,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*7+4,'" Dimensions="',nz,ny,nx,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'                  <DataItem  Format="Binary" DataType="Float" Precision="4" Endian="little"'
    write(nfil,*)'                     Seek="',nx1*ny1*nz1*4*6+nx*ny*nz*4*8+4,'" Dimensions="',nz,ny,nx,' 1">'
    write(nfil,*)'                     '//trim(a)//chits
    write(nfil,*)'                  </DataItem>'
    write(nfil,*)'               </DataItem>'
    write(nfil,*)'            </Attribute>'

    endif


    write(nfil,*)'        </Grid>'

  enddo
  write(nfil,'(/)')
  write(nfil,*)'    </Grid>'
  write(nfil,*)'</Domain>'
  write(nfil,'(A7)')'</Xdmf>'
  close(nfil)

  end subroutine visu
