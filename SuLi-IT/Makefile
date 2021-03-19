#======================================================================
# Makefile for SuLi
#=======================================================================

# GNU
FC = gfortran
OPTFC = -O3 -fdefault-real-8 -fdefault-double-8 -mcmodel=large # -fcheck=all # -floop-optimize -funroll-loops -ftree-fre -finline-functions -fpeel-loops -funsafe-math-optimizations -fprefetch-loop-arrays 

CC = gcc
FLAGS = -O3

# List of source files
SRC =   0_modules.f90 1_principal.f90 2_iniciais.f90 3_contorno.f90 4_plot.f90 5_viscosidade.f90 6_convdiff.f90 7_advectivo.f90 8_derivadas.f90 9_grad.f90 10_pos.f90 11_diversos.f90 12_paraview_3D.f90 13_mms.f90

#-----------------------------------------------------------------------
# Normally no need to change anything below

OBJ =	$(SRC:.f90=.o)

rodar_3d : $(OBJ)
	$(FC) -O3 -o $@ $(OBJ) $(LIBS)

%.o : %.f90
	$(FC) $(OPTFC) $(OPTIONS) -c $<

.PHONY: clean
clean:
	rm -f *.o *.mod 

allclean:
	rm -f *.o *.mod rodar_3d
	rm -rf arquivos
	rm -rf dados
	rm -f saida

