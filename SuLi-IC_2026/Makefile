#======================================================================
# Makefile for SuLi
#=======================================================================

# GNU
ifdef USE_GPU ## está desativado. retirei neste código as chamadas, mas tem no github uma versão com elas.
	FC = nvfortran
	OPTFC = -O3 -mcmodel=medium -fopenmp -DUSE_GPU -MP -cpp #-mp=gpu -fast -gpu=mem:separate # Ativa pré-processador
else
	FC = gfortran
	OPTFC = -O3 -fdefault-real-8 -fdefault-double-8 -mcmodel=medium -fopenmp -cpp # -fcheck=all 
	#-fdefault-real-8 -fdefault-double-8
endif

CC = #scorep


# Intel
#FC = ifort
#OPTFC =  -fpp -O3 -autodouble -mcmodel=medium -qopenmp
#CC = icc

FLAGS = -O3

# FFT-based Poisson solver for t_press=2 (requires FFTW3 dev package)
# Most distros install fftw3.f03 in /usr/include
FFTW_INC ?= /usr/include
OPTFC += -I$(FFTW_INC)
FFTW_LIBS = -lfftw3 -lm
 
# List of source files
SRC =   0_modules.f90 1_principal.f90 2_iniciais.f90 3_contorno.f90 4_plot.f90 5_viscosidade.f90 6_convdiff.f90 7_advectivo.f90 8_derivadas.f90 9_grad.f90 10_pos.f90 11_diversos.f90 12_paraview_3D.f90  14_level_set.f90 15_ibm.f90 # 13_mms.f90 MMS DESATUALIZADO

#-----------------------------------------------------------------------
# Normally no need to change anything below

OBJ =	$(SRC:.f90=.o)

rodar_3d : $(OBJ)
	$(CC) $(FC) -O3 -o $@ $(OBJ) $(OPTFC) $(FFTW_LIBS)

%.o : %.f90
	$(CC) $(FC) $(OPTFC) $(OPTIONS) -c $<

.PHONY: clean
clean:
	rm -f *.o *.mod 

allclean:
	rm -f *.o *.mod rodar_3d
	rm -rf arquivos
	rm -rf dados
	rm -f saida
	rm -rf fort.*
