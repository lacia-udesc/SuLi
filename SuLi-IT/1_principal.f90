! SULI 1.0
! Caso você utilize o SuLi para algum trabalho ou publicação por favor citar os seguintes trabalhos:
!MONTEIRO, L. R.; SCHETTINI, E. B. C. . Comparação entre a aproximação hidrostática e a não-hidrostática na simulação numérica de escoamentos com superfície livre. Revista Brasileira de Recursos Hídricos, v. 20, p. 1051-1062, 2015.

!MONTEIRO, L. R. Simulação numérica de escoamentos com superfície livre com aproximação não-hidrostática. 2014. 94f. Dissertação de Mestrado, Programa de Pós-Graduação em Engenharia de Recursos Hídricos e Saneamento Ambiental da Universidade Federal do Rio Grande do Sul, 2014.

! Não nos responsabilizamos pelos resultados gerados pelo código.
  
  
  
!!!! Cálculo da dinâmica de fluidos (escoamento ou propagação de ondas) !!!!
! Baseado na Dissertação de Mestrado de Leonardo Romero Monteiro (abril de 2014)

!!!! Descrição do método !!!!
! Aproximação de pressão não-hidrostática utilizando o método de fracionamento do tempo (1ª Parte = Hidrostática, 2ª Parte = Não-Hidrostática);
! Ambiente tri-dimensional, x, y (horizontais) e z (vertical);
! Método semi-implícito e com correção quasi-implícita (theta) em diferenças finitas com a grade deslocada;

!!! Implementação 15/11/2014
! Leonardo Romero Monteiro

!!! Última modificação
! Leonardo Romero Monteiro - 08/08/2020

PROGRAM PNH
!===================================================================================================================
!!!! Declaração de Variáveis !!!!
USE disc

IMPLICIT NONE

!Condições iniciais
CALL iniciais()

!Adicionar os contornos na plotagem inicial
CALl contorno()

!Solução manufaturada
if (mms_t > 0) CALL mms()


!Plotagens iniciais
CALL plot_i()
!CALL analise_i()
!===================================================================================================================
!RESOLUÇÃO DO PROBLEMA
!===================================================================================================================
!Parte 1: visco() - cálculo da viscosidade turbulenta
!Parte 2: convdiff() - Cálculo das velocidades não corrigidas
!Parte 3: convdiff() - Cálculo do Desnível pelo Gradiente Conjugado (não corrigido)
!Parte 4: convdiff() - Cálculo das novas velocidades e desnível (hidrosático) 
!Parte 5: graddin e posdin - Cálculo das pressão dinâmica e correção do desnível
!Parte 6: Plotagens (plot_f)

do it = 1, ts
t = it * dt

	! Termo fonte para o método da sulução manufaturada (MMS)
	if ((mms_t == 1) .and. (it == 1)) call termo_fonte1()
	if (mms_t == 2) call termo_fonte2()

	CALL visco()
	
	CALL convdiff()

	if (wave_t > 0) call boundary_waves() ! for wave propagation
	!Condições de Contorno para a parte Hidrostática
	CALl contorno()

	if (mms_t .eq. 0) then
	 CALL graddin()
	 CALL posdin()
	else
	 !Solução manufaturada; cálculo do erro
	 CALL mms()
	endif

	CALL plot_f()
	!CALL analise_f()

enddo

CALL plot_atrib()

!em plot.f90
 close (unit=100001)
 close (unit=100002)
 close (unit=200000)
 close (unit=9999991)

!==================================================================================================================
End program PNH

