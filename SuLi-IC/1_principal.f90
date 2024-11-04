!SULI 1.0
!Caso você utilize o SuLi para algum trabalho ou publicação por favor citar os seguintes trabalhos:
!MONTEIRO, L. R.; SCHETTINI, E. B. C. . Comparação entre a aproximação hidrostática e a não-hidrostática na simulação numérica de escoamentos com superfície livre. Revista Brasileira de Recursos Hídricos, v. 20, p. 1051-1062, 2015.

!MONTEIRO, L. R. Simulação numérica de escoamentos com superfície livre com aproximação não-hidrostática. 2014. 94f. Dissertação de Mestrado, Programa de Pós-Graduação em Engenharia de Recursos Hídricos e Saneamento Ambiental da Universidade Federal do Rio Grande do Sul, 2014.

!Não nos responsabilizamos pelos resultados gerados
  
!Cálculo da dinâmica de fluidos (escoamento ou propagação de ondas)
!Baseado na Dissertação de Mestrado de Leonardo Romero Monteiro (abril de 2014)

!Descrição do método
!Aproximação de pressão não-hidrostática utilizando o método de fracionamento do tempo (1ª Parte = Hidrostática, 2ª Parte = Não Hidrostática);
!Ambiente tri-dimensional, x, y (horizontais) e z (vertical);
!Método semi implícito e com correção quasi-implícita (theta) em diferenças finitas com a grade deslocada;

!!!Implementação em 15/11/2014
!Leonardo Romero Monteiro

!!!Modificações
!Leonardo Romero Monteiro em 15/05/2015

PROGRAM PNH

!Declaração de Variáveis!
	!$ USE omp_lib
	USE disc, only : it, ts, t, dt, t_press, wave_t, mms_t, a_dt, numb_threads, tt, ntt
	USE restart, only : irest, interv_rest

	IMPLICIT NONE
	
	!Condições iniciais
	CALL parametros()
	
	if (irest.eq.0) then
		CALL iniciais()
	else
		CALL restart_ini()
	endif
	

	!$ CALL OMP_set_num_threads(numb_threads)
	!$OMP PARALLEL
	numb_threads = OMP_GET_NUM_THREADS()
	!$OMP END PARALLEL
	!$ write(*,'(A,I2)') "CÓDIGO EM PARALELO, quantidade de núcleos utilizados: ", numb_threads
	!$ CALL OMP_set_dynamic(.FALSE.)
	!$ CALL OMP_set_nested(.FALSE.)

	!Adicionar os contornos na plotagem inicial
	CALL contorno(1)
	CALL contorno(2)
	CALL contorno(3)

	!Solução manufaturada
	if (mms_t > 0) CALL mms()

	!Plotagens iniciais
	CALL PLOT()

	!RESOLUÇÃO DO PROBLEMA

	!Parte 1: Função distância; level_set()
	!Parte 2: Viscosidade Turbulenta; visco()
	!Parte 3: Passo preditor; convdiff() e tempo()
	!Parte 4 : Condições de contorno; call boundary_waves() e contorno()
	!Parte 5: Passo corretor; graddin() e posdin()

	do it = 1, ts
		t = t + dt

		!write(*,*) it

		!Termo fonte para o método da sulução manufaturada (MMS)
		if ((mms_t == 1) .and. (it == 1)) call termo_fonte1()
		if (mms_t == 2) call termo_fonte2()


		CALL level_set()
		CALL contorno(3)
		
		if ((t_press .eq. 2) .and. (mms_t .eq. 0)) CALL exp_press()
				
		CALL visco()
			
		do tt = 1, ntt  !inicio loop RK
			dt = a_dt(tt)
			
			CALL convdiff()
			CALL complementos()
			CALL tempo()

			if (wave_t > 0) call boundary_waves() !For wave propagation

			CALL contorno(2)

			if (mms_t .eq. 0) then
				if (t_press .eq. 0) then
					CALL pressh()	!Condições de Contorno para a parte Hidrostática
					CALL posdin()
				elseif (t_press .eq. 1) then
					CALL graddin()
					CALL posdin()
				endif
				
				CALl contorno(1)
						
			endif			
		 enddo !fim loop RK

		if ((t_press .eq. 2) .and. (mms_t .eq. 0)) CALL graddin()
		

		!Solução manufaturada; cálculo do erro
		if (mms_t > 0) CALL mms()

		!Plotagens por passo de tempo
		CALL PLOT()

		if (mod(it,ceiling(interv_rest/dt)).eq.0) then
			CALL restart_salva()
		endif

	enddo

End program PNH

