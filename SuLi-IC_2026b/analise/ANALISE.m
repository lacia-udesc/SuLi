clc
clear
close all

% --- Configurações do Usuário -----------------------------------------
fileName = '../dados/sondaxH51'; % Nome do arquivo de dados    

data = readmatrix(fileName);

start_col_vel = 2;            % Ponto inicial da coluna de velocidade. O ponto 1 se refere ao tempo
end_col_vel   = size(data, 2);           % Ponto final da coluna de velocidade

dist = 2
coef = 10
dzmin = 0.001
dzmax = 0.002
for k = start_col_vel-1: end_col_vel-1
	hss(k) = 0.5 * (1.+ tanh((-dist + k/coef)));
	dz(k) = dzmin * (1.-hss(k)) + dzmax * hss(k);
		
end

nz(1) = dzmin/2-0.016;
for k = start_col_vel: end_col_vel-1
    nz(k) =  nz(k-1)+dz(k);
end

% A coluna de tempo é a primeira
tempo = data(2:end, 1);

% Extrai os dados da sonda (velocidade) - todas as linhas a partir da segunda, colunas 2 até 'aux'
velocidade = data(2:end, start_col_vel:end_col_vel);


% --- Dados Experimentais ----------------------------------------------

config.us = 0.043;   % Velocidade de atrito (m/s)
config.H  = 0.127;    % Profundidade do canal (m) (Zampiron et al., 2022)

%exp_uu = [7.132 8.221 9.265 10.216 11.031 11.755 12.388 13.066 13.562 14.058 14.554 14.913]; % Velocidades adimensionais
%exp_zH = [0.09706 0.12704 0.1768 0.24139 0.30842 0.39029 0.46965 0.55398 0.64073 0.72749 0.8192 0.91334]; % Posições verticais adimensionais

exp_uu = [7.3679 8.4784 9.5884 10.4759 11.2642 12.003 12.6182 13.1099 13.5766 13.9448 14.2881 14.5327]; % Velocidades adimensionais
exp_zH = [0.09825 0.12651 0.17093 0.22746 0.30013 0.37685 0.46299 0.55047 0.64603 0.73755 0.83311 0.92732]; % Posições verticais adimensionais


% Converter para unidades físicas
exp_u = exp_uu * config.us;  % Velocidades em m/s
exp_H = exp_zH * config.H;   % Profundidades em metros (altura a partir do fundo)

% --- Marcacao do Ibm ----------------------------------------------

ibm(1,1) = 0
ibm(2,1) = 1
ibm(1:2,2) = 0.0

% --- Importar e Processar os Dados ------------------------------------

% Verifica as dimensões
disp(['Dimensões dos tempos: ', num2str(size(tempo))]);
disp(['Dimensões dos pontos de profundidade (nz): ', num2str(size(nz))]);
disp(['Dimensões dos dados da sonda (velocidade): ', num2str(size(velocidade))]);

% Validação crucial: o número de colunas de 'velocidade' DEVE ser igual ao número de 'nz'
if size(velocidade, 2) ~= length(nz)
    % Esta validação deve passar agora se o cálculo de 'start_col_velocidade' e 'end_col_velocidade' estiver correto.
    error('Discrepância nas dimensões após extração. O número de colunas de velocidade lidas do arquivo (%d) não corresponde ao número de pontos definidos em ''nz'' (%d).', size(velocidade, 2), length(nz));
end
if size(velocidade, 1) ~= length(tempo)
    warning('O número de linhas em ''velocidade'' não corresponde ao número de pontos em ''tempo''. Verifique a extração de dados.');
end


% --- Cálculo das Médias por Intervalo de Tempo ------------------------

% Definir os intervalos de tempo para calcular as médias
intervalos_tempo = [
    0,0;
    0, 10;  % Intervalo 1
    10, 40;  % Intervalo 2
    40, 50;  % Intervalo 3
    50, 80;  % Intervalo 4
    80, 100;  % Intervalo 5
    100, 160;  % Intervalo 6 
    160, 200;  % Intervalo 7
    200, 240;  % Intervalo 8
    240, 280;  % Intervalo 9
    280, 320;  % Intervalo 10
    tempo(end-1,1)-1, tempo(end,1); % último perfil em 1 seg
];

num_intervalos = size(intervalos_tempo, 1);
medias_velocidade = zeros(num_intervalos, length(nz)); % Matriz para armazenar as médias
tempos_centrais_intervalos = mean(intervalos_tempo, 2);

fprintf('\nCalculando médias de velocidade para os intervalos:\n');
for k = 1:num_intervalos
    t_inicio = intervalos_tempo(k, 1);
    t_fim = intervalos_tempo(k, 2);
    idx_intervalo = find(tempo >= t_inicio & tempo <= t_fim);

    if isempty(idx_intervalo)
        warning('Nenhum dado encontrado para o intervalo de tempo %.2f s a %.2f s. Verifique os limites ou dt.', t_inicio, t_fim);
        medias_velocidade(k, :) = NaN;
    else
        velocidade_intervalo = velocidade(idx_intervalo, :);
        medias_velocidade(k, :) = mean(velocidade_intervalo, 1);
        fprintf('  Média no intervalo %.2f s - %.2f s calculada.\n', t_inicio, t_fim);
    end
end




%%% cálculo de erros 
if ~isnan(medias_velocidade(num_intervalos, length(nz)))
    nzz = length(nz)
elseif ~isnan(medias_velocidade(num_intervalos, length(nz)-1))
    nzz = length(nz)-1
elseif ~isnan(medias_velocidade(num_intervalos, length(nz)-2))
    nzz = length(nz)-2
elseif ~isnan(medias_velocidade(num_intervalos, length(nz)-3))
    nzz = length(nz)-3
else
    nzz = num_intervalos
end
%% fórmula gerada no calc, usando os pontos medidos e teve um R² de 0,99973
u_exp = @(y) -6202513*y^6 + 2540626*y^5 - 414670*y^4 + 34437.4*y^3 - 1549.42*y^2 + 39.383*y;

%%% começa no 10 porque é onde começa a ter dados experimentais
for k = 7:nzz
    erro(k-6) = medias_velocidade(7, k) - u_exp(nz(k));
end

erromedio = mean(erro)
erromedioq = sqrt(mean(erro.^2))


%% esse é apenas o erro inicial, teria que pensar para fazer o erro no final....
erro_prof=(nz(end_col_vel-1)-0.127)/0.127 * 100

% --- Visualização das Médias -----------------------------------------

% % 1. Velocidade Média (X) vs. Profundidade (Y)
figure('Name', 'Perfis de Velocidade Média (Simulado vs. Experimental)');
hold on;

% --- Plotar Dados Simulados (Médias dos Intervalos) ---
cores_simulado = parula(num_intervalos); % Cores para as linhas simuladas

for k = 1:num_intervalos
    % Plota o perfil médio de velocidade (eixo X) vs. profundidade (eixo Y)
    if ~all(isnan(medias_velocidade(k, :)))
        plot(medias_velocidade(k, :), nz, '-o', ...
             'Color', cores_simulado(k,:), ...
             'LineWidth', 1.5, ...
             'MarkerSize', 5, ...
             'DisplayName', sprintf('Simulado Média: %.0f - %.0f s', ...
             intervalos_tempo(k, 1), intervalos_tempo(k, 2)));
    end
end

% --- Plotar Dados Experimentais ---
% Os dados experimentais geralmente são plotados com marcadores e uma linha tracejada
% ou sólida, mas visivelmente diferente dos simulados.
plot(exp_u, exp_H, 'b--o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'w', 'DisplayName', 'Experimental');

plot(ibm(:,1),ibm(:,2), 'r--', 'LineWidth', 1.5, 'DisplayName', 'IBM');
hold off;

xlim([0 1.])
ylim([0 0.13])
xlabel('Velocidade Média (m/s)', 'FontSize', 12);
ylabel('Profundidade (m)', 'FontSize', 12);
title(sprintf('Perfis de Velocidade Média: Simulação vs. Experimental\nRMSE = %.4f m/s | Erro Rel. de prof. = %.2f%%', ...
      erromedioq, erro_prof), 'FontSize', 14);
%title('Perfis de Velocidade Média: Simulação vs. Experimental', 'FontSize', 14);

% Ajustar a legenda para incluir todos os elementos plotados
legend('Location', 'bestoutside');
grid on;
box on;

saveas(gcf, fullfile('','perfil_velocidade.png'));

fileName2 = '../dados/umed.txt';   
data2 = readmatrix(fileName2);
end_col_vel2   = size(data2, 2); 
tempo2 = data2(2:end, 1);
u2 = data2(2:end, 2);
plot(tempo2, u2, 'LineWidth',1.0);
xlabel('tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Velocidade média por tempo', 'FontSize', 14);
legend('Velocidade média');
grid on;
box on;

saveas(gcf, fullfile('','vel_med.png'));

fileName3 = '../dados/SXcentered';   
data3 = readmatrix(fileName3);
end_col_vel3   = size(data3, 2); 
tempo3 = data3(2:end, 1);
u3 = data3(2:end, 2);
v3 = data3(2:end, 3);
w3 = data3(2:end, 4);
plot(tempo3, u3, 'LineWidth',1.0);
xlabel('tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Velocidade x por tempo', 'FontSize', 14);
legend('Velocidade');
grid on;
box on;
saveas(gcf, fullfile('','vel_meio.png'));
plot(tempo3, v3,tempo3, w3, 'LineWidth',1.0);
xlabel('tempo (s)', 'FontSize', 12);
ylabel('Velocidade (m/s)', 'FontSize', 12);
title('Velocidade y e z por tempo', 'FontSize', 14);
legend('Velocidade y', 'Velocidade z');
grid on;
box on;
saveas(gcf, fullfile('','vel_meioyz.png'));