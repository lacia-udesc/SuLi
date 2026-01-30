clc
clear
close all

% --- Configurações do Usuário -----------------------------------------

fileName = 'dados/sondaxH51'; % Nome do arquivo de dados                   
dt = 0.002;                   % Intervalo de tempo entre as amostras no arquivo (informado pelo usuário)
%nz_full = 0.008:0.008:0.127;  % Posições verticais simuladas (altura a partir do fundo, em metros)
%nz = nz_full(4:end);          % Ignora os 3 primeiros pontos
start_col_vel = 2;            % Coluna inicial de velocidade
end_col_vel   = 28;           % Última coluna de velocidade

dist = 2
coef = 10
dzmin = 0.002
dzmax = 0.01
for k = start_col_vel: end_col_vel
	hss(k) = 0.5 * (1.+ tanh((-dist + k/coef)));
	dz(k) = dzmin * (1.-hss(k)) + dzmax * hss(k);
		
end

nz(start_col_vel) = 0;
for k = (start_col_vel+1): end_col_vel-1
    nz(k) =  nz(k-1)+dz(k);
end

% Verificação
try
    data = readmatrix(fileName);
catch ME
    error('Erro ao ler o arquivo %s. Verifique o formato ou o caminho. Erro: %s', fileName, ME.message);
end

if end_col_vel > size(data, 2)
    error('As colunas de velocidade a serem lidas (%d a %d) excedem o número real de colunas no arquivo (%d). Ajuste ''nz'' ou o número de colunas a ignorar.', start_col_velocidade, end_col_velocidade, size(full_data, 2));
end

% A coluna de tempo é a primeira
tempo = data(2:end, 1);

% Extrai os dados da sonda (velocidade) - todas as linhas a partir da segunda, colunas 2 até 'aux'
velocidade = data(2:end, start_col_vel:end_col_vel);


% --- Dados Experimentais ----------------------------------------------

config.us = 0.043;   % Velocidade de atrito (m/s)
config.H  = 0.127;    % Profundidade do canal (m) (Zampiron et al., 2022)

exp_uu = [7.132 8.221 9.265 10.216 11.031 11.755 12.388 13.066 13.562 14.058 14.554 14.913]; % Velocidades adimensionais
exp_zH = [0.09706 0.12704 0.1768 0.24139 0.30842 0.39029 0.46965 0.55398 0.64073 0.72749 0.8192 0.91334]; % Posições verticais adimensionais

% Converter para unidades físicas
exp_u = exp_uu * config.us;  % Velocidades em m/s
exp_H = exp_zH * config.H;   % Profundidades em metros (altura a partir do fundo)


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
    0, 10;  % Intervalo 1: 60s a 70
    10, 20;  % Intervalo 1: 60s a 70s
    20, 30;  % Intervalo 1: 60s a 70s
    60, 70;  % Intervalo 1: 60s a 70s
    70, 80;  % Intervalo 2: 70s a 80s
    80, 90;  % Intervalo 3: 80s a 90s
    90, 100;  % Intervalo 4: 90s a 100s
    100, 110;  % Intervalo 5: 100s a 110s
    110, 120;  % Intervalo 6: 110s a 120s
    120, 130;  % Intervalo 7: 120s a 130s
    130, 140;  % Intervalo 8: 130s a 140s
    140, 150;  % Intervalo 9: 140s a 145s
    100, 150;  % Intervalo Geral
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
if ~isnan(medias_velocidade(14, length(nz)))
    nzz = length(nz)
elseif ~isnan(medias_velocidade(14, length(nz)-1))
    nzz = length(nz)-1
elseif ~isnan(medias_velocidade(14, length(nz)-2))
    nzz = length(nz)-2
elseif ~isnan(medias_velocidade(14, length(nz)-3))
    nzz = length(nz)-3
else
    nzz = 9
end
%% fórmula gerada no calc, usando os pontos medidos e teve um R² de 0,99973
u_exp = @(y) -6202513*y^6 + 2540626*y^5 - 414670*y^4 + 34437.4*y^3 - 1549.42*y^2 + 39.383*y;

%%% começa no 9 porque é onde começa a ter dados experimentais
for k = 9:nzz
    u_e(k) = u_exp(nz(k));
    erro(k-8) = medias_velocidade(14, k) - u_e(k);
end

erromedio = mean(erro)
erromedioq = sqrt(mean(erro.^2))


%% esse é apenas o erro inicial, teria que pensar para fazer o erro no final....
erro_prof=(nz(nzz)-exp_H(12))/exp_H(12) * 100

% --- Visualização das Médias -----------------------------------------

% % 1. Velocidade Média (X) vs. Profundidade (Y)
figure('Name', 'Perfis de Velocidade Média (Simulado vs. Experimental)');
hold on;

% --- Plotar Dados Simulados (Médias dos Intervalos) ---
cores_simulado = parula(num_intervalos); % Cores para as linhas simuladas

for k = 1:num_intervalos
    % Plota o perfil médio de velocidade (eixo X) vs. profundidade (eixo Y)
    if ~all(isnan(medias_velocidade(k, :)))
        plot(medias_velocidade(k, :), nz, ...
             'Color', cores_simulado(k,:), 'LineWidth', 1.5, ...
             'DisplayName', sprintf('Simulado Média: %.0f - %.0f s', intervalos_tempo(k, 1), intervalos_tempo(k, 2)));
    end
end

% --- Plotar Dados Experimentais ---
% Os dados experimentais geralmente são plotados com marcadores e uma linha tracejada
% ou sólida, mas visivelmente diferente dos simulados.
plot(exp_u, exp_H, 'k--o', 'LineWidth', 1.5, 'MarkerSize', 6, ...
     'MarkerFaceColor', 'k', 'DisplayName', 'Experimental');

hold off;

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

% % 2. Tempo (X) vs. Velocidade Média (Y)
% figure('Name', 'Evolução da Velocidade Média por Profundidade ao Longo do Tempo');
% hold on;
% 
% num_dz_points = length(dz);
% cores = jet(num_dz_points); % Cores para cada ponto de profundidade
% 
% % Itera sobre cada ponto de profundidade (cada coluna de medias_velocidade)
% for i = 1:num_dz_points
%     % Plotamos os tempos centrais dos intervalos no Eixo X
%     % e as médias de velocidade para aquele dz específico no Eixo Y
%     plot(tempos_centrais_intervalos, medias_velocidade(:, i), ...
%          'Color', cores(i,:), 'LineWidth', 1.5, ...
%          'Marker', 'o', 'MarkerSize', 5, ... % Adiciona marcadores para cada ponto de intervalo
%          'DisplayName', sprintf('z = %.2f', dz(i))); % Legenda para cada dz
% end
% 
% hold off;
% xlabel('Tempo Central do Intervalo (s)', 'FontSize', 12);
% ylabel('Velocidade Média (m/s)', 'FontSize', 12);
% title('Evolução da Velocidade Média ao Longo dos Intervalos de Tempo para Cada Profundidade', 'FontSize', 14);
% legend('Location', 'bestoutside', 'NumColumns', 2); % Melhora a visualização da legenda
% grid on;
% box on;


% --- Visualização dos Dados -------------------------------------------

% % 1. Heatmap (Mapa de Cores) da Profundidade vs. Tempo
% figure('Name', 'Heatmap da Velocidade'); % Adiciona um nome à figura
% imagesc(dz, tempo, velocidade); % dz no eixo X, tempo no eixo Y
% colorbar; % Adiciona barra de cores
% xlabel('Profundidade (dz) [unidades]', 'FontSize', 12); % Unidades apropriadas
% ylabel('Tempo (s)', 'FontSize', 12);
% title('Distribuição da Velocidade ao Longo da Profundidade e Tempo', 'FontSize', 14);
% set(gca, 'YDir', 'normal'); % Inverte o eixo Y para que o tempo cresça para cima (mais intuitivo)
% colormap parula; % Colormap padrão do MATLAB (muitas vezes melhor que 'jet')
% % Você também pode experimentar: 'viridis', 'magma', 'inferno', 'plasma', 'hot'
% axis tight; % Ajusta os limites dos eixos para os dados

% % 2. Tempo (X) vs. Velocidade (Y) para Todos os Pontos de Profundidade ---
% 
% figure('Name', 'Velocidade vs. Tempo para cada Profundidade');
% hold on;
% 
% num_dz_points = length(dz);
% cores = jet(num_dz_points); % Usa um colormap para cores distintas para cada linha
% 
% for i = 1:num_dz_points
%     % Plotamos a velocidade (coluna i) no eixo X e o tempo no eixo Y
%     plot(tempo, velocidade(:, i), 'Color', cores(i,:), 'LineWidth', 1.2, ...
%          'DisplayName', sprintf('z = %.2f', dz(i)));
% end
% 
% hold off;
% xlabel('Tempo (s)', 'FontSize', 12); % Eixo X agora é velocidade
% ylabel('Velocidade (m/s)', 'FontSize', 12); % Eixo Y agora é tempo
% title('Velocidade ao Longo do Tempo para Diferentes Profundidades', 'FontSize', 14);
% legend('Location', 'bestoutside', 'NumColumns', 2); % Melhora a legenda para muitos pontos
% grid on;
% box on;
