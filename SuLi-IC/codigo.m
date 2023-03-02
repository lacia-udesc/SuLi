%Leitura de sondas para validação dos efeitos da rugosidade no SuLi
%Referência: Zampiron et al. (2022)

%Implementação: Mariana de Cesaro
%Modificações: Bruna Soares

clc
clear

%ENTRADA - arquivo de leitura da sonda
[a] = dlmread("sondaxH16");

l=size(a,1); %numero de linhas    #primeira linha: cabeçalho
c=size(a,2); %numero de colunas   #primeira coluna: tempo

%media=mean(a(12:l,2:c));

%GRAFICOS
figure(1)

%ax(1)=subplot(221);
%x=2:5;
subplot(221);
%set (ax(1), "tag", "1");
plot(a(2:52,1),a(2:52,2),"m", a(2:52,1),a(2:52,3),"r",a(2:52,1),a(2:52,4),"b",a(2:52,1),a(2:52,5),"g");
legend("Sonda 1", "Sonda 2", "Sonda 3", "Sonda 4", "location", "northwestoutside")
grid on;
xlabel ("t (s)");
ylabel ("Vmedia(m/s)");

subplot(222);
plot(a(2:52,1),a(2:52,6),"c",a(2:52,1),a(2:52,7),"k",a(2:52,1),a(2:52,8),"r",a(2:52,1),a(2:52,9),"g")
legend("Sonda 5", "Sonda 6", "Sonda 7", "Sonda 8", "location", "northeastoutside");
grid on;
xlabel ("t (s)");
ylabel ("Vmedia(m/s)");

subplot(223);
plot(a(2:52,1),a(2:52,10),"m",a(2:52,1),a(2:52,11),"r",a(2:52,1),a(2:52,12),"b",a(2:52,1),a(2:52,13),"g");
legend("Sonda 9", "Sonda 10", "Sonda 11", "Sonda 12", "location", "southwestoutside");
grid on;
xlabel ("t (s)");
ylabel ("Vmedia(m/s)");

subplot(224);
plot(a(2:52,1),a(2:52,14),"c",a(2:52,1),a(2:52,15),"k",a(2:52,1),a(2:52,16),"r",a(2:52,1),a(2:52,17),"g")
legend("Sonda 13", "Sonda 14", "Sonda 15", "Sonda 16", "location", "southeastoutside");
grid on;
xlabel ("t (s)");
ylabel ("Vmedia(m/s)");

print -f1 figure1.pdf

%figure(2)
%subplot(221);
%subplot(222);
%subplot(223);
%subplot(224);
%print -f1 figure2.pdf

%figure(3)
%subplot(221);
%subplot(222);
%plot();
%subplot(223);
%plot();
%subplot(224);
%plot();
%print -f1 figure3.pdf

%figure(4)
%subplot(221);
%;
%subplot(222);
%plot();
%subplot(223);
%plot();
%subplot(224);
%plot();
%print -f1 figure4.pdf





