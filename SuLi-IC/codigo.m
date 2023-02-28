%Leitura de sondas para validação dos efeitos da rugosidade no SuLi
%Referência: Zampiron et al. (2022)

%Implementação: Mariana de Cesaro

clc
clear

%ENTRADA - arquivo de leitura da sonda
[a] = dlmread("sondaxH16");

l=size(a,1);
c=size(a,2);

media=mean(a(12:l,2:c));

%GRAFICOS
figure(1)

%ax(1)=subplot(221);
x=2:5;
subplot(221);
%set (ax(1), "tag", "1");
plot(a(2:52,1),a(2:52,2),"m", a(2:52,1),a(2:52,3),"r",a(2:52,1),a(2:52,4),"b",a(2:52,1),a(2:52,5),"g");
legend("Sonda 1", "Sonda 2", "Sonda 3", "Sonda 4", "location", "northwestoutside")
grid on;
xlabel ("t (s)")
ylabel ("$\overline{u}$");


subplot(222);
plot(a(2:52,1),a(2:52,6),"c;sonda 5;",a(2:52,1),a(2:52,7),"k;sonda 6;",a(2:52,1),a(2:52,8),"r;sonda 7;",a(2:52,1),a(2:52,9),"g;sonda 8;")
legend("Sonda 5", "Sonda 6", "Sonda 7", "Sonda 8", "location", "northeastoutside");
grid on;
xlabel ("t (s)")
ylabel vmédia(m/s)

subplot(223);
plot(a(2:52,1),a(2:52,10),"m;sonda 9;",a(2:52,1),a(2:52,11),"r;sonda 10;",a(2:52,1),a(2:52,12),"b;sonda 11;",a(2:52,1),a(2:52,13),"g;sonda 12;");
legend("Sonda 9", "Sonda 10", "Sonda 11", "Sonda 12", "location", "southwestoutside");
grid on;
xlabel ("t (s)")
ylabel vmédia(m/s)

subplot(224);
plot(a(2:52,1),a(2:52,14),"c;sonda 13;",a(2:52,1),a(2:52,15),"k;sonda 14;",a(2:52,1),a(2:52,16),"r;sonda 15;",a(2:52,1),a(2:52,17),"g;sonda 16;")
legend("Sonda 13", "Sonda 14", "Sonda 15", "Sonda 16", "location", "southeastoutside");
grid on;
xlabel ("t (s)")
ylabel vmédia(m/s)

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





