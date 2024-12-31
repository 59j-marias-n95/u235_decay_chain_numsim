%clear
%clc
format longe;

%tiempo0 = clock();
#Definase la matriz M
lambda=[9.846e-10, 1.786e-3, 2.116e-5, 3.184e-2, 1.353e+1, 1.656e+4, 2.213e+1, 5.520e+6, 1.227e+10, 1.009e+4, 2.186e+11, 1.702e+5, 4.236e+7, 7.638e+4, 0];

#Las fracciones de desintegración para los canales alfa y beta en las etapas de
#desintegración correspondientes al Actinio 227, Polonio 2015 y Bismuto 211.
ratio=[0.013800, 0.98620, 0.9999977, 0.0000023, 0.99724, 0.00276];

#La matriz del sistema de ecuaciones de orden 15 está compuesto por el siguiente
#arreglo conformado por los valores del vector lambda y oportunamente incluyendo
#los del vector ratio
M=zeros(15,15);
M(1,1)=-lambda(1); M(2,2)=-lambda(2); M(3,3)=-lambda(3); M(4,4)=-lambda(4);
M(5,5)=-lambda(5); M(6,6)=-lambda(6); M(7,7)=-lambda(7); M(8,8)=-lambda(8);
M(9,9)=-lambda(9); M(10,10)=-lambda(10); M(11,11)=-lambda(11);
M(12,12)=-lambda(12); M(13,13)=-lambda(13); M(14,14)=-lambda(14);

M(2,1)=lambda(1); M(3,2)=lambda(2); M(4,3)=lambda(3); M(5,4)=ratio(1)*lambda(4);
M(6,4)=ratio(2)*lambda(4); M(7,5)=lambda(5); M(7,6)=lambda(6);
M(8,7)=lambda(7); M(9,8)=lambda(8); M(10,9)=ratio(3)*lambda(9);
M(11,9)=ratio(4)*lambda(9); M(12,10)=lambda(10); M(12,11)=lambda(11);
M(13,12)=ratio(5)*lambda(12); M(14,12)=ratio(6)*lambda(12); M(15,13)=lambda(13);
M(15,14)=lambda(14);

%El vector de transiciones
DeltaN = [[-1, 1, zeros(1, 13)]', [0, -1, 1, zeros(1, 12)]', [zeros(1, 2), -1, 1,zeros(1, 11)]', [zeros(1, 3), -1, 1, 1,zeros(1, 9)]',[zeros(1, 4), -1, 0, 1, zeros(1, 8)]', [zeros(1, 5), -1, 1, zeros(1, 8)]', [zeros(1, 6), -1, 1, zeros(1, 7)]', [zeros(1, 7), -1, 1, zeros(1, 6)]',[zeros(1, 8), -1, 1, 1, zeros(1, 4)]', [zeros(1, 9), -1, 0, 1, zeros(1, 3)]', [zeros(1, 10), -1, 1, zeros(1, 3)]', [zeros(1, 11), -1, 1, 1, 0]',[zeros(1, 12), -1, 0, 1]', [zeros(1, 13), -1, 1]', zeros(1, 15)'];

%El numero inicial de particulas
%N1=(7.2331)*(1/238.03)*(6.022e+23);
N = [[(7.2331e+03)*(1/238.03)*(6.022e+23), zeros(1,14)]', zeros(15, 1)];
P = zeros(15, 1);
Sigmat = zeros(15, 15);

%La escala de tiempo
m = 2.00e+03;
deltat = 2.0e+02;
t = 0:deltat:m*deltat;

%La matriz de evolucion del sistema
for a = 1:m
  %Probabilidad de que la kesima sustancia se desintegre
  for k = 1:15
    P(k,1) = lambda(k)*N(k,a)*deltat;
  endfor

  %Calculo de la covarianza matriz
  for k = 1:15
    Sigmat = Sigmat + P(k,1)*DeltaN(1:15, k)*DeltaN(1:15, k)';
  endfor
  [sigmat, X, Y] = lu(Sigmat);
  sigmat = sqrt(deltat)*sigmat;

  %El ruido blanco
  DeltaW = sqrt(deltat)*(rand(15,1));

  for k = 1:a
    N(1:15, k+1) = (eye(15) + deltat*M)*N(1:15, k) + sigmat*DeltaW; N(N<=0)=0;
  endfor
endfor

N2 = [[(7.2331e+03)*(1/238.03)*(6.022e+23), zeros(1,14)]', zeros(15, 1)];
P2 = zeros(15, 1);
Sigmat = zeros(15, 15);
for a = 1:m
  %Probabilidad de que la kesima sustancia se desintegre
  for k = 1:15
    P2(k,1) = lambda(k)*N2(k,a)*deltat;
  endfor

  %Calculo de la covarianza matriz
  for k = 1:15
    Sigmat = Sigmat + P2(k,1)*DeltaN(1:15, k)*DeltaN(1:15, k)';
  endfor
  [sigmat, X, Y] = lu(Sigmat);
  sigmat = sqrt(deltat)*sigmat;

  %El ruido blanco
  DeltaW = sqrt(deltat)*(rand(15,1));

  for k = 1:a
    N2(1:15, k+1) = (eye(15) + deltat*M)*N2(1:15, k) + sigmat*DeltaW;
    N2(N2<=0)=0;
  endfor
endfor

N3 = [[(7.2331e+03)*(1/238.03)*(6.022e+23), zeros(1,14)]', zeros(15, 1)];
P3 = zeros(15, 1);
Sigmat = zeros(15, 15);
for a = 1:m
  %Probabilidad de que la kesima sustancia se desintegre
  for k = 1:15
    P3(k,1) = lambda(k)*N3(k,a)*deltat;
  endfor

  %Calculo de la covarianza matriz
  for k = 1:15
    Sigmat = Sigmat + P3(k,1)*DeltaN(1:15, k)*DeltaN(1:15, k)';
  endfor
  [sigmat, X, Y] = lu(Sigmat);
  sigmat = sqrt(deltat)*sigmat;

  %El ruido blanco
  DeltaW = sqrt(deltat)*(rand(15,1));

  for k = 1:a
    N3(1:15, k+1) = (eye(15) + deltat*M)*N3(1:15, k) + sigmat*DeltaW;
    N3(N3<=0)=0;
  endfor
endfor

%Figura 1
figure(1,"position",[0, 0, 1000, 800])
%hold on
%subplot(2,1,1)
semilogx(N(1, 1:m),'Color', [1 0 0],'LineWidth',3);
hold on;
semilogx(N2(1, 1:m), 'Color', [1 1 0], 'LineWidth', 2);
semilogx(N3(1, 1:m), 'Color', [0 1 1], 'LineWidth', 1);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{1}: U235')

%{
subplot(2,1,2)
semilogx(t(1:m), N(2, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{2}: Th231')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 1 y 2.', 'fontsize', 14);

%Figura 2
figure(2,"position",[0, 0, 1200, 1000])
hold on
subplot(2,1,1)
semilogx(t(1:m), N(3, 1:m), 'Color', [0 0 1], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 14)
ylabel('N(t)/núcleos', 'fontsize', 14)
legend('N_{3}: Pa231')

subplot(2,1,2)
semilogx(t(1:m), N(4, 1:m), 'Color', [1 0 1], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{4}: Ac227')
figure2 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 3 y 4.', 'fontsize', 14);
hold off

%Figura 3
figure(3,"position",[0, 0, 1200, 1000])
hold on
subplot(2,1,1)
semilogx(t(1:m), N(5, 1:m), 'Color', [1 0 1], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{5}: Th227')

subplot(2,1,2)
semilogx(t(1:m), N(6, 1:m), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{6}: Fr223')
figure3 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 5 y 6.', 'fontsize', 14);
hold off

%Figura 4
figure(4,"position",[0, 0, 1200, 1000])
hold on
subplot(2,1,1)
semilogx(t(1:m), N(7, 1:m), 'Color', [0 0 0], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 14)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{7}: Ra223')

subplot(2,1,2)
semilogx(t(1:m), N(8, 1:m), 'Color', [0.6350 0.0780 0.1840], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 14)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{8}: Rn219')
figure4 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 7 y 8.', 'fontsize', 14);
hold off

%Figura 5
figure(5,"position",[0, 0, 1200, 1000])
hold on
subplot(2,1,1)
semilogx(t(1:m), N(9, 1:m), 'Color', [0.3010 0.7450 0.9330], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 14)
ylabel('N(t)/núcleos', 'fontsize', 14)
legend('N_{9}: Po215')

subplot(2,1,2)
semilogx(t(1:m), N(10, 1:m), 'Color', [0.4660 0.6740 0.1880], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{10}: Pb211')
figure5 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 9 y 10.', 'fontsize', 14)
hold off

%Figura 6
figure(6,"position",[0, 0, 1200, 1000])
hold on
subplot(2,1,1)
semilogx(t(1:m), N(11, 1:m), 'Color', [0.4940 0.1840 0.5560], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{11}: At215')

subplot(2,1,2)
semilogx(t(1:m), N(12, 1:m), 'Color', [0.9290 0.6940 0.1250], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{12}: Bi211')
figure6 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 11 y 12.', 'fontsize', 14)
hold off

%Figura 7
figure(7,"position",[0, 0, 1200, 1000])
hold on
subplot(2,1,1)
semilogx(t(1:m), N(13, 1:m), 'Color', [0.8500 0.3250 0.0980], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{13}: Po211')

subplot(2,1,2)
semilogx(t(1:m), N(14, 1:m), 'Color', [0 0.4470 0.7410], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{14}: Tl207')
figure7 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 13 y 14.', 'fontsize', 14)
hold off

%Figura 8
figure(8,"position",[0, 0, 1200, 1000])
subplot(1,1,1)
semilogx(t(1:m), N(15, 1:m), 'Color', [1 0 0], 'LineWidth', 2)
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{15}: Pb207')
figure8 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para N_{15}: Pb207.', 'fontsize', 14)
}%
