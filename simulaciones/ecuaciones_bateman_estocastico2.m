clear
clc
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

%La escala de tiempo
m = 2.00e+03;
deltat = 5.0e+01;

%El numero inicial de particulas: (7.2331e+03)*(1/238.03)*(6.022e+23)
N0 = 1e+9;
N = [[N0, zeros(1,14)]', zeros(15, 1)];
P = zeros(15, 1);
Sigmat = zeros(15, 15);
DeltaW = zeros(15, 1);
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
  for k = 1:a
    N(1:15, k+1) = pinv(eye(15) - deltat*M)*(N(1:15, k) + sigmat*DeltaW);
    DeltaW = sqrt(deltat)*(rand(15,1));
  endfor
endfor

N2 = [[N0, zeros(1,14)]', zeros(15, 1)];
P2 = zeros(15, 1);
Sigmat = zeros(15, 15);
DeltaW = zeros(15, 1);
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
  for k = 1:a
    N2(1:15, k+1) = pinv(eye(15) - deltat*M)*(N2(1:15, k) + sigmat*DeltaW);
    DeltaW = sqrt(deltat)*(rand(15,1));
  endfor
endfor

N3 = [[N0, zeros(1,14)]', zeros(15, 1)];
P3 = zeros(15, 1);
Sigmat = zeros(15, 15);
DeltaW = zeros(15, 1);
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
  for k = 1:a
    N3(1:15, k+1) = pinv(eye(15) - deltat*M)*(N3(1:15, k) + sigmat*DeltaW);
    DeltaW = sqrt(deltat)*(rand(15,1));
  endfor
endfor

%Figura 1
figure(1,"position",[0, 0, 800, 600])
%hold on
subplot(2,1,1)
semilogx(N(1,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(1,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(1,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{1}: U235', 'N2_{1}: U235', 'N3_{1}: U235')

%
subplot(2,1,2)
semilogx(N(2, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
hold on;
semilogx(N2(2,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(2,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{2}: U235', 'N2_{2}: U235', 'N3_{2}: U235')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 1 y 2.', 'fontsize', 14);

%Figura 2
figure(2,"position",[0, 0, 800, 600])
%hold on
subplot(2,1,1)
semilogx(N(3,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(3,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(3,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{3}: U235', 'N2_{3}: U235', 'N3_{3}: U235')

%
subplot(2,1,2)
semilogx(N(4, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
hold on;
semilogx(N2(4,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(4,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{4}: U235', 'N2_{4}: U235', 'N3_{4}: U235')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 3y 4.', 'fontsize', 14);

%Figura 3
figure(3,"position",[0, 0, 800, 600])
%hold on
subplot(2,1,1)
semilogx(N(5,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(5,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(5,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{5}: U235', 'N2_{5}: U235', 'N3_{5}: U235')

%
subplot(2,1,2)
semilogx(N(6, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
hold on;
semilogx(N2(6,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(6,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{6}: U235', 'N2_{6}: U235', 'N3_{6}: U235')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras5 y 6.', 'fontsize', 14);

%Figura 4
figure(4,"position",[0, 0, 800, 600])
%hold on
subplot(2,1,1)
semilogx(N(7,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(7,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(7,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{7}: U235', 'N2_{7}: U235', 'N3_{7}: U235')

%
subplot(2,1,2)
semilogx(N(8, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
hold on;
semilogx(N2(8,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(8,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{8}: U235', 'N2_{8}: U235', 'N3_{8}: U235')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras7 y 8.', 'fontsize', 14);

%Figura 5
figure(5,"position",[0, 0, 800, 600])
%hold on
subplot(2,1,1)
semilogx(N(9,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(9,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(9,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{9}: U235', 'N2_{9}: U235', 'N3_{9}: U235')

%
subplot(2,1,2)
semilogx(N(10, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
hold on;
semilogx(N2(10,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(10,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{10}: U235', 'N2_{10}: U235', 'N3_{10}: U235')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras9 y 10.', 'fontsize', 14);

%Figura 6
figure(6,"position",[0, 0, 800, 600])
%hold on
subplot(2,1,1)
semilogx(N(11,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(11,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(11,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{11}: U235', 'N2_{11}: U235', 'N3_{11}: U235')

%
subplot(2,1,2)
semilogx(N(12, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
hold on;
semilogx(N2(12,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(12,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{12}: U235', 'N2_{12}: U235', 'N3_{12}: U235')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras11 y 12.', 'fontsize', 14);

%Figura 7
figure(7,"position",[0, 0, 800, 600])
%hold on
subplot(2,1,1)
semilogx(N(13,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(13,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(13,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{13}: U235', 'N2_{13}: U235', 'N3_{13}: U235')

%
subplot(2,1,2)
semilogx(N(14, 1:m), 'Color', [0 1 0], 'LineWidth', 2)
hold on;
semilogx(N2(14,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(14,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{14}: U235', 'N2_{14}: U235', 'N3_{14}: U235')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras13 y 14.', 'fontsize', 14);

%Figura 7
figure(8,"position",[0, 0, 800, 600])
%hold on
semilogx(N(15,:),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(15,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(15,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N_{15}: U235', 'N2_{15}: U235', 'N3_{15}: U235')

