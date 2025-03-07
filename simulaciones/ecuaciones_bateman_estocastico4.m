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
m = 1.0e3;
deltat = 2.0e2;

%El numero inicial de particulas: (7.2331e+03)*(1/238.03)*(6.022e+23)
%N0 = (7.2331e+03)*(1/238.03)*(6.022e+23);
N0 = 9.9999e+14;
N = [[N0, zeros(1,14)]', zeros(15, 1)];
N2 = [[N0, zeros(1,14)]', zeros(15, 1)];
N3 = [[N0, zeros(1,14)]', zeros(15, 1)];
DeltaW = zeros(15, 2);
DeltaW2 = zeros(15, 2);
DeltaW3 = zeros(15, 2);

%La matriz de evolucion del sistema
%La matriz de evolucion del sistema
for a = 1:m
  %Probabilidad de que la kesima sustancia se desintegre
  P(:, a) = (lambda(:).*N(:,a).*deltat)./(norm(lambda(:).*N(:,a).*deltat));
  P2(:, a) = (lambda(:).*N2(:,a).*deltat)./(norm(lambda(:).*N2(:,a).*deltat));
  P3(:, a) = (lambda(:).*N3(:,a).*deltat)./(norm(lambda(:).*N3(:,a).*deltat));

  Sigmat = zeros(15, 15);
  Sigmat2 = zeros(15, 15);
  Sigmat3 = zeros(15, 15);

  %Calculo de la covarianza matriz
  for k = 1:15
    Sigmat = Sigmat + P(k, a)*DeltaN(:, k)*DeltaN(:, k)';
    Sigmat2 = Sigmat2 + P2(k, a)*DeltaN(:, k)*DeltaN(:, k)';
    Sigmat3 = Sigmat3 + P3(k, a)*DeltaN(:, k)*DeltaN(:, k)';
  endfor
  [X, sigmat] = lu(Sigmat);
  [X, sigmat2] = lu(Sigmat2);
  [X, sigmat3] = lu(Sigmat3);

  sigmat = sqrt(deltat)*sigmat;
  sigmat2 = sqrt(deltat)*sigmat2;
  sigmat3 = sqrt(deltat)*sigmat3;

  Y1(:,a) = N(:, a) + sigmat*(DeltaW(:,2)-DeltaW(:,1));
  Y2(:,a) = N2(:, a) + sigmat2*(DeltaW2(:,2)-DeltaW2(:,1));
  Y3(:,a) = N3(:, a) + sigmat3*(DeltaW3(:,2)-DeltaW3(:,1));

  %Probabilidad de que la kesima sustancia se desintegre
  PY(:, a) = (lambda(:).*Y1(:,a).*deltat)./(norm(lambda(:).*Y1(:,a).*deltat));
  PY2(:, a) = (lambda(:).*Y2(:,a).*deltat)./(norm(lambda(:).*Y2(:,a).*deltat));
  PY3(:, a) = (lambda(:).*Y3(:,a).*deltat)./(norm(lambda(:).*Y3(:,a).*deltat));

  SigmatY = zeros(15, 15);
  SigmatY2 = zeros(15, 15);
  SigmatY3 = zeros(15, 15);

  %Calculo de la covarianza matriz
  for k = 1:15
    SigmatY = SigmatY + PY(k, a)*DeltaN(:, k)*DeltaN(:, k)';
    SigmatY2 = SigmatY2 + PY2(k, a)*DeltaN(:, k)*DeltaN(:, k)';
    SigmatY3 = SigmatY3 + PY3(k, a)*DeltaN(:, k)*DeltaN(:, k)';
  endfor

  [X, sigmatY] = lu(SigmatY);
  [X, sigmatY2] = lu(SigmatY2);
  [X, sigmatY3] = lu(SigmatY3);

  sigmatY = sqrt(deltat)*sigmatY;
  sigmatY2 = sqrt(deltat)*sigmatY2;
  sigmatY3 = sqrt(deltat)*sigmatY3;

  %El ruido blanco
N(:,a+1)=pinv(eye(15)-deltat*M)*(N(:,a)+sigmat*(DeltaW(:,2)-DeltaW(:,1))+0.5*sqrt(deltat)*(sigmatY-sigmat)*(((DeltaW(:,2)-DeltaW(:,1))/sqrt(deltat)).^2-1));
  DeltaW(:, 1) = DeltaW(:, 2);
  DeltaW(:, 2) = sqrt(deltat)*rand(15,1);

N2(:,a+1)=pinv(eye(15)-deltat*M)*(N2(:,a)+sigmat2*(DeltaW2(:,2)-DeltaW2(:,1))+0.5*sqrt(deltat)*(sigmatY2-sigmat2)*(((DeltaW2(:,2)-DeltaW2(:,1))/sqrt(deltat)).^2-1));
  DeltaW2(:, 1) = DeltaW2(:, 2);
  DeltaW2(:, 2) = sqrt(deltat)*rand(15,1);

N3(:,a+1)=pinv(eye(15)-deltat*M)*(N3(:,a)+sigmat3*(DeltaW3(:,2)-DeltaW3(:,1))+0.5*sqrt(deltat)*(sigmatY3-sigmat3)*(((DeltaW3(:,2)-DeltaW3(:,1))/sqrt(deltat)).^2-1));
  DeltaW3(:, 1) = DeltaW3(:, 2);
  DeltaW3(:, 2) = sqrt(deltat)*rand(15,1);
endfor
%
%Figura 1
figure(1,"position",[0, 0, 1200, 1000])
subplot(2,1,1)
semilogx(N(1, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(1, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(1, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo1, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{1}: U235', 'N^{2}_{1}: U235', 'N^{3}_{1}: U235', 'Control')

%
subplot(2,1,2)
semilogx(N(2, :), 'Color', [1 0 0], 'LineWidth', 2)
hold on;
semilogx(N2(2, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(2, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo2, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{2}: Th231', 'N^{2}_{2}: Th231', 'N^{3}_{2}: Th231', 'Control')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 1 y 2.', 'fontsize', 14);

%Figura 2
figure(2,"position",[0, 0, 1200, 1000])
subplot(2,1,1)
semilogx(N(3, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(3, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(3, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo3, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{3}: Pa231', 'N^{2}_{3}: Pa231', 'N^{3}_{3}: Pa231', 'Control')

%
subplot(2,1,2)
semilogx(N(4, :), 'Color', [1 0 0], 'LineWidth', 2)
hold on;
semilogx(N2(4, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(4, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo4, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{4}: Ac227', 'N^{2}_{4}: Ac227', 'N^{3}_{4}: Ac227', 'Control')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 3y 4.', 'fontsize', 14);

%Figura 3
figure(3,"position",[0, 0, 1200, 1000])
subplot(2,1,1)
semilogx(N(5, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(5, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(5, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo5, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{5}: Th227', 'N^{2}_{5}: Th227', 'N^{3}_{5}: Th227', 'Control')

%
subplot(2,1,2)
semilogx(N(6, :), 'Color', [1 0 0], 'LineWidth', 2)
hold on;
semilogx(N2(6, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(6, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo6, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{6}: Fr223', 'N^{2}_{6}: Fr223', 'N^{3}_{6}: Fr223')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras5 y 6.', 'fontsize', 14);

%Figura 4
figure(4,"position",[0, 0, 1200, 1000])
subplot(2,1,1)
semilogx(N(7, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(7, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(7, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo7, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{7}: Ra223', 'N^{2}_{7}: Ra223', 'N^{3}_{7}: Ra223', 'Control')

%
subplot(2,1,2)
semilogx(N(8, :), 'Color', [1 0 0], 'LineWidth', 2)
hold on;
semilogx(N2(8, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(8, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo8, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{8}: Rn219', 'N^{2}_{8}: Rn219', 'N^{3}_{8}: Rn219', 'Control')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras7 y 8.', 'fontsize', 14);

%Figura 5
figure(5,"position",[0, 0, 1200, 1000])
subplot(2,1,1)
semilogx(N(9, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(9, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(9, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo9, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{9}: Po215', 'N^{2}_{9}: Po215', 'N^{3}_{9}: Po215', 'Control')

%
subplot(2,1,2)
semilogx(N(10, :), 'Color', [1 0 0], 'LineWidth', 2)
hold on;
semilogx(N2(10, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(10, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo10, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{10}: Pb211', 'N^{2}_{10}: Pb211', 'N^{3}_{10}: Pb211', 'Control')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras9 y 10.', 'fontsize', 14);

%Figura 6
figure(6,"position",[0, 0, 1200, 1000])
subplot(2,1,1)
semilogx(N(11, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(11, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(11, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo11, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{11}: At215', 'N^{2}_{11}: At215', 'N^{3}_{11}: At215', 'Control')

%
subplot(2,1,2)
semilogx(N(12, :), 'Color', [1 0 0], 'LineWidth', 2)
hold on;
semilogx(N2(12,:), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(12,:), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo12, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{12}: Bi211', 'N^{2}_{12}: Bi211', 'N^{3}_{12}: Bi211', 'Control')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 11 y 12.', 'fontsize', 14);

%Figura 7
figure(7,"position",[0, 0, 1200, 1000])
subplot(2,1,1)
semilogx(N(13, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(13, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(13, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo13, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{13}: Po211', 'N^{2}_{13}: Po211', 'N^{3}_{13}: Po211', 'Control')

%
subplot(2,1,2)
semilogx(N(14, :), 'Color', [1 0 0], 'LineWidth', 2)
hold on;
semilogx(N2(14, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(14, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo14, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{14}: Tl207', 'N^{2}_{14}: Tl207', 'N^{3}_{14}: Tl207', 'Control')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para las muestras 13 y 14.', 'fontsize', 14);

%Figura 8
figure(8,"position",[0, 0, 1200, 1000])
semilogx(N(15, :),'Color', [1 0 0],'LineWidth', 2);
hold on;
semilogx(N2(15, :), 'Color', [1 0.5 0.1], 'LineWidth', 2);
semilogx(N3(15, :), 'Color', [0.5 0.2 1], 'LineWidth', 2);
%semilogx(nucleo15, 'Color', [0 0 0], 'LineWidth', 2)
hold off;
xlabel('t/a', 'fontsize', 12)
ylabel('N(t)/núcleos', 'fontsize', 12)
legend('N^{1}_{15}: Pb207', 'N^{2}_{15}: Pb207', 'N^{3}_{15}: Pb207', 'Control')
figure1 = axes('visible', 'off', 'title', 'Gráficas del número de núcleos en el tiempo para lasmuestras 15.', 'fontsize', 14);
%
