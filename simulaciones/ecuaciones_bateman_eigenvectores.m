clc
clear

atomosU235=(7.2331*10e+3)*(1/238.03)*(6.022e+23);
t=0:100:9e9;
lambda=[9.7216e-10, 2.2793e+2, 2.0208e-5, 3.1507e-2, 1.3602e+1, 1.7348e+4, 2.2589e+1, 5.5763e+6, 1.1945e+10, 1.0092e+4, 1.6867e5, 4.2037e+7, 7.6058e+4, 0];

ratio=zeros(1,4);
ratio(2)=0.013800;
ratio(1)=1-ratio(2);
ratio(4)=0.99724;
ratio(3)=1-ratio(4);

M=zeros(14,14);
M(1,1)=-lambda(1); M(2,2)=-lambda(2); M(3,3)=-lambda(3); M(4,4)=-lambda(4);
M(5,5)=-lambda(5); M(6,6)=-lambda(6); M(7,7)=-lambda(7); M(8,8)=-lambda(8);
M(9,9)=-lambda(9); M(10,10)=-lambda(10); M(11,11)=-lambda(11); M(12,12)=-lambda(12);
M(13,13)=-lambda(13);

M(2,1)=lambda(1);M(3,2)=lambda(2); M(4,3)=lambda(3); M(5,4)=ratio(1)*lambda(4);
M(6,4)=ratio(2)*lambda(4); M(7,5)=lambda(5); M(7,6)=lambda(6);
M(8,7)=lambda(7); M(9,8)=lambda(8); M(10,9)=lambda(9); 
M(11,10)=lambda(10); M(12,11)=ratio(3)*lambda(11); M(13,11)=ratio(4)*lambda(11);
M(14,12)=lambda(12); M(14,13)=lambda(13);

[vectores,valores]=eig(M);

condicionesIniciales = zeros(14,1);
condicionesIniciales(1,1)=1;
Coeficientes=linsolve(vectores,condicionesIniciales);

%
vectores(:,1)=Coeficientes(1)*vectores(:,1);
vectores(:,2)=Coeficientes(2)*vectores(:,2);
vectores(:,3)=Coeficientes(3)*vectores(:,3);
vectores(:,4)=Coeficientes(4)*vectores(:,4);
vectores(:,5)=Coeficientes(5)*vectores(:,5);
vectores(:,6)=Coeficientes(6)*vectores(:,6);
vectores(:,7)=Coeficientes(7)*vectores(:,7);
vectores(:,8)=Coeficientes(8)*vectores(:,8);
vectores(:,9)=Coeficientes(9)*vectores(:,9);
vectores(:,10)=Coeficientes(10)*vectores(:,10);
vectores(:,11)=Coeficientes(11)*vectores(:,11);
vectores(:,12)=Coeficientes(12)*vectores(:,12);
vectores(:,13)=Coeficientes(13)*vectores(:,13);
vectores(:,14)=Coeficientes(14)*vectores(:,14);

%
exponenciales=zeros(length(t), 14);
exponenciales(:,1)=exp(valores(1,1)*t);
exponenciales(:,2)=exp(valores(2,2)*t);
exponenciales(:,3)=exp(valores(3,3)*t);
exponenciales(:,4)=exp(valores(4,4)*t);
exponenciales(:,5)=exp(valores(5,5)*t);
exponenciales(:,6)=exp(valores(6,6)*t);
exponenciales(:,7)=exp(valores(7,7)*t);
exponenciales(:,8)=exp(valores(8,8)*t);
exponenciales(:,9)=exp(valores(9,9)*t);
exponenciales(:,10)=exp(valores(10,10)*t);
exponenciales(:,11)=exp(valores(11,11)*t);
exponenciales(:,12)=exp(valores(12,12)*t);
exponenciales(:,13)=exp(valores(13,13)*t);
exponenciales(:,14)=exp(valores(14,14)*t);

nucleo1 = atomosU235*exponenciales*vectores(1,:)'; nucleo1(nucleo1<=0)=0;
nucleo2 = atomosU235*exponenciales*vectores(2,:)'; nucleo2(nucleo2<=0)=0;
nucleo3 = atomosU235*exponenciales*vectores(3,:)'; nucleo3(nucleo3<=0)=0;
nucleo4 = atomosU235*exponenciales*vectores(4,:)'; nucleo4(nucleo4<=0)=0;
nucleo5 = atomosU235*exponenciales*vectores(5,:)'; nucleo5(nucleo5<=0)=0;
nucleo6 = atomosU235*exponenciales*vectores(6,:)'; nucleo6(nucleo6<=0)=0;
nucleo7 = atomosU235*exponenciales*vectores(7,:)'; nucleo7(nucleo7<=0)=0;
nucleo8 = atomosU235*exponenciales*vectores(8,:)'; nucleo8(nucleo8<=0)=0;
nucleo9 = atomosU235*exponenciales*vectores(9,:)'; nucleo9(nucleo9<=0)=0;
nucleo10 = atomosU235*exponenciales*vectores(10,:)'; nucleo10(nucleo10<=0)=0;
nucleo11 = atomosU235*exponenciales*vectores(11,:)'; nucleo11(nucleo11<=0)=0;
nucleo12 = atomosU235*exponenciales*vectores(12,:)'; nucleo12(nucleo12<=0)=0;
nucleo13 = atomosU235*exponenciales*vectores(13,:)'; nucleo13(nucleo13<=0)=0;
nucleo14 = atomosU235*exponenciales*vectores(14,:)'; nucleo14(nucleo14<=0)=0;
%} 
 
%
figure(1,"position",[0, 0, 1800, 1500])
title('Gráficas del número de núcleos en el tiempo para las muestras 1 y 2.')
subplot(1,2,1)
semilogx(nucleo1, 'Color', [1 0 0], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{1}')

subplot(1,2,2)
semilogx(nucleo2, 'Color', [0 1 0], 'LineWidth',2)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{2}')

figure(2,"position",[0, 0, 1800, 1500])
title('Gráficas del número de núcleos en el tiempo para las muestras 3 y 4.')
subplot(1,2,1)
semilogx(nucleo3, 'Color', [0 0 1], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{3}')

subplot(1,2,2)
semilogx(nucleo4, 'Color', [0 1 1], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{4}')

figure(3,"position",[0, 0, 1800, 1500])
title('Gráficas del número de núcleos en el tiempo para las muestras 5 y 6.')
subplot(1,2,1)
semilogx(nucleo5, 'Color', [1 0 1], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{5}')

subplot(1,2,2)
semilogx(nucleo6, 'Color', [1 1 0], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{6}')

figure(4,"position",[0, 0, 1800, 1500])
title('Gráficas del número de núcleos en el tiempo para las muestras 7 y 8.')
subplot(1,2,1)
semilogx(nucleo7, 'Color', [0 0 0], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{7}')

subplot(1,2,2)
semilogx(nucleo8, 'Color', [0.6350 0.0780 0.1840], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{8}')

figure(5,"position",[0, 0, 1800, 1500])
title('Gráficas del número de núcleos en el tiempo para las muestras 9 y 10.')
subplot(1,2,1)
semilogx(nucleo9, 'Color', [0.3010 0.7450 0.9330], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{9}')

subplot(1,2,2)
semilogx(nucleo10, 'Color', [0.4660 0.6740 0.1880], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{10}')

figure(6,"position",[0, 0, 1800, 1500])
title('Gráficas del número de núcleos en el tiempo para las muestras 11 y 12.')
subplot(1,2,1)
semilogx(nucleo11, 'Color', [0.4940 0.1840 0.5560], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{11}')

subplot(1,2,2)
semilogx(nucleo12, 'Color', [0.9290 0.6940 0.1250], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{12}')

figure(7,"position",[0, 0, 1800, 1500])
title('Gráficas del número de núcleos en el tiempo para las muestras 13 y 14.')
subplot(1,2,1)
semilogx(nucleo13, 'Color', [0.8500 0.3250 0.0980], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{13}')

subplot(1,2,2)
semilogx(nucleo14, 'Color', [0 0.4470 0.7410], 'LineWidth',1.5)
xlabel('t/a')
ylabel('N(t)/núcleos')
title('N_{14}')

%}