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
m = 1.00e+04;
deltat = 2.0e+02;

%El numero inicial de particulas: (7.2331e+03)*(1/238.03)*(6.022e+23)
%atomosu235 = (7.2331e+03)*(1/238.03)*(6.022e+23);
N0 = 1e+15;
N = [[N0, zeros(1,14)]', zeros(15, 1)];
N2 = [[N0, zeros(1,14)]', zeros(15, 1)];
N3 = [[N0, zeros(1,14)]', zeros(15, 1)];
DeltaW = zeros(15, 2);
DeltaW2 = zeros(15, 2);
DeltaW3 = zeros(15, 2);

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
    Sigmat = Sigmat + P(k, a)*DeltaN(1:15, k)*DeltaN(1:15, k)';
    Sigmat2 = Sigmat2 + P2(k, a)*DeltaN(1:15, k)*DeltaN(1:15, k)';
    Sigmat3 = Sigmat3 + P3(k, a)*DeltaN(1:15, k)*DeltaN(1:15, k)';
  endfor
  [X, sigmat] = lu(Sigmat);
  sigmat = sqrt(deltat)*sigmat;

  [X, sigmat2] = lu(Sigmat2);
  sigmat2 = sqrt(deltat)*sigmat2;

  [X, sigmat3] = lu(Sigmat3);
  sigmat3 = sqrt(deltat)*sigmat3;

  %El ruido blanco
  N(:, a + 1) = pinv(eye(15) - deltat*M)*(N(1:15, a) + sigmat*(DeltaW(1:15, 2) - DeltaW(1:15, 1)));
  DeltaW(:, 1) = DeltaW(1:15, 2);
  DeltaW(:, 2) = sqrt(deltat)*(rand(15,1));

  N2(:, a + 1) = pinv(eye(15) - deltat*M)*(N2(1:15, a) + sigmat2*(DeltaW2(1:15, 2) - DeltaW2(1:15, 1)));
  DeltaW2(:, 1) = DeltaW2(1:15, 2);
  DeltaW2(:, 2) = sqrt(deltat)*(rand(15,1));

  N3(:, a + 1) = pinv(eye(15) - deltat*M)*(N3(1:15, a) + sigmat3*(DeltaW3(1:15, 2) - DeltaW3(1:15, 1)));
  DeltaW3(:, 1) = DeltaW3(1:15, 2);
  DeltaW3(:, 2) = sqrt(deltat)*(rand(15,1));
endfor

for a = 1:m
  EN(:, a) = (N(:,a) + N2(:,a) + N3(:,a))./(3);
endfor

#Definición del intervalo de tiempo sobre la cual se simulará la desintegración.
#Se corre la variable t, tiempo, desde 0 años hasta 9e+9 años en pasos de 1000
#años.
t=0:deltat:m*deltat;

#El cálculo de los eigenvectores y eigenvalores correspondientes para la matriz
#del sistema surgen de la siguiente evaluación numérica.
[vectores,valores]=eig(M);

#La solución del sistema se compone de una combinación lineal de los vectores
#propios del sistema, por lo que se desea encontrar los coeficientes que hacen
#que la solución satisfaga a las condiciones iniciales
condicionesIniciales = [N0; zeros(14,1)];
Coeficientes=linsolve(vectores,condicionesIniciales);

for i = 1:15
  Vectores(:,i)=Coeficientes(i)*vectores(:,i);
endfor

for i = 1:15
  exponenciales(:,i)=exp(valores(i,i)*t);
endfor

for k = 1:15
  nucleos(k, :) = exponenciales*Vectores(k, :)';
  Errores(k) = norm(nucleos(k,1:m) - EN(k,1:m))/norm(nucleos(k,1:m))*100;
endfor
Errores'

%{
figure(1)
semilogx(EN(1, :))
hold on
semilogx(nucleos(1, :))
hold off

figure(2)
semilogx(EN(2, :))
hold on
semilogx(nucleos(2, :))
hold off

figure(3)
semilogx(EN(3, :))
hold on
semilogx(nucleos(3, :))
hold off

figure(4)
semilogx(EN(4, :))
hold on
semilogx(nucleos(4, :))
hold off

figure(5)
semilogx(EN(5, :))
hold on
semilogx(nucleos(5, :))
hold off

figure(6)
semilogx(EN(6, :))
hold on
semilogx(nucleos(6, :))
hold off

figure(7)
semilogx(EN(7, :))
hold on
semilogx(nucleos(7, :))
hold off

figure(8)
semilogx(EN(8, :))
hold on
semilogx(nucleos(8, :))
hold off

figure(9)
semilogx(EN(9, :))
hold on
semilogx(nucleos(9, :))
hold off

figure(10)
semilogx(EN(10, :))
hold on
semilogx(nucleos(10, :))
hold off

figure(11)
semilogx(EN(11, :))
hold on
semilogx(nucleos(11, :))
hold off

figure(12)
semilogx(EN(12, :))
hold on
semilogx(nucleos(12, :))
hold off

figure(13)
semilogx(EN(13, :))
hold on
semilogx(nucleos(13, :))
hold off

figure(14)
semilogx(EN(14, :))
hold on
semilogx(nucleos(14, :))
hold off

figure(15)
semilogx(EN(15, :))
hold on
semilogx(nucleos(15, :))
hold off}%
