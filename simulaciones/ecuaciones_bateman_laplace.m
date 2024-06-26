clc
clear

#Simulación numérica de la cadena de desintegración del uranio 235 mediante
#el método de Transformada de Laplace. Esta simulación contempla las fracciones de 
#desintegración en el proceso. 

#Definición del número de átomos en una muestra de 7.2331 kg de U235.
#atomosU235=(7.2331*10e+3)*(1/238.03)*(6.022e+23);

#El número de átomos en atomosU235 se obtiene mediante la obtención de la cantidad
#de sustancia presente en la muestra. El Uranio 235 tiene un peso molecular de 
#238.03 g/mol. El número 6.022e+23 corresponde al número de Avogrado de átomos de
#U235 en una mol de sustancia.

#Definición del intervalo de tiempo sobre la cual se simulará la desintegración.
#Se corre la variable tspan, intervalo de tiempo, desde 0 años hasta 9e+9 años 
#en pasos de 1000 años.  
#tspan=0:1000:9e+9;

#Definición del vectores de coeficientes de desintegración nuclear para cada
#elemento de la cadena incluyendo el isotopo estable PB207 al que se le atribuye
#una constante de decaimiento de 0 1/a. Suponiendo una período de semidesintegración
#infinito, dada su estabilidad nuclear.  
lambda = [9.846e-10, 1.786e-3, 2.116e-5, 3.184e-2, 1.353e+1, 1.656e+4, 2.213e+1, 5.520e+6, 1.227e+10, 1.009e+6, 2.186e+11, 1.702e+5, 4.236e+7, 7.638e+4, 0];

#Las fracciones de desintegración para los canales alfa y beta en las etapas de
#desintegración correspondientes al Actinio 227, Polonio 2015 y Bismuto 211.
ratio = [0.013800, 0.98620, 0.9999977, 0.0000023, 0.99724, 0.00276];

#Se definen los símbolos que constituirán las funciones de interés transformadas por Laplace
#así como el parámetro de tiempo y el parámetro de Laplace
syms t s Y01 Y02 Y03 Y04 Y05 Y06 Y07 Y08 Y09 Y10 Y11 Y12 Y13 Y14 Y15

#Los símbolos Y01 al Y15 quedan comprimidos en una matriz de símbolos
Y = [Y01, Y02, Y03, Y04, Y05, Y06, Y07, Y08, Y09, Y10, Y11, Y12, Y13, Y14, Y15];

#Definición del vector de condiciones iniciales del sistema de ecuaciones
y00 = [1.00; zeros(14,1)];

#Definición del lado derecho del sistema de ecuaciones "right hand side" 
RHS = zeros(15, 1);
%{
RHS = [0;
       lambda(1) * Y01;
       lambda(2) * Y02;
       lambda(3) * Y03;
       ratio(2) * lambda(4) * Y04;
       ratio(1) * lambda(4) * Y04;
       lambda(5) * Y05 + lambda(6) * Y06;
       lambda(7) * Y07;
       lambda(8) * Y08;
       ratio(3) * lambda(9) * Y09;
       ratio(4) * lambda(9) * Y09;
       lambda(10) * Y10 + lambda(11) * Y11;
       ratio(5) * lambda(12) * Y12;
       ratio(6) * lambda(12) * Y12;
       lambda(13) * Y13 + lambda(14) * Y14
       ];
%}

#Definición del lado izquierdo del sistema de ecuaciones conteniendo todas las operaciones
#y combinaciones de las variables transformadas bajo laplace
%
LHS = [s * Y01 - y00(1) + lambda(1) * Y01;
       s * Y02 - y00(2) + lambda(2) * Y02 - lambda(1) * Y01;
       s * Y03 - y00(3) + lambda(3) * Y03 - lambda(2) * Y02;
       s * Y04 - y00(4) + lambda(4) * Y04 - lambda(3) * Y03;
       s * Y05 - y00(5) + lambda(5) * Y05 - ratio(2) * lambda(4) * Y04;
       s * Y06 - y00(6) + lambda(6) * Y06 - ratio(1) * lambda(4) * Y04;
       s * Y07 - y00(7) + lambda(7) * Y07 - lambda(5) * Y05 - lambda(6) * Y06;
       s * Y08 - y00(8) + lambda(8) * Y08 - lambda(7) * Y07;
       s * Y09 - y00(9) + lambda(9) * Y09 - lambda(8) * Y08;
       s * Y10 - y00(10) + lambda(10) * Y10 - ratio(3) * lambda(9) * Y09;
       s * Y11 - y00(11) + lambda(11) * Y11 - ratio(4) * lambda(9) * Y09;
       s * Y12 - y00(12) + lambda(12) * Y12 - lambda(10) * Y10 - lambda(11) * Y11;
       s * Y13 - y00(13) + lambda(13) * Y13 - ratio(5) * lambda(12) * Y12;
       s * Y14 - y00(14) + lambda(14) * Y14 - ratio(6) * lambda(12) * Y12;
       s * Y15 - y00(15) + lambda(15) * Y15 - lambda(13) * Y13 - lambda(14) * Y14
      ];      
%
%{
LHS = [s * Y01 - y00(1) + lambda(1) * Y01;
       s * Y02 - y00(2) + lambda(2) * Y02;
       s * Y03 - y00(3) + lambda(3) * Y03;
       s * Y04 - y00(4) + lambda(4) * Y04;
       s * Y05 - y00(5) + lambda(5) * Y05;
       s * Y06 - y00(6) + lambda(6) * Y06;
       s * Y07 - y00(7) + lambda(7) * Y07;
       s * Y08 - y00(8) + lambda(8) * Y08;
       s * Y09 - y00(9) + lambda(9) * Y09;
       s * Y10 - y00(10) + lambda(10) * Y10;
       s * Y11 - y00(11) + lambda(11) * Y11;
       s * Y12 - y00(12) + lambda(12) * Y12;
       s * Y13 - y00(13) + lambda(13) * Y13;
       s * Y14 - y00(14) + lambda(14) * Y14;
       s * Y15 - y00(15) + lambda(15) * Y15
      ];
%}

#Solución del sistema de ecuaciones en el espacio de Laplace
for k = 1:15
  SOL(k) = solve(LHS(k) - RHS(k), Y(k));
endfor

for k = 1:15
  for j = 1:k-1
    SOL(k) = simplify(subs(SOL(k), Y(1:k-j), SOL(1:k-j)));
  endfor
endfor
SOL = transpose(SOL);

#Recuperación de las soluciones en el espacio de tiempo
[sol, polinomio_sol, SOL_fracparc] = my_ilaplace(SOL, lambda);

#Al mismo tiempo se hace la ecuperación de los coeficientes de Bateman
for k = 1:15
  coeficientes_bateman(k,1:length(coeffs(polinomio_sol(k), 'all'))) = double(coeffs(polinomio_sol(k), 'all'));
endfor
