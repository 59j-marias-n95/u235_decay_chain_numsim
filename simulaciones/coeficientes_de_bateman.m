clc
clear

lambda=[9.7216e-10, 2.2793e2, 2.0208e-5, 3.1507e-2, 1.3602e1, 1.7348e4, 2.2589e1, 5.5763e6, 1.1945e10, 1.0092e4, 1.6867e5, 4.2037e7, 7.6058e4, 0];

n=14;
numerador=ones(n,1);
denominador=ones(n,n);

for k=1:n-1
  numerador(k+1,1)=numerador(k,1)*lambda(k); 
endfor

for k=3:n
  for j=1:n-k
    denominador(j+1,k)=denominador(j,k)*(lambda(j)-lambda(k));  
  endfor
endfor

for k=2:n
  j=1;
  while j<k
    denominador(k,k)=denominador(k-1,k)*(lambda(j)-lambda(k)); 
    j++;
  endwhile
endfor

for k=1:n
  for j=1:n-k
    denominador(j+k,k)=denominador(j+k-1,k)*(lambda(j+k)-lambda(k));  
  endfor
endfor

coeffbateman=numerador./denominador;