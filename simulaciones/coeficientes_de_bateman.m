%clc
%clear

lambda=[9.846e-10, 1.786e-3, 2.116e-5, 3.184e-2, 1.353e+1, 1.656e+4, 2.213e+1, 5.520e+6, 1.227e+10, 1.009e+4, 2.186e+11, 1.702e+5, 4.236e+7, 7.638e+4, 0];

lambdacanales = [0.4394e-3, 3.140e-3, 1.227e10, 2.822e4, 2.180e11, 6.033e8];

n=15;
numerador=zeros(n,n);
numerador(1,1)=1;
denominador=ones(n,n);

%{
for k=1:n-1
    if k+1==5
      numerador(k+1,1)=numerador(k,1)*lambdacanales(2);
    elseif k+1==6
      numerador(k+1,1)=numerador(k,1)*lambdacanales(1);
    elseif k+1==10
      numerador(k+1,1)=numerador(k,1)*lambdacanales(3);
    elseif k+1==11
      numerador(k+1,1)=numerador(k,1)*lambdacanales(4);
    elseif k+1==13
      numerador(k+1,1)=numerador(k,1)*lambdacanales(6);
    elseif k+1==14
      numerador(k+1,1)=numerador(k,1)*lambdacanales(5);
    else
      numerador(k+1,1)=numerador(k,1)*lambda(k);
    endif
endfor
%}

for k=1:n-1
    numerador(k+1,1)=numerador(k,1)*lambda(k);
endfor

for j=2:n
  for k=j:n
    numerador(k,j)=numerador(k,j-1);
  endfor
endfor

%{
for k=1:n
  for j=1:n
    if j<k
      if j+1==5
        denominador(j+1,k)=denominador(j,k)*(lambdacanales(2)-lambda(k));
      elseif j+1==6
        denominador(j+1,k)=denominador(j,k)*(lambdacanales(1)-lambda(k));
      elseif j+1==10
        denominador(j+1,k)=denominador(j,k)*(lambdacanales(3)-lambda(k));
      elseif j+1==11
        denominador(j+1,k)=denominador(j,k)*(lambdacanales(4)-lambda(k));
      elseif j+1==13
        denominador(j+1,k)=denominador(j,k)*(lambdacanales(6)-lambda(k));
      elseif j+1==14
        denominador(j+1,k)=denominador(j,k)*(lambdacanales(5)-lambda(k));
      else
        denominador(j+1,k)=denominador(j,k)*(lambda(j)-lambda(k));
        continue
      endif
    elseif j==k 
      continue 
    else
      if k==5
        denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambdacanales(2));
      elseif k==6
        denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambdacanales(1));
      elseif k==10
        denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambdacanales(3));
      elseif k==11
        denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambdacanales(4));
      elseif k==13
        denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambdacanales(6));
      elseif k==14
        denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambdacanales(5));
      else
        denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambda(k));
      endif      
    endif
  endfor
endfor
%}

for k=1:n
  for j=1:n
    if j<k
       denominador(j+1,k)=denominador(j,k)*(lambda(j)-lambda(k));
       continue
    elseif j==k 
      continue 
    else
      denominador(j,k)=denominador(j-1,k)*(lambda(j)-lambda(k));      
    endif
  endfor
endfor

coeffbateman=fliplr(numerador./denominador);
