function [sol, polinomio_sol, sol_espacio_laplace] = my_ilaplace (SOL, lambda)
  syms s t y
  for k=1:length(SOL)
    sol_espacio_laplace(k) = SOL(k);
    [N, D] = numden(sol_espacio_laplace(k));
    N = coeffs(N, 'all');  
    D = coeffs(D, 'all');
    [a, b, c, d] = residue(double(N), double(D));
%{
    if b(1) == 0
      b(1) = -lambda(1)
    endif
%}
    sol_espacio_laplace(k) = sum(a./(s-b).^sym(d));
    sol(k) = subs(ilaplace(sol_espacio_laplace(k), s, t), heaviside(t), 1);
    polinomio_sol(k) = subs(sol(k), exp(-t*lambda(1:k)), y.^[1:k]);
  endfor
endfunction
