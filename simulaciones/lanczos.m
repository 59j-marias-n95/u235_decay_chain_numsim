function [Q, T] = lanczos(A, k)
  % A es la matriz de entrada (simétrica)
  % k es el número de iteraciones (número de vectores propios a aproximar)
  % Q es la matriz de vectores ortogonales
  % T es la matriz tridiagonal resultante

  n = size(A, 1);    % Tamaño de la matriz A
  Q = zeros(n, k);   % Matriz de vectores ortogonales
  T = zeros(k, k);   % Matriz tridiagonal

  % Inicializar el primer vector aleatorio
  v = rand(n, 1);
  v = v / norm(v);   % Normalizar el primer vector

  % Primer valor propio
  beta = 0;

  for j = 1:k
    % Guardar el vector ortogonal
    Q(:, j) = v;

    % Multiplicar A por el vector v
    w = A * v;

    % Eliminar la componente en la dirección de Q(:, j-1)
    if j > 1
      w = w - beta * Q(:, j-1);
    end

    % Eliminar la componente en la dirección de Q(:, j)
    alpha = w' * Q(:, j);
    w = w - alpha * Q(:, j);

    % Calcular el valor propio (alpha) y el siguiente beta
    beta = norm(w);

    % Actualizar la matriz tridiagonal T
    T(j, j) = alpha;
    if j < k
      T(j, j+1) = beta;
      T(j+1, j) = beta;
    end

    % Normalizar el nuevo vector v
    if j < k
      v = w / beta;
    end
  end
end
