format long
% System of linear equations with n = 50.
%{
    for i = 1
        x(i) + x(i+1) = 1.50
    for i = 2 : n/2
        x(i-1) + 4x(i) + x(i+1) = 1.00
    for i = (n/2) + 1 : n - 1
        x(i-1) + 5x(i) + x(i+1) = 2.00
    for i = n
        x(i-1) + x(i) = 3.00
%}

n = 50;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using direct method of Gaussian Elimination with partial pivoting
%%%%
A(1,1) = 1.00;
A(1,2) = 1.00;
b(1,1) = 1.50;

for i = 2 : n/2
    A(i,i-1) = 1.00;
    A(i,i) = 4.00;
    A(i,i+1) = 1.00;
    b(i,1) = 1.00;
end

for i = (n/2) + 1 : n - 1
    A(i,i-1) = 1.00;
    A(i,i) = 5.00;
    A(i,i+1) = 1.00;
    b(i,1) = 2.00;
end

A(n,n-1) = 1.00;
A(n,n) = 1.00;
b(n,1) = 3.00;

%%%%
fprintf("\n-----------------------------------------------\n");
fprintf("\n - Método Eliminação Gaussiana com Pivotação Parcial - \n");
[x, flt_op] = fGauss(A,b);
fprintf("\n - x(1) = %f e x(n) = %f\n", x(1), x(n));
rmax = fresiduo(A, b, x);
fprintf("\n - resíduo máximo = %f\n", rmax);
fprintf("\n - Número de operações em ponto flutuante: %i\n", flt_op);

fprintf("\n-----------------------------------------------\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using optimized method for tridiagonal matrix
%%%%

t(1) = 0.00;
r(1) = 1.00;
d(1) = 1.00;

for i = 2 : n/2
    t(i) = 1.00;
    r(i) = 4.00;
    d(i) = 1.00;
end

for i = (n/2) + 1 : n - 1
    t(i) = 1.00;
    r(i) = 5.00;
    d(i) = 1.00;
end

t(n) = 1.00;
r(n) = 1.00;
d(n) = 0.00;

%%%%
fprintf("\n - Método Otimizado Matriz Tridiagonal - \n");

[y, flt_op] = fGaussTRD(n, t, r, d, b);
fprintf("\n - x(1) = %f e x(n) = %f\n", y(1), y(n));
rmax2 = fresiduoTRD(n, t, r, d, b, y);
fprintf("\n - resíduo máximo = %f\n", rmax2);
fprintf("\n - Número de operações em ponto flutuante: %i\n", flt_op);

fprintf("\n-----------------------------------------------\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using iteractive method of Gauss-Seidel
%%%%

fprintf("\n - Método Iterativo de Gauss-Seidel - \n");

fprintf(['\n - O sistema tem convergência garantida, pois possui uma diagonal\n', ... 
'dominante, ou seja, o elemento da diagonal principal é maior ou igual\n', ...
'à soma de seus elementos vizinhos e, em pelo menos uma linha, esse\n', ...
'elemento é maior que a soma dos seus vizinhos.\n']);

fprintf(['\nForam testados os seguintes fatores de relaxação:\n', ...
'lambda = 1     ->      k = 10\n', ...
'lambda = 0.8   ->      k = 15\n', ...
'lambda = 0.9   ->      k = 12\n', ...
'lambda = 1.1   ->      k = 9\n', ...
'lambda = 1.2   ->      k = 11\n', ...
'Como lambda = 1.1 gerou o menor número de iterações, ele foi o escolhido.\n']);

fprintf("\n- Valores com tolerância 1e-4 -\n");
tolerancia = 1e-4;
[z, flt_op] = fGaussSeidel(n, A, b, tolerancia);
fprintf("\n - x(1) = %f e x(n) = %f\n", z(1), z(n));
fprintf("\n - Número de operações em ponto flutuante: %i\n", flt_op);

fprintf("\n- Valores com tolerância (1e-4)^2 -\n");
[ze, flt_op] = fGaussSeidel(n, A, b, tolerancia^2);
fprintf("\n - x(1) = %f e x(n) = %f\n", ze(1), ze(n));
fprintf("\n - Número de operações em ponto flutuante: %i\n", flt_op);

Errotruncamento=max(abs(z-ze));
fprintf("\n- Erro de truncamento: %.20f\n", Errotruncamento);

fprintf("\n-----------------------------------------------\n");
