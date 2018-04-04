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
flt_op = 0;

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
fprintf("\n - x(1) = %f e x(2) = %f\n", x(1), x(n));
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
flt_op = 0;
[y, flt_op] = fGaussTRD(n, t, r, d, b);
fprintf("\n - x(1) = %f e x(2) = %f\n", y(1), y(n));
rmax2 = fresiduoTRD(n, t, r, d, b, y);
fprintf("\n - resíduo máximo = %f\n", rmax2);
fprintf("\n - Número de operações em ponto flutuante: %i\n", flt_op);

fprintf("\n-----------------------------------------------\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using iteractive method of Gauss-Seidel
%%%%

fprintf("\n - Método Iterativo de Gauss-Seidel - \n");
flt_op = 0;

fprintf(['\n - O sistema tem convergência garantida, pois possui uma diagonal\n', ... 
'dominante, ou seja, o elemento da diagonal principal é maior ou igual\n', ...
'à soma de seus elementos vizinhos e, em pelo menos uma linha, esse\n', ...
'elemento é maior que a soma dos seus vizinhos.\n']);

[z, flt_op] = fGaussSeidel(n, A, b);
fprintf("\n - x(1) = %f e x(2) = %f\n", z(1), z(n));
fprintf("\n - Número de operações em ponto flutuante: %i\n", flt_op);

fprintf("\n-----------------------------------------------\n");
