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
fprintf("\n - M�todo Elimina��o Gaussiana com Pivota��o Parcial - \n");
[x, flt_op] = fGauss(A,b);
fprintf("\n - x(1) = %.10f e x(n) = %.10f\n", x(1), x(n));
rmax = fresiduo(A, b, x);
fprintf("\n - res�duo m�ximo = %.10f\n", rmax);
fprintf("\n - N�mero de opera��es em ponto flutuante: %i\n", flt_op);

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
fprintf("\n - M�todo Otimizado Matriz Tridiagonal - \n");

[y, flt_op] = fGaussTRD(n, t, r, d, b);
fprintf("\n - x(1) = %.10f e x(n) = %.10f\n", y(1), y(n));
rmax2 = fresiduoTRD(n, t, r, d, b, y);
fprintf("\n - res�duo m�ximo = %.10f\n", rmax2);
fprintf("\n - N�mero de opera��es em ponto flutuante: %i\n", flt_op);

fprintf("\n-----------------------------------------------\n");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Using iteractive method of Gauss-Seidel
%%%%

fprintf("\n - M�todo Iterativo de Gauss-Seidel - \n");

fprintf(['\n - O sistema tem converg�ncia garantida, pois possui uma diagonal\n', ... 
'dominante, ou seja, todo elemento da diagonal principal � maior ou igual\n', ...
'� soma de seus elementos vizinhos e, em pelo menos uma linha, esse\n', ...
'elemento � maior que a soma dos seus vizinhos.\n']);

fprintf(['\nForam testados os seguintes fatores de relaxa��o:\n', ...
'lambda = 1     ->      k = 10\n', ...
'lambda = 0.8   ->      k = 15\n', ...
'lambda = 0.9   ->      k = 12\n', ...
'lambda = 1.1   ->      k = 9\n', ...
'lambda = 1.2   ->      k = 11\n', ...
'lambda = 1.05  ->      k = 8\n', ...
'lambda = 1.04  ->      k = 8\n', ...
'lambda = 1.06  ->      k = 8\n', ...
'Como lambda = 1.06 gerou o menor n�mero de itera��es, ele foi o escolhido.\n', ...
'Note que ao aumentar ou diminuir o valor de lambda, o n�mero de itera��es\n', ...
'apenas cresce.\n', ...
'O valor lambda = 1.06 tamb�m resultou em um menor n�mero de itera��es com\n', ...
'a toler�ncia^2, com k = 14, enquanto 1.05 e 1.07, geraram k = 15.\n']);

fprintf("\n- Valores com toler�ncia 1e-4 -\n");
tolerancia = 1e-4;
[z, flt_op] = fGaussSeidel(n, A, b, tolerancia);
fprintf("\n - x(1) = %.10f e x(n) = %.10f\n", z(1), z(n));
rmax3 = fresiduo(A, b, z);
fprintf("\n - res�duo m�ximo = %.10f\n", rmax3);
fprintf("\n - N�mero de opera��es em ponto flutuante: %i\n", flt_op);

fprintf("\n- Valores com toler�ncia (1e-4)� -\n");
[ze, flt_op] = fGaussSeidel(n, A, b, tolerancia^2);
fprintf("\n - x(1) = %.10f e x(n) = %.10f\n", ze(1), ze(n));
rmax4 = fresiduo(A, b, ze);
fprintf("\n - res�duo m�ximo = %.10f\n", rmax4);
fprintf("\n - N�mero de opera��es em ponto flutuante: %i\n", flt_op);

Errotruncamento=max(abs(z - ze));
fprintf("\n- Erro de truncamento: %.20f\n", Errotruncamento);

fprintf("\n-----------------------------------------------\n");
