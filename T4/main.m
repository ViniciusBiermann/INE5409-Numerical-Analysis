format long

%---------------------------------------------------------------------------
% Questao 7.1

m = 6;
T = [ 0.2 0.4 0.6 0.8 0.9 1.0 ];
V = [ 0.04 0.14 0.30 0.45 0.61 0.69 ];

% V(T) = ln(a + b*T^2)

ai = [1; 1];
a = fNewtonSisDerivadaNum(ai, T, V);

Tp = min(T) : 0.01 : max(T);
Vp = log(a(1) .+ a(2).*Tp.*Tp);

desvio = log(a(1) .+ a(2).*T.*T) .- V;
d2 = sum(desvio.*desvio);

n = m-1;

yp = fPnGregoryNewton(n, T, V, Tp);

%plot(T, V, '*', Tp, Vp, '-b', Tp, yp, '-g');xlim([0.2, 1.1]);

fprintf("\n-----------------------------------------------\n");
fprintf("\n - Questão 7.1 - \n");
fprintf("\n equação 1: sum( (log(a + b * T^2 ) - V ) * 1/(a + b * T^2 )) \n")
fprintf("\n equação 2: sum( (log(a + b * T^2 ) - V ) * T^2/(a + b * T^2 )) \n")
fprintf("\n Coeficientes: a = %.15f e b = %.15f \n", a(1), a(2))
fprintf("\n Desvios: \n")
desvio
fprintf("\n Desvio quadrático: %.15f \n", d2)
fprintf("\n-----------------------------------------------\n");

%---------------------------------------------------------------------------
% Questao 7.6

m = 7;
T = [ 13.9 37.0 67.8 79.0 85.5 93.1 99.2 ];
V = [ 1.04 1.18 1.29 1.35 1.28 1.21 1.06 ];

%---------------------------------------------------------------------------
% Questao 7.7

m = 4;
T = [ 0.00 0.39 0.78 1.18 ];
V = [ 0.99 0.92 0.71 0.38 ];
% V = a * T + b * cos(T), T em radianos

a = fCalcCoef(m, T, V);

Tp = min(T) : 0.01 : max(T);

Vp = a(1).*Tp .+ a(2).*cos(Tp);

d2 = a(1).*T .+ a(2).*cos(T) .- V;
D2 = sum(d2.*d2);

%plot(T, V, '*', Tp, Vp, '-b');xlim([0, 1.4]);
%plot(T, d2, '*')

fprintf("\n-----------------------------------------------\n");
fprintf("\n - Questão 7.7 - \n");

fprintf("\n Desvios: \n")
d2

fprintf(['Podemos notar que todos os desvios estão na ordem de 10^-3 \n', ...
         'e o gráfico da função teve um comportamento bom, assim, a função \n', ...
         'fornecida é adequada à amostra. \n']);
         
fprintf("\n-----------------------------------------------\n");

%---------------------------------------------------------------------------
% Questao 8.6

fprintf("\n-----------------------------------------------\n");
fprintf("\n - Questão 8.6 - \n");

a = 0;
b = 1;

% ////////////////
% Metodo dos trapezios
n = 130;

Tn = fTn(n, a, b);

% ErroTn = |Tn - Ie|
Ie = erf(b) - erf(a);

% Erro exato
ErroTn = abs(Tn - Ie);

% Erro aproximado(estimado)
Tn2 = fTn(2*n, a, b);
ErroTnApprox = abs(Tn - Tn2);

fprintf("\n - Método dos trapézios - \n");
fprintf("\n n = %i \n", n)
fprintf("\n Tn = %e \n", Tn)
fprintf("\n Erro Exato: %e \n", ErroTn)
fprintf("\n Erro Estimado: %e \n", ErroTnApprox)

% ////////////////
% Metodo de Simpson
n = 8;

Sn = fSn(n, a, b);

% Erro exato Sn
ErroSn = abs(Sn - Ie);

% Erro aproximado(estimado) Sn
Sn2 = fSn(2*n, a, b);
ErroSnApprox = abs(Sn - Sn2);

fprintf("\n - Método de Simpson - \n");
fprintf("\n n = %i \n", n)
fprintf("\n Sn = %e \n", Sn)
fprintf("\n Erro Exato: %e \n", ErroSn)
fprintf("\n Erro Estimado: %e \n", ErroSnApprox)

% ////////////////
% Metodo de Gauss-Legendre
m = 4;

f1 = @(x) 2/(sqrt(pi)).*exp(-x.^2);
Gm = fGm(m, a, b, f1);

% erro exato Gm
ErroGm = abs(Gm - Ie);

% Erro aproximado(estimado) Gm
% m + 1 já é suficiente para capturar o erro
Gm2 = fGm(m+1, a, b, f1);
ErroGmApprox = abs(Gm - Gm2);

fprintf("\n - Método de Gauss-Legendre - \n");
fprintf("\n m = %i \n", m)
fprintf("\n Gm = %e \n", Gm)
fprintf("\n Erro Exato: %e \n", ErroGm)
fprintf("\n Erro Estimado: %e \n", ErroGmApprox)


fprintf("\n-----------------------------------------------\n");

%---------------------------------------------------------------------------
% Questao 8.7

a = -1;
b = 1;

m1 = 10

fG = @(x) (log(1.+x))./sqrt(1.-x.*x);
Gm = fGm(m1, a, b, fG);

Ie = -pi*log(2)

ErroGm = abs(Gm - Ie);

% ////////////////
% Utilizando Gauss-Tchebyshev
m2  = 1000 %10;
fT = @(x) log(1.+x);
GTm = fGTm(m2, fT);

erroGTm = abs(GTm - Ie);

fprintf("\n-----------------------------------------------\n");
fprintf("\n - Questão 8.7 - \n");
fprintf("\n Via Gauss-Legendre \n")
fprintf("\n m = %i \n", m1)
fprintf("\n Gm = %.15f \n", Gm)
fprintf("\n erroGm = %.15f \n", ErroGm)

fprintf("\n Via Gauss-Tchebyshev \n")
fprintf("\n m = %i \n", m2)
fprintf("\n GTm = %.15f \n", GTm)
fprintf("\n erroGm = %.15f \n", erroGTm)

fprintf(['\nAmbos os métodos não se aplicam bem para a função dada,  devido à \n', ...
         'descontinuidade da função em x = -1.\n \n'])

