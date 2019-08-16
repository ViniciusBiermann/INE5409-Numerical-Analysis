
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%

%{
    f(x) = x * tg(x) - 1
%}

A = 0; B = 5*pi;
x = A : 0.1 : B;
y = f(x);

xi = fLocalizacao(A, B);

for i = 1 : 5
    roots1(i) = fNewton(xi(i));
end


fprintf("\n-----------------------------------------------\n");
fprintf('Zeros da função f(x) = x*tg(x) - 1:\n')
funcZeros = transpose(roots1)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
%%%%

a=[1 -3.4 +2.35 4.32 -7.1685 +1.56006 +3.287061 -2.480058 0.531441];

% Segundo parametro é para definição da multiplicidade. Se for 1, a multiplicidade
% é 1, se for qualquer outro número, a multiplicidade é estimada.
roots2 = raizes(a, 1);

roots3 = raizes(a, 2);

roots4 = roots(a);

fprintf("\n-----------------------------------------------\n");
fprintf(['Zeros do polinômio ', polyout(a, 'x'), ' :\n'])
fprintf('\nUtilizando método de Newton tradicional. Com multiplicidade 1:\n')
M1Roots = roots2

fprintf('\nUtilizando método de Newton com estimativa de multiplicidade:\n')
EstMRoots = roots3

fprintf(['\nPolinômio fatorado: ((x - 0.9)^6)*((x + 1)^2)\n'])

fprintf('\nRaízes obtidas pela função roots:\n')
OctaveRoots = roots4

fprintf('\nRaízes obtidas pelo Wolfram|Alpha:\n')
WolfRoots = [ '-1.000001 - 0.000922i'; '-1.000001 + 0.000922i';
              ' 0.789954 - 0.053450i'; ' 0.789954 + 0.053450i';
              ' 0.886049 - 0.136679i'; ' 0.886049 + 0.136679i';
              ' 1.023998 - 0.078400i'; ' 1.023998 + 0.078400i']

