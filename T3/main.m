%---------------------------------------------------------------------
% Interpolação Polinomial
%intervalo
a = -2;
b = 2;

%grau
n = 16;

h = (b-a)/n;

x = a : h : b;

y = exp(-x.^2);

xpi = a : h/20 : b;

yp = fPnGregoryNewton(n, x, y, xpi);

ye = exp(-xpi.*xpi);

erroInter = abs(yp .- ye);

erroMaxInterpolacao = max(erroInter);

%plot(x,y,xp,yp);
%legend()
%xlim([-2 2.5])

fprintf("-Interpolação Polinomial Gregory-Newton-\n");
fprintf("Grau: %d \n", n);
fprintf("Erro Máximo: %e \n", erroMaxInterpolacao);
%fprintf("Coeficientes: \n");
%yp
fprintf("\n ------//------ \n");

%---------------------------------------------------------------------
% Interpolação Spline cúbica

n = 20;

h2 = (b-a)/n;

xspli = a : h2 : b;

yspli = exp(-xspli.^2);

% grau
n = length(xspli) - 1;
% número de pontos
m = length(xspli);

as=[];bs=[]; %zera memoria

[as bs c d]=fSplineb(m,xspli,yspli); % Splines quadráticas NAS PONTAS da figura S1=S2 e Sn+1=Sn

np=8; % 8 sub-divisões internas em cada n sub-intervalo entre x(i) e x(i+1)

xpp = []; ypp = [];

for i = 1 : n
    xs = xspli(i):(xspli(i+1)-xspli(i))/np:xspli(i+1);
    for k = 1 : np + 1
        ys(k) = as(i)*(xs(k)-xspli(i))*(xs(k)-xspli(i))*(xs(k)-xspli(i))+bs(i)*(xs(k)-xspli(i))*(xs(k)-xspli(i))+c(i)*(xs(k)-xspli(i))+d(i);
    end
    xpp = [xpp xs];
    ypp = [ypp ys];
end

yep = exp(-xpp.*xpp);

erroSpline = abs(ypp.-yep);
erroMaxSpline = max(erroSpline);

%plot(xspli,yspli,"*",xpp,ypp,'--b','LineWidth',1, xpp,yep,'-r','LineWidth',1)
%legend()
%xlim([-2 2])

fprintf("\n-Interpolação Spline Cúbica-\n");
fprintf("Grau: %d \n", n);
fprintf("Erro Máximo: %e \n", erroMaxSpline);
%fprintf("Coeficientes: \n");
%ypp
fprintf("\n ------//------ \n");

%---------------------------------------------------------------------
% Curvas de Bezier

x1 = [-2 -1.5 -1.3 -1]; y1 = [0.01831 0.055399 0.17452 0.3679];
x2 = [-1 -0.7 -0.3 0]; y2 = [0.3679 0.6170 0.9930 1];
x3 = [0 0.3 0.7 1]; y3 = [1 0.9930 0.6170 0.3679];
x4 = [1 1.3 1.5 2]; y4 = [0.3679 0.17452 0.055399 0.01831];

m=4; %numero de pontos:

[xx1, yy1] = fBezier(m,x1,y1);
[xx2, yy2] = fBezier(m,x2,y2);
[xx3, yy3] = fBezier(m,x3,y3);
[xx4, yy4] = fBezier(m,x4,y4);

xf = [xx1 xx2 xx3 xx4];
yf = [yy1 yy2 yy3 yy4];

%plot(xp,yp,xf,yf)
%plot(xp,yp,xx1,yy1)

ybez = exp(-xf.*xf);

erroMaxBezier = max(abs(yf.-ybez));

%plot(x,y,'x',x,y,'--k',x2,y2,'x',x2,y2,'--k','linewidth',1,xf,yf,'-k','linewidth',2)
%plot(x1,y1,'*', xf,yf,'--b')
%grid on
%xlim([-2.5 2.5])
%ylim([0 1])

fprintf("\n-Curvas de Bezier-\n");
fprintf("Erro Máximo: %e \n", erroMaxBezier);
%fprintf("Coeficientes: \n");
%yf
fprintf("\n ------//------ \n");

%---------------------------------------------------------------------
% Série de Maclaurin
tp = -1 : 0.1 : 1;

xp = 0.5.*(b-a).*tp .+ 0.5.*(b+a);

ye = exp(-4*tp.^2);

grau = 30;

c = fCoefMaclaurin(grau);

yMc = fPn(tp, grau, c);

%plot(xp, ye, '-b', xp, yMc, '-r');

% erro do método de Maclaurin
erroMaclaurin = abs(yMc .- ye);
erroMaxMc = max(erroMaclaurin);

%plot(xp, ye)

fprintf("\n-Série de Maclaurin-\n");
fprintf("Grau: %d \n", grau);
fprintf("Erro Máximo: %e \n", erroMaxMc);
fprintf("Coeficientes: \n");
c
fprintf("\n ------//------ \n");


%---------------------------------------------------------------------
% Série de Maclaurin-Tchebyschev

grau = 2;
coef = [0 0 4];
yMcTc = fPn(tp, grau, coef);

%plot(xp, ye, '-b', xp, yMc, '-r', xp, yMcTc, '-g');

erroMcTc = abs(yMcTc .- ye);
erroMaxMcTc = max(erroMcTc);

%plot(xp, erroMaclaurin, xp, erroMcTc, 2.2, 0)

fprintf("\n-Série de Maclaurin-Tchebychev-\n");
fprintf("Grau: %d \n", grau);
fprintf("Erro Máximo: %e \n", erroMaxMcTc);

fprintf("Encontrando os coeficientes:\n");
fprintf("\n1º Pegar a série de Maclaurin até grau 4\n")
fprintf("Mc4(t) = 1*t^0 + 0*t^1 -4*t^2 + 0*t^3 + 8*t^4\n")
fprintf("\n2º Trocar ti por Ti:\n")
fprintf("Mc4(T) = 1*T0 + 0*T1 -4*(T2 + T0 / 2) + 0*(T3 + 3T1 / 4) + 8(T4 + 4T2 + 3T0 / 8)\n")
fprintf("\n3º Agrupar os termos:\n")
fprintf("Tc3(T) = (1 - 2 + 3)T0 + (0 + 0)T1 + (-2 +4)T2 + (0)T3 + (1)\n")
fprintf("\n4º Eliminando o termo de maior grau (T4) e termos nulos:\n")
fprintf("Tc3(T) = 2T0 + 2T2\n")
fprintf("\n5º Voltando para t:\n")
fprintf("Tc3(t) = 2*t^0 + 2*(2*t^2) = 2 + 4*t^2 -2 = 4*t^2\n")

fprintf("\nCoeficientes: \n");
coef
fprintf("\n ------//------ \n");

%---------------------------------------------------------------------
% Tchebychev
grau = 12;
bt = fcoefTchebychev(grau);
yTc = fcalculaTchebychev(grau,bt,tp);

%plot(xp, ye, '-b', xp, yMc, '-r', xp, yMcTc, '-g', xp, yTc, '-m');

erroTc = abs(yTc .- ye);
erroMaxTc = max(erroTc);

%plot(xp, erroMaclaurin, xp, erroMcTc, xp, erroTc, 2.2, 0)

fprintf("\n-Tchebychev-\n");
fprintf("Grau: %d \n", grau);
fprintf("Erro Máximo: %e \n", erroMaxTc);
fprintf("Coeficientes: \n");
bt
fprintf("\n ------//------ \n");

%---------------------------------------------------------------------
% Série de Racional de Padé
% grau n e m (recomendado serem iguais ou próximos);
n = 9;
m = 8;

% M = n + m
M = n + m;

c = fCoefMaclaurin(M);

[a b] = fPade(n, m, c);

yPade = fPn(tp, n, a) ./ fPn(tp, m, b);

%erro Padé
erroPade = abs(yPade .- ye);
erroMaxPade = max(erroPade);

%plot(xp, ye, '-b', xp, yPade, '-r');

fprintf("\n-Série Racional de Padé-\n");
fprintf("Grau: n = %d , m = %d, M = %d\n", n, m, M);
fprintf("Erro Máximo: %e \n", erroMaxPade);
fprintf("Coeficientes: \n");
c
fprintf("\n ------//------ \n");

%---------------------------------------------------------------------
% Plots e Conclusoes
warning('off','all');

plot(xpi, erroInter, '--m', xpp, erroSpline, '--g', xp, erroMaclaurin, '--b', xp, erroTc, '--r', xp, erroPade, '--k');
xlim([-2 2.5])
legend({'Gregory Newton', 'Spline', 'Maclaurin', 'Tchebychev', 'Pade'},'NumColumns',2,'Location','North');


%fprintf("\nPodemos notar que a melhor aproximação foi a racional de Padé, por ter o menor grau para atingir o erro O(10^-4).\n");
fprintf("\nPodemos notar que Tchebychev foi a aproximação que atingiu o erro O(10^-4) com o menor valor de grau.\n");
