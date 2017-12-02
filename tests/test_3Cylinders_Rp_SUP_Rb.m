clear; close all; clc;


h = 150;
rho = 1025;
g = 9.81;

WECstructure.Rbo = 15;          % diamètre exterieur - dDe
WECstructure.Rbi = 10;          % diamètre interieur - dDi
WECstructure.Rc = 9.95;%8.3         % diamètre de la colonne - dCeh  --> 8
WECstructure.b = 10;           % tirant d'eau flotteur
WECstructure.e1 = 25.5;          % tirant d'eau embase-colonne
WECstructure.e2 = WECstructure.e1+10;        % hauteur de l'embase heB

WECstructure.Rp=20;%

% Numerical truncation 
options.Truncate.Ni = 55;
options.Truncate.Nl = 150;
options.Truncate.Nn = 150;
options.Truncate.Nj = 150;

options.Zc = [0,0]';

Omega=(.1:.1:3)'; nw=length(Omega);

% 
[Fe1, A1, B1, A_inf1, Fe_Haskind1] = threeCylinders_Rp_SUP_Rbo__Rbi_EQ_Rc(Omega, h, WECstructure, options);
[Fe2, A2, B2, A_inf2, Fe_Haskind2] = threeCylinders_Rp_SUP_Rbo__Rbi_SUP_Rc(Omega, h, WECstructure, options);

i = 2; j = 2;

figure;
grid('on'), hold('on');
p1 = plot(Omega, [abs(Fe1(:, i)), abs(Fe2(:, i))]);
p2 = plot(Omega, [abs(Fe_Haskind1(:, i)), abs(Fe_Haskind2(:, i))], '+');
p2(1).Color = p1(1).Color;
p2(2).Color = p1(2).Color;
xlabel('$$\omega$$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('Excitation')

figure;
grid('on'), hold('on');
p1 = plot(Omega, [reshape(A1(i, j, :), [nw, 1]), reshape(A2(i, j, :), [nw, 1])]);
p2 = plot(Omega, ones(nw, 1)*[A_inf1(i, j), A_inf2(i, j)], '-.');
p2(1).Color = p1(1).Color;
p2(2).Color = p1(2).Color;
xlabel('$$\omega$$ [rad/s]', 'Interpreter', 'latex', 'FontSize', 12)
ylabel('added mass [kg]')










