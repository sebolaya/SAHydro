function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = oneCylinder_SPAR(Omega, depth, WECstructure, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the hydrodynamic coefficients for the two coaxial 
% cylinder structure depicted below.
% 
% Inputs :
% - Omega				: Vector of wave frequencies (rad/s)
% - depth				: Water depth (m), 0 for deep water 
% - WECstructure		: WEC description (see below)
% - options				: Some options
%		* Trunc : Numerical truncation of the infinite series (see Ref)
%		* Zc    : Center of rotation for pitch moment evaluation 
%				  Use CoG for exemple, default value is Zc=[0;0];
%		* Dir   : DATA subdirectory for hydrodynamic coefficients saving
%
% Outputs :
% - Fe : Matrix length(Omega)x3x2 of exciation forces (complex values)
% - A  : Matrix (3x2)x(3x2)xlength(Omega) of added mass coefficients
% - B  : Matrix (3x2)x(3x2)xlength(Omega) of radiation damping coefficients
% - A_inf : Matrix (3x2)x(3x2)xlength(Omega) of infinite added mass
% - Fe_Haskind : Fe computation with the Haskind's theorem 
%				 (for verification purposes)
% - Fe_FK : Froude-Krilov efforts -> Fe without diffraction phenomena
%
% *************************************************************************
%					WEC Description
%
%
%           |-----> (Rc)
%           .                                                              
% --------  +-----+  -------------------------------- (z=0)
%           .ooooo|				^      ^        ^
%           |ooooo|				|      |        |
%           .ooooo|				|      |        |
%           |ooooo|				|      |        |
%           .ooooo|				|      |        |
%           |ooooo|				|      |        |
%           .ooooo|				|      |        |
%           |ooooo|				|      |        |
%           .ooooo|				|      |        |
%           |ooooo|				|      |        |
%           .ooooo|				|      |        |
%           |ooooo|				|      |        |
%           .ooooo|				|      |        |
%           |ooooo|				|      |        |
%           .ooooo|				|      |        |
%           |ooooo|				|      |        |
%           .ooooo|				v (e1) |        |
%           |ooooo+------------+	   |        |
%           .ooooo|oooooooooooo|       |        |
%           |ooooo|oooooooooooo|       |        |
%           +-----+------------+       v (e2)   .
%           .                                   .
%           |------------------> (Rp)           .
%                                               |
%                                               |
%                                               |
%                                               v (h)
%
%
%
% Written by S. OLAYA, IRDL/ENIB*, olaya@enib.fr
% * ENIB - École Nationale d'Ingénieurs de Brest, France
% Date : 2016/15/03
% Please include the following citation when you use this code :
%
% Olaya, S., Bourgeot, J.-M., & Benbouzid, M. E.-H. (2015). Hydrodynamic 
% Coefficient Computation for a Partially Submerged Wave Energy Converter.
% IEEE Journal of Oceanic Engineering, 40(3), 522–535. 
% http://doi.org/10.1109/JOE.2014.2344951
%
% Revisions :
%
% Copyright (C) 2016 Sébastien OLAYA
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the 
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if nargin<4
	options=[];
end

% DATA sub-directory
if ~isfield(options,'DIR')
	DIR=pwd;
else
	DIR = options.DIR;
end

if ~isdir([DIR,filesep,'DATA'])
	mkdir([DIR,filesep,'DATA']);
	addpath(genpath(DIR));
	addpath(genpath([DIR,filesep,'DATA']));
end

g = 9.81; % gravity
rho = 1025; % water density
h = depth; % water depth

Rc = WECstructure.Rc;
Rp = WECstructure.Rp;

e1 = WECstructure.e1;
e2 = WECstructure.e2;

if ~isfield(options,'Trunc')
	Ni=80; Nn=100; Nj=100;
else
	Ni = options.Trunc.Ni;
	Nn = options.Trunc.Nn;
	Nj = options.Trunc.Nj;
end


if ~isfield(options,'Zc')
	Zc=0;
else
	Zc=options.Zc;
end

epsilon = 1e-7;

%%% Initialisation des sorties
Zr = zeros(3,3,length(Omega)); %%% matrice d'impédance --> A + (1i/omega)*B
A = zeros(3,3,length(Omega));
B = zeros(3,3,length(Omega));
Fe = zeros(length(Omega),3); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),3); %%% effort d'excitation
Fe_FK = zeros(length(Omega),3);%%% Froude-Krilov

j = (1:1:Nj-1)';
lambda_j = (pi/(h-e2)).*j;

for w=1:length(Omega)+1

clc;
disp([num2str(round(w*100/(length(Omega)+1))),'%']);

if w<=length(Omega)
	omega = Omega(w);
	lambda_i=fcn_waveDispersion(Omega(w),h,Ni);
	lambda_n=fcn_waveDispersion(Omega(w),e1,Nn); lambda_n=conj(lambda_n)';

	k0 = -imag(lambda_i(1));
else
	omega=Inf;
	i=1:1:Ni;
	n=(1:1:Nn)';

	lambda_i = .5*(2*i-1)*pi/h;
	lambda_n = .5*(2*n-1)*pi/e1;
end

N_lambda_i = .5.*(1 + sin(2.*lambda_i.*h)./(2.*lambda_i.*h));
N_lamnda_n = .5.*(1 + sin(2.*lambda_n.*e1)./(2.*lambda_n.*e1));

L_jtau = zeros(Nj,Ni);
L_jtau(1,:) = N_lambda_i.^(-.5).*sin(lambda_i.*(h-e2))./( (h-e2).*lambda_i );
for i = 1 : Ni
	for p = 1 : Nj-1
		if (lambda_j(p)/lambda_i(i) < 1+epsilon)&&(lambda_j(p)/lambda_i(i) > 1-epsilon)
			L_jtau(p+1,i) = .5*sqrt(2).*N_lambda_i(i)^(-.5);
		else
			L_jtau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*(sin((lambda_j(p)-lambda_i(i))*(h-e2))/(lambda_j(p)-lambda_i(i))...
				+ sin((lambda_j(p)+lambda_i(i))*(h-e2))/(lambda_j(p)+lambda_i(i)))/(h-e2);
		end
	end
end

K_ltau = zeros(Nn,Ni);
for i = 1 : Ni
	for p = 1 : Nn
		if (lambda_n(p)/lambda_i(i) < 1+epsilon)&&(lambda_n(p)/lambda_i(i) > 1-epsilon)
			K_ltau(p,i) = (1/e1)*(N_lambda_i(i)^(-.5)*N_lamnda_n(p).^(-.5))*...
				(.25*(sin(lambda_i(i)*(h+e1)) - sin(lambda_i(i)*(h-e1)))/lambda_i(i) + e1*.5*cos(lambda_i(i)*(h-e1)));
		else
			K_ltau(p,i) = (1/e1)*(N_lambda_i(i)^(-.5)*N_lamnda_n(p).^(-.5))*...
				.5*( (sin(lambda_n(p)*e1-lambda_i(i)*h)-sin(lambda_n(p)*e1-lambda_i(i)*h - e1*(lambda_n(p)-lambda_i(i))))/(lambda_n(p)-lambda_i(i))...
						+ (sin(lambda_n(p)*e1+lambda_i(i)*h)-sin(lambda_n(p)*e1+lambda_i(i)*h - e1*(lambda_n(p)+lambda_i(i))))/(lambda_n(p)+lambda_i(i)));
		end
	end
end

Gamma_0j(1,1) = 0;
Gamma_0j(2:Nj,1) = lambda_j.*besseli(1,lambda_j.*Rp)./besseli(0,lambda_j.*Rp);

Gamma_1j(1,1) = 1/Rp;
Gamma_1j(2:Nj,1) = lambda_j.*besseli(2,lambda_j.*Rp)./besseli(1,lambda_j.*Rp) + 1/Rp;

Delta_0i(1,1:Ni) = -lambda_i.*(besselk(1,lambda_i.*Rp) ./ besselk(0,lambda_i.*Rp));
Delta_1i(1,1:Ni) = -lambda_i.*(besselk(2,lambda_i.*Rp)./besselk(1,lambda_i.*Rp)) + 1/Rp;

[S_0l_prim_Rp,S_0l_tild_prim_Rp] = Sprim_func(0,Rp,Rp,Rc,lambda_n);
[S_0l_prim_Rc,S_0l_tild_prim_Rc] = Sprim_func(0,Rc,Rp,Rc,lambda_n);

[S_1l_prim_Rp,S_1l_tild_prim_Rp] = Sprim_func(1,Rp,Rp,Rc,lambda_n);
[S_1l_prim_Rc,S_1l_tild_prim_Rc] = Sprim_func(1,Rc,Rp,Rc,lambda_n);

[S_0l_int,S_0l_tild_int] = Sint_func(0,Rp,Rc,lambda_n);
[S_1l_int,S_1l_tild_int] = Sint_func(1,Rp,Rc,lambda_n);

I1 = N_lambda_i.^(-.5).*( sin(lambda_i.*(h-e1)) - sin(lambda_i.*(h-e2)) ) ./ lambda_i;		%% Int_(-e2)^(-e1){ Z_i(z) dz }
I2 = N_lamnda_n.^(-.5).*sin(lambda_n.*e1)./lambda_n;											%% Int_(-e1)^(0){ Z_l(z) dz }
I3 = N_lambda_i.^(-.5).*f2(lambda_i,-e2,-e1,-Zc,h)';										%% Int_(-e2)^(-e1){ (z-Zc)*Z_i(z) dz }
I4 = N_lamnda_n.^(-.5).*f2(lambda_n,-e1,0,0,e1);												%% Int_(-e1)^(0){ z*Z_l }
I5 = N_lamnda_n.^(-.5).*f2(lambda_n,-e1,0,-Zc,e1);											%% Int_(-e1)^(0){ (z-Zc)*Z_l }
I6 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-e2,h,h)';											%% Int_(-h)^(-e2){ (z+h)^2*Z_i(z) dz }

%% Définition des matrices intrinsèques à la structure

D_r3_tau = zeros(Ni+Nn,Ni+Nn);
%%% interface r=Rp
D_r3_tau(1:Ni,1:Ni) = diag(Delta_0i) - ((h-e2)/h).*(L_jtau'*((Gamma_0j*ones(1,Ni)).*L_jtau)) - (e1/h).*(K_ltau'*(diag(S_0l_prim_Rp)*K_ltau));
D_r3_tau(1:Ni,Ni+1:end) = -(e1/h).*( diag(S_0l_tild_prim_Rp)*K_ltau )';
%%% interface r=Rc
D_r3_tau(1+Ni:end,1:Ni) = diag(S_0l_prim_Rc)*K_ltau;
D_r3_tau(1+Ni:end,1+Ni:end) = diag(S_0l_tild_prim_Rc);

D_r1_tau = zeros(Ni+Nn,Ni+Nn);
%%% interface r=Rp
D_r1_tau(1:Ni,1:Ni) = diag(Delta_1i) - ((h-e2)/h).*(L_jtau'*((Gamma_1j*ones(1,Ni)).*L_jtau)) - (e1/h).*(K_ltau'*(diag(S_1l_prim_Rp)*K_ltau));
D_r1_tau(1:Ni,Ni+1:end) = -(e1/h).*( diag(S_1l_tild_prim_Rp)*K_ltau )';
%%% interface r=Rc
D_r1_tau(Ni+1:end,1:Ni) = diag(S_1l_prim_Rc)*K_ltau;
D_r1_tau(Ni+1:end,Ni+1:end) = diag(S_1l_tild_prim_Rc);

D_r5_tau = D_r1_tau;


%% Problème de RADIATION
%% SURGE MODE

H_r1_tau = zeros(Ni+Nn,1);
%%% Interface r = Rp
H_r1_tau(1:Ni,1) = (1/h).*I1';
%%% Interface r = Rc
H_r1_tau(Ni+1:end,1) = (1/e1).*I2;

coeff = D_r1_tau\H_r1_tau;

A_r1 = coeff(1:Ni,1);
C_r1 = coeff(Ni+1:end,1);

B_r1 = K_ltau*A_r1;
F_r1 = L_jtau*A_r1;


Zr(1,1,w) = -pi*rho*(Rp*I1*A_r1 + Rc*I2'*C_r1);

tmp1 = .25*F_r1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*F_r1(2:end);
tmp2 = sum(N_lamnda_n.^(-.5).*(S_1l_int.*B_r1 + S_1l_tild_int.*C_r1));
Zr(1,3,w) = - pi*rho*( Rp*I3*A_r1 + Rc*I5'*C_r1 + tmp1 - tmp2);

% SURGE Wave excitation force --> Haskind's theorem
Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*A_r1(1)) / (cosh(k0*h)*besselh(1,k0*Rp));

%% HEAVE MODE

P_3j(1,1) = ( (h-e2) / 12 ) * ( 2 - 3*(Rp/(h-e2))^2);
P_3j(2:Nj,1) = (sqrt(2)*(h-e2).*(-1).^j) ./ (j.^2 .* pi^2);
P_3l = ( I4 + (g/omega^2).*N_lamnda_n.^(-.5).*sin(lambda_n.*e1)./lambda_n ) ./ e1;
H_r3_tau = zeros(Ni+Nn,1);
%%% Interface r = Rp
H_r3_tau(1:Ni) = -((h-e2)/h).*(L_jtau'*(Gamma_0j.*P_3j)) - (e1/h).*(K_ltau'*(S_0l_prim_Rp.*P_3l)) - (.5*Rp/h).*L_jtau(1,:)';
%%% Interface r = Rc
H_r3_tau(Ni+1:end) = S_0l_prim_Rc.*P_3l;

%%% Calcul des coefficients A_mi, C_mn
coeff = D_r3_tau\H_r3_tau;

A_r3 = coeff(1:Ni);
C_r3 = coeff(Ni+1:end);

B_r3 = K_ltau*A_r3 - P_3l;
F_r3 = L_jtau*A_r3 - P_3j;

tmp1 = .5*Rp^2*F_r3(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*F_r3(2:end);
tmp2 = .5*(h-e2)*( .5*Rp^2 - .125*Rp^4/(h-e2)^2 );
tmp3 = sum(N_lamnda_n.^(-.5).*S_0l_int.*B_r3) + sum(N_lamnda_n.^(-.5).*S_0l_tild_int.*C_r3);
tmp4 = .5*(g/omega^2 - e1)*(Rp^2-Rc^2);
Zr(2,2,w) = 2*pi*rho*(tmp1 + tmp2 - (tmp3 + tmp4));

% HEAVE Wave excitation force --> Haskind's theorem
Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*A_r3(1)) / (cosh(-imag(lambda_i(1))*h)*besselh(0,-imag(lambda_i(1))*Rp));



%%%                         PITCH MODE

P_5j(1,1) = (1/(24*(h-e2))) * ( 3*Rp^3 - 4*Rp*(h-e2)^2);
P_5j(2:Nj,1) = -(sqrt(2)*Rp*(h-e2).*(-1).^j) ./ (j.*pi).^2;


P_5l =  - (Rp/e1)*N_lamnda_n.^(-.5).*f2(lambda_n,-e1,0,g/omega^2,e1);

H_r5_tau = zeros(Ni+Nn,1);
%%% Interface r = Rp
tmp1 = N_lambda_i.^(-.5).*f2(lambda_i,-e1,0,g/omega^2,h)';
H_r5_tau(1:Ni,1) =  -((h-e2)/h).*(L_jtau'*(Gamma_1j.*P_5j)) - (e1/h).*(K_ltau'*(S_1l_prim_Rp.*P_5l))...
								+ ((3*Rp^2)/(8*h)).*L_jtau(1,:)' - (.5/(h*(h-e2))).*I6' - (1/h).*tmp1' + (1/h).*I3';
%%% Interface r = Rc
H_r5_tau(Ni+1:end,1) = S_1l_prim_Rc.*P_5l + (1/e1).*(I4 + (g/omega^2)*I2) + (1/e1).*I5;

coeff = D_r5_tau\H_r5_tau;

A_r5 = coeff(1:Ni,1);
C_r5 = coeff(Ni+1:end,1);

B_r5 = K_ltau*A_r5 - P_5l;
F_r5 = L_jtau*A_r5 - P_5j;

tmp1 = Rc^2*e1*(-(1/3)*e1^2 + .5*e1*(g/omega^2-Zc) + Zc*g/omega^2);
tmp2 = .25*F_r5(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*F_r5(2:end);
tmp3 = ( .125/(h-e2) ) * ( Rp^6/6 - Rp^4*(h-e2)^2);
tmp4 = sum(N_lamnda_n.^(-.5).*S_1l_int.*B_r5) + sum(N_lamnda_n.^(-.5).*S_1l_tild_int.*C_r5);
tmp5 = -.25*(g/omega^2 - e1)*(Rp^4-Rc^4);
Zr(3,3,w) = - pi*rho*( Rp*I3*A_r5 +Rc*I5'*C_r5 + tmp1 + (tmp2 + tmp3) - (tmp4 + tmp5));

tmp1 = Rc^2*(.5*e1^2 - (g/omega^2)*e1);
Zr(3,1,w) = -pi*rho*(Rp*I1*A_r5 + Rc*I2'*C_r5 + tmp1);

% PITCH Wave excitation force --> Haskind's theorem
Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*A_r5(1)) / (cosh(k0*h)*besselh(1,k0*Rp));


%% Problème de DIFFRACTION

a0 = 1i*omega;

B0 = (-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));
B1 = 2*1i*B0;



H_d_0tau = zeros(Ni+Nn,1);
H_d_0tau(1:Ni,1) = B0.*( besselj(0,k0*Rp).*( ((h-e2)/h).*(L_jtau'*(Gamma_0j.*L_jtau(:,1)))...
		+ (e1/h).*(K_ltau'*(S_0l_prim_Rp.*K_ltau(:,1)))) + k0*besselj(1,k0*Rp)*[1;zeros(Ni-1,1)]);
H_d_0tau(1+Ni:end,1) = -B0*besselj(0,k0*Rp).*(S_0l_prim_Rc.*K_ltau(:,1));

coeff = D_r3_tau\H_d_0tau;

a_0i = coeff(1:Ni,1);
c_0l = coeff(Ni+1:end,1);

b_0l = K_ltau*a_0i + B0.*besselj(0,k0*Rp).*K_ltau(:,1);
f_0j = L_jtau*a_0i + B0.*besselj(0,k0*Rp).*L_jtau(:,1);    

tmp1 = .5*Rp^2*f_0j(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*f_0j(2:end);
tmp2 = sum(N_lamnda_n.^(-.5).*S_0l_int.*b_0l) + sum(N_lamnda_n.^(-.5).*S_0l_tild_int.*c_0l);
Fe(w,2) = 2*pi*a0*rho*(tmp1 - tmp2);

Fe_FK(w,2) = 2*pi*a0*rho*B0*N_lambda_i(1)^(-.5)*(cosh(k0*(h-e2))*(Rp*besselj(1,k0*Rp)) - cosh(k0*(h-e1))*(Rp*besselj(1,k0*Rp) - Rc*besselj(1,k0*Rc)))/k0;



%% m=1
H_d_1tau = zeros(Ni+Nn,1);
H_d_1tau(1:Ni) = B1.*( besselj(1,k0*Rp).*( ((h-e2)/h).*(L_jtau'*(Gamma_1j.*L_jtau(:,1)))...
		+ (e1/h).*(K_ltau'*(S_1l_prim_Rp.*K_ltau(:,1)))) - (k0*besselj(0,k0*Rp)-besselj(1,k0*Rp)/Rp).*[1;zeros(Ni-1,1)]);
H_d_1tau(1+Ni:end) = -B1*besselj(1,k0*Rp).*(S_1l_prim_Rc.*K_ltau(:,1));

coeff = D_r1_tau\H_d_1tau;

a_1i = coeff(1:Ni);
c_1l = coeff(Ni+1:end);

b_1l = K_ltau*a_1i + B1.*besselj(1,k0*Rp).*K_ltau(:,1);
f_1j = L_jtau*a_1i + B1.*besselj(1,k0*Rp).*L_jtau(:,1);    


Fe(w,1) = - pi*a0*rho*( Rp*(B1*besselj(1,k0*Rp)*I1(1) + I1*a_1i) + Rc*I2'*c_1l);

Fe_FK(w,1) = -pi*a0*rho*B1*(Rp*besselj(1,k0*Rp)*I1(1)...
	+ Rc*besselj(1,k0*Rc)*N_lambda_i(1)^(-.5)*(sin(lambda_i(1)*h) - sin(lambda_i(1)*(h-e1)))/lambda_i(1));


tmp1 = .25*f_1j(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*f_1j(2:end);
tmp2 = sum(N_lamnda_n.^(-.5).*S_1l_int.*b_1l) + sum(N_lamnda_n.^(-.5).*S_1l_tild_int.*c_1l);
Fe(w,3) = - pi*a0*rho*( Rp*(B1*besselj(1,k0*Rp)*I3(1) + I3*a_1i) + Rc*(I5'*c_1l) + tmp1 - tmp2);

tmp1 = Rp*besselj(1,k0*Rp)*N_lambda_i(1)^(-.5)*f2(lambda_i(1),-e2,-e1,-Zc,h);
tmp2 = Rc*besselj(1,k0*Rc)*N_lambda_i(1)^(-.5)*f2(lambda_i(1),-e1,0,-Zc,h);
tmp3 = N_lambda_i(1)^(-.5)*cos(lambda_i(1)*(h-e2))*(Rp^2*besselj(2,k0*Rp))/k0;
tmp4 = N_lambda_i(1)^(-.5)*cos(lambda_i(1)*(h-e1))*(Rp^2*besselj(2,k0*Rp)-Rc^2*besselj(2,k0*Rc))/k0;
Fe_FK(w,3) = -pi*a0*rho*B1*(tmp1 + tmp2 + tmp3 - tmp4);


end%% END OF FOR() LOOP
%% EXIT
for i=1:3
	for j=1:3
		A(i,j,:) = real(Zr(i,j,1:end-1));
		B(i,j,:) = squeeze(imag(Zr(i,j,1:end-1))).*Omega;
	end
end
A_inf = real(Zr(:,:,end));
% ATTENTION -> ici on reprend une dépendance en temps positive exp(-iwt) --> exp(iwt)
Fe = conj(Fe(1:end-1,:));
Fe_Haskind = conj(Fe_Haskind(1:end-1,:));
Fe_FK = conj(Fe_FK(1:end-1,:));


save([DIR,filesep,'DATA/hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');

end%% END OF FUNCTION

function out = f1( alpha, a, b, c, d )
%%% Int_(a)^(b){(z+c)^2*cos(alpha*(z+d))}
epsilon = 1e-8;
for i=1:length(alpha)
if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
	out(i,1) = ((b+c)^3 - (a+c)^3) / 3;
else
	out(i,1) = (((b+c)^2 - 2./alpha(i).^2).*sin(alpha(i).*(b+d)) - ((a+c)^2 - 2./alpha(i).^2).*sin((i).*(a+d)))./alpha(i)...
				+ 2.*((b+c).*cos(alpha(i).*(b+d)) - (a+c).*cos(alpha(i).*(a+d)))./alpha(i).^2;
end
end
end

function out = f2(alpha, a, b, c, d )
%%% Int_(a)^(b){(z+c)*cos(alpha*(z+d))}
epsilon = 1e-8;
for i=1:length(alpha)
if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
	out(i,1) = .5*((b+c)^2 - (a+c)^2);
else
	out(i,1) = ((b+c).*sin(alpha(i).*(b+d)) - (a+c).*sin(alpha(i).*(a+d)))./alpha(i)  + (cos(alpha(i).*(b+d)) - cos(alpha(i).*(a+d)))./alpha(i).^2;
end
end
end

function [S_prim,S_tild_prim] = Sprim_func(m,r,R1,R2,alpha)
	
cst = besseli(m,alpha.*R1).*besselk(m,alpha.*R2) - besseli(m,alpha.*R2).*besselk(m,alpha.*R1);

S_prim = (alpha./cst).*( besselk(m,alpha.*R2).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r))+...
						  besseli(m,alpha.*R2).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));

S_tild_prim = -(alpha./cst).*( besselk(m,alpha.*R1).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r))+...
						  besseli(m,alpha.*R1).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));
end

function [S_int,S_tild_int] = Sint_func(m,R1,R2,alpha)
	cst = besseli(m,alpha.*R1).*besselk(m,alpha.*R2) - besseli(m,alpha.*R2).*besselk(m,alpha.*R1);

	S_int = ( besselk(m,alpha.*R2).*(R1^(m+1).*besseli(m+1,alpha.*R1)-R2^(m+1).*besseli(m+1,alpha.*R2)) +...
			  besseli(m,alpha.*R2).*(R1^(m+1).*besselk(m+1,alpha.*R1)-R2^(m+1).*besselk(m+1,alpha.*R2)) )./(alpha.*cst);

	S_tild_int = -( besselk(m,alpha.*R1).*(R1^(m+1).*besseli(m+1,alpha.*R1)-R2^(m+1).*besseli(m+1,alpha.*R2)) +...
					besseli(m,alpha.*R1).*(R1^(m+1).*besselk(m+1,alpha.*R1)-R2^(m+1).*besselk(m+1,alpha.*R2)) )./(alpha.*cst);    
end



