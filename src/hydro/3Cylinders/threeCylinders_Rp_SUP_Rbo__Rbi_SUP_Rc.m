function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = threeCylinders_Rp_SUP_Rbo__Rbi_SUP_Rc(Omega, depth, WECstructure, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the hydrodynamic coefficients for the three cylinder
% structure depicted below (Rp>Rbo).
% 
% Inputs :
% - Omega				: Vector of wave frequencies (rad/s)
% - depth				: Water depth (m)
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
%           |--------------------------> (Rbo)
%           |--------> (Rbi)
%           |-----> (Rc)
%           .                                                              
% --------  +-----+  +-----------------+ --------------------------- (z=0)
%           .ooooo|  |ccccccccccccccccc|  ^      ^      ^        ^
%           |ooooo|  |ccccccccccccccccc|  |      |      |        |
%           .ooooo|  |ccccccccccccccccc|  |      |      |        |
%           |ooooo|  +-----------------+  v (b)  |      |        |
%           .ooooo|                              |      |        |
%           |ooooo|                              |      |        |
%           .ooooo|                              |      |        |
%           |ooooo|                              |      |        |
%           .ooooo|                              |      |        |
%           |ooooo|                              |      |        |
%           .ooooo|                              |      |        |
%           |ooooo|                              |      |        |
%           .ooooo|                              |      |        |
%           |ooooo|                              |      |        |
%           .ooooo|                              |      |        |
%           |ooooo|                              |      |        |
%           .ooooo|                              v (e1) |        |
%           |ooooo+-----------------------------+       |        |
%           .ooooo|ooooooooooooooooooooooooooooo|       |        |
%           |ooooo|ooooooooooooooooooooooooooooo|       |        |
%           +-----+-----------------------------+       v (e2)   .
%           .                                                    .
%           |-----------------------------------> (Rp)           .
%                                                                |
%                                                                |
%                                                                |
%                                                                v (h)
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

g = 9.81;
rho = 1025;
h = depth; %%% water depth

Rbo = WECstructure.Rbo; %%% rayon exterieur flotteur
Rbi = WECstructure.Rbi; %%% rayon interieur flotteur
Rc = WECstructure.Rc;
Rp = WECstructure.Rp;

b = WECstructure.b;
e1 = WECstructure.e1;
e2 = WECstructure.e2;

% alpha= WECstructure.alpha;

if ~isfield(options,'Trunc')
	Ni=80; Nl=80; Nn=100; Nj=100;
else
	Ni = options.Trunc.Ni;
	Nl = options.Trunc.Nl;
	Nn = options.Trunc.Nn;
	Nj = options.Trunc.Nj;
end

if ~isfield(options,'Zc')
	zG=[0;0];
else
	zG=options.Zc;
end

epsilon = 1e-7;

%%% Initialisation des sorties
Zr = zeros(6,6,length(Omega)); %%% matrice d'impédance --> A + (1i/omega)*B
A = zeros(6,6,length(Omega));
A_inf = zeros(6,6);
B = zeros(6,6,length(Omega));
Fe = zeros(length(Omega),6); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),6); %%% effort d'excitation
Fe_FK = zeros(length(Omega),6);%%% Froude-Krilov	


n = (1:1:Nn-1)';
alpha_n = (pi/(e1-b)).*n;
j = (1:1:Nj-1)';
beta_j = (pi/(h-e2)).*j;

for w=1:length(Omega)+1

	clc;
	disp([num2str(round(w*100/(length(Omega)+1))),'%']);
	
	if w<=length(Omega)
		omega = Omega(w);
		lambda_i=fcn_waveDispersion(Omega(w),h,Ni);
		gamma_l=fcn_waveDispersion(Omega(w),e1,Nl); gamma_l=conj(gamma_l)';

		k_e = -imag(lambda_i(1));
	else
		omega = Inf;% omega >> 1
		i=1:1:Ni;
		l=(1:1:Nl)';
		lambda_i = .5*(2*i-1)*pi/h;
		gamma_l = .5*(2*l-1)*pi/e1;
	end
	
    N_lambda_i = .5.*( 1 + sin(2.*lambda_i.*h)./(2.*lambda_i.*h));
    N_gamma_l = .5.*( 1 + sin(2.*gamma_l.*e1)./(2.*gamma_l.*e1));
	
	L_jtau = zeros(Nj,Ni);
	L_jtau(1,:) = N_lambda_i.^(-.5).*sin(lambda_i.*(h-e2))./( (h-e2).*lambda_i );
	for i = 1 : Ni
		for p = 1 : Nj-1
		    if (beta_j(p)/lambda_i(i) < 1+epsilon)&&(beta_j(p)/lambda_i(i) > 1-epsilon)
		        L_jtau(p+1,i) = .5*sqrt(2).*N_lambda_i(i)^(-.5);
		    else
		        L_jtau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*(sin((beta_j(p)-lambda_i(i))*(h-e2))/(beta_j(p)-lambda_i(i))...
					+ sin((beta_j(p)+lambda_i(i))*(h-e2))/(beta_j(p)+lambda_i(i)))/(h-e2);
		    end
		end
	end

    K_ltau = zeros(Nl,Ni);
	for i = 1 : Ni
		for p = 1 : Nl
			if (gamma_l(p)/lambda_i(i) < 1+epsilon)&&(gamma_l(p)/lambda_i(i) > 1-epsilon)
                K_ltau(p,i) = (1/e1)*(N_lambda_i(i)^(-.5)*N_gamma_l(p).^(-.5))*...
					(.25*(sin(lambda_i(i)*(h+e1)) - sin(lambda_i(i)*(h-e1)))/lambda_i(i) + e1*.5*cos(lambda_i(i)*(h-e1)));
			else
                K_ltau(p,i) = (1/e1)*(N_lambda_i(i)^(-.5)*N_gamma_l(p).^(-.5))*...
					.5*( (sin(gamma_l(p)*e1-lambda_i(i)*h)-sin(gamma_l(p)*e1-lambda_i(i)*h - e1*(gamma_l(p)-lambda_i(i))))/(gamma_l(p)-lambda_i(i))...
							+ (sin(gamma_l(p)*e1+lambda_i(i)*h)-sin(gamma_l(p)*e1+lambda_i(i)*h - e1*(gamma_l(p)+lambda_i(i))))/(gamma_l(p)+lambda_i(i)));
			end
		end
	end

	M_ntau = zeros(Nn,Nl);
	M_ntau(1,:) = N_gamma_l.^(-.5).*sin(gamma_l.*(e1-b))./((e1-b).*gamma_l);
	for i = 1 : Nl
		for p = 1 : Nn-1
		    if (alpha_n(p)/gamma_l(i) < 1+epsilon)&&(alpha_n(p)/gamma_l(i) > 1-epsilon)
		        M_ntau(p+1,i) = .5*sqrt(2).*N_gamma_l(i)^(-.5);
		    else
		        M_ntau(p+1,i) = .5*sqrt(2)*N_gamma_l(i)^(-.5)*...
					( sin((alpha_n(p)-gamma_l(i))*(e1-b))/(alpha_n(p)-gamma_l(i))...
					+ sin((alpha_n(p)+gamma_l(i))*(e1-b))/(alpha_n(p)+gamma_l(i)) ) / (e1-b);
		    end
		end
	end


	Gamma_0j(1,1) = 0;
	Gamma_0j(2:Nj,1) = beta_j.*besseli(1,beta_j.*Rp)./besseli(0,beta_j.*Rp);

	Gamma_1j(1,1) = 1/Rp;
	Gamma_1j(2:Nj,1) = beta_j.*besseli(2,beta_j.*Rp)./besseli(1,beta_j.*Rp) + 1/Rp;

	Delta_0i(1,1:Ni) = -lambda_i.*(besselk(1,lambda_i.*Rp) ./ besselk(0,lambda_i.*Rp));
	Delta_1i(1,1:Ni) = -lambda_i.*(besselk(2,lambda_i.*Rp)./besselk(1,lambda_i.*Rp)) + 1/Rp;

    [S_0l_prim_Rp,S_0l_tild_prim_Rp] = Sprim_func(0,Rp,Rp,Rbo,gamma_l);
    [S_0l_prim_Rbo,S_0l_tild_prim_Rbo] = Sprim_func(0,Rbo,Rp,Rbo,gamma_l);

    [S_1l_prim_Rp,S_1l_tild_prim_Rp] = Sprim_func(1,Rp,Rp,Rbo,gamma_l);
    [S_1l_prim_Rbo,S_1l_tild_prim_Rbo] = Sprim_func(1,Rbo,Rp,Rbo,gamma_l);

    [T_0n_prim_Rbo,T_0n_tild_prim_Rbo] = Tprim_func(0,Rbo,Rbo,Rbi,alpha_n);
    [T_0n_prim_Rbi,T_0n_tild_prim_Rbi] = Tprim_func(0,Rbi,Rbo,Rbi,alpha_n);
	
    [T_1n_prim_Rbo,T_1n_tild_prim_Rbo] = Tprim_func(1,Rbo,Rbo,Rbi,alpha_n);
	[T_1n_prim_Rbi,T_1n_tild_prim_Rbi] = Tprim_func(1,Rbi,Rbo,Rbi,alpha_n);

	
    [N_0l_prim_Rbi,N_0l_tild_prim_Rbi] = Sprim_func(0,Rbi,Rbi,Rc,gamma_l);
    [N_0l_prim_Rc,N_0l_tild_prim_Rc] = Sprim_func(0,Rc,Rbi,Rc,gamma_l);

    [N_1l_prim_Rbi,N_1l_tild_prim_Rbi] = Sprim_func(1,Rbi,Rbi,Rc,gamma_l);
    [N_1l_prim_Rc,N_1l_tild_prim_Rc] = Sprim_func(1,Rc,Rbi,Rc,gamma_l);
	
	
	%%% Intégrale
	[S_0l_int,S_0l_tild_int] = Sint_func(0,Rp,Rbo,gamma_l);
    [S_1l_int,S_1l_tild_int] = Sint_func(1,Rp,Rbo,gamma_l);

    [T_0n_int,T_0n_tild_int] = Ti_func(0,Rbi,Rbo,alpha_n);
    [T_1n_int,T_1n_tild_int] = Ti_func(1,Rbi,Rbo,alpha_n);

	[N_0l_int,N_0l_tild_int] = Sint_func(0,Rbi,Rc,gamma_l);
    [N_1l_int,N_1l_tild_int] = Sint_func(1,Rbi,Rc,gamma_l);

	
    I1 = N_gamma_l.^(-.5).*( sin(gamma_l.*e1)-sin(gamma_l.*(e1-b)) ) ./ gamma_l;					%% Int_(-b)^(0){ Z_l(z) dz }
    I2 = N_lambda_i.^(-.5).*( sin(lambda_i.*(h-e1)) - sin(lambda_i.*(h-e2)) ) ./ lambda_i;			%% Int_(-e2)^(-e1){ Z_i^I }
	I10 = N_gamma_l.^(-.5).*sin(gamma_l.*e1)./gamma_l;												%% Int_(-e1)^(0){ Z_l(z) dz }
    I3 = N_gamma_l.^(-.5).*f2(gamma_l,-b,0,-zG(1),e1);												%% Int_(-b)^(0){ (z-Zc)*Z_l^II }
    I4 = N_lambda_i.^(-.5).*f2(lambda_i,-e2,-e1,-zG(2),h)';											%% Int_(-e2)^(-e1){ (z-Zc)*Z_i^I }
    I5 = N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,0,e1);													%% Int_(-e1)^(0){ z*Z_l^II }
    I6 = N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,-zG(2),e1);												%% Int_(-e1)^(0){ (z-Zc)*Z_l(z) dz }
    I7_k1 = N_gamma_l.^(-.5).*f1(gamma_l,-e1,-b,e1,e1);												%% Int_(-e1)^(-b){ (z+e1)^2*Z_l^II }
    I7_k2 = N_gamma_l.^(-.5).*f1(gamma_l,-e1,-b,b,e1);												%% Int_(-e1)^(-b){ (z+b)^2*Z_l^II }
	I9 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-e2,h,h);												%% Int_(-h)^(-e2){ (z+h)^2*Z_i^I }

	%%% ***********************************************************************
    %%% Définition des matrices intrinsèques à la structure
	%%% ***********************************************************************
	
	D_m0 = zeros(Ni+3*Nl,Ni+3*Nl);
	
	%%% Interface r = Rp
    D_m0(1:Ni,1:Ni) = diag(Delta_0i) - ((h-e2)/h).*(L_jtau'*((Gamma_0j*ones(1,Ni)).*L_jtau)) - (e1/h).*(K_ltau'*(diag(S_0l_prim_Rp)*K_ltau));
    D_m0(1:Ni,1+Ni:Ni+Nl) = -(e1/h).*( diag(S_0l_tild_prim_Rp)*K_ltau )';
    D_m0(1:Ni,1+Ni+Nl:Ni+2*Nl) = zeros(Ni,Nl);
	D_m0(1:Ni,1+Ni+2*Nl:end) = zeros(Ni,Nl);
	
	%%% Interface r = Rbo
    D_m0(1+Ni:Ni+Nl,1:Ni) = diag(S_0l_prim_Rbo)*K_ltau;
    D_m0(1+Ni:Ni+Nl,1+Ni:Ni+Nl) = diag(S_0l_tild_prim_Rbo) - ((e1-b)/e1).*(M_ntau'*(diag(T_0n_prim_Rbo)*M_ntau));
    D_m0(1+Ni:Ni+Nl,1+Ni+Nl:Ni+2*Nl) = - ((e1-b)/e1).*(M_ntau'*(diag(T_0n_tild_prim_Rbo)*M_ntau));
	D_m0(1+Ni:Ni+Nl,1+Ni+2*Nl:end) = zeros(Nl,Nl);

	%%% Interface r = Rbi
    D_m0(1+Ni+Nl:Ni+2*Nl,1:Ni) = zeros(Nl,Ni);
    D_m0(1+Ni+Nl:Ni+2*Nl,1+Ni:Ni+Nl) = -((e1-b)/e1).*(M_ntau'*(diag(T_0n_prim_Rbi)*M_ntau));
    D_m0(1+Ni+Nl:Ni+2*Nl,1+Ni+Nl:Ni+2*Nl) = diag(N_0l_prim_Rbi) - ((e1-b)/e1).*(M_ntau'*(diag(T_0n_tild_prim_Rbi)*M_ntau));
	D_m0(1+Ni+Nl:Ni+2*Nl,1+Ni+2*Nl:end) = diag(N_0l_tild_prim_Rbi);
	
	%%% Interface r = Rc
    D_m0(1+Ni+2*Nl:end,1:Ni) = zeros(Nl,Ni);
    D_m0(1+Ni+2*Nl:end,1+Ni:Ni+Nl) = zeros(Nl,Nl);
    D_m0(1+Ni+2*Nl:end,1+Ni+Nl:Ni+2*Nl) = diag(N_0l_prim_Rc);
	D_m0(1+Ni+2*Nl:end,1+Ni+2*Nl:end) = diag(N_0l_tild_prim_Rc);
	

    % ***
	D_m1 = zeros(Ni+3*Nl,Ni+3*Nl);
	
	%%% Interface r = Rp
    D_m1(1:Ni,1:Ni) = diag(Delta_1i) - ((h-e2)/h).*(L_jtau'*((Gamma_1j*ones(1,Ni)).*L_jtau)) - (e1/h).*(K_ltau'*(diag(S_1l_prim_Rp)*K_ltau));
    D_m1(1:Ni,1+Ni:Ni+Nl) = -(e1/h).*( diag(S_1l_tild_prim_Rp)*K_ltau )';
    D_m1(1:Ni,1+Ni+Nl:Ni+2*Nl) = zeros(Ni,Nl);
	D_m1(1:Ni,1+Ni+2*Nl:end) = zeros(Ni,Nl);
	
	%%% Interface r = Rbo
    D_m1(1+Ni:Ni+Nl,1:Ni) = diag(S_1l_prim_Rbo)*K_ltau;
    D_m1(1+Ni:Ni+Nl,1+Ni:Ni+Nl) = diag(S_1l_tild_prim_Rbo) - ((e1-b)/e1).*(M_ntau'*(diag(T_1n_prim_Rbo)*M_ntau));
    D_m1(1+Ni:Ni+Nl,1+Ni+Nl:Ni+2*Nl) = - ((e1-b)/e1).*(M_ntau'*(diag(T_1n_tild_prim_Rbo)*M_ntau));
	D_m1(1+Ni:Ni+Nl,1+Ni+2*Nl:end) = zeros(Nl,Nl);

	%%% Interface r = Rbi
    D_m1(1+Ni+Nl:Ni+2*Nl,1:Ni) = zeros(Nl,Ni);
    D_m1(1+Ni+Nl:Ni+2*Nl,1+Ni:Ni+Nl) = - ((e1-b)/e1).*(M_ntau'*(diag(T_1n_prim_Rbi)*M_ntau));
    D_m1(1+Ni+Nl:Ni+2*Nl,1+Ni+Nl:Ni+2*Nl) = diag(N_1l_prim_Rbi) - ((e1-b)/e1).*(M_ntau'*(diag(T_1n_tild_prim_Rbi)*M_ntau));
	D_m1(1+Ni+Nl:Ni+2*Nl,1+Ni+2*Nl:end) = diag(N_1l_tild_prim_Rbi);
	
	%%% Interface r = Rc
    D_m1(1+Ni+2*Nl:end,1:Ni) = zeros(Nl,Ni);
    D_m1(1+Ni+2*Nl:end,1+Ni:Ni+Nl) = zeros(Nl,Nl);
    D_m1(1+Ni+2*Nl:end,1+Ni+Nl:Ni+2*Nl) = diag(N_1l_prim_Rc);
	D_m1(1+Ni+2*Nl:end,1+Ni+2*Nl:end) = diag(N_1l_tild_prim_Rc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  Problème de RADIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         SURGE MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
%%% ***********************************************************************                  
%%%     k = 1       Buoy oscillating, column fixed
%%% ***********************************************************************

	H_q1_k1 = zeros(Ni+3*Nl,1);
    H_q1_k1(1:Ni) = zeros(Ni,1);
    H_q1_k1(1+Ni:Ni+Nl) = (1/e1).*I1;
    H_q1_k1(1+Ni+Nl:Ni+2*Nl) = (1/e1).*I1;
	H_q1_k1(1+Ni+2*Nl:Ni+3*Nl) = zeros(Nl,1);
  
    coeff = D_m1\H_q1_k1;

    A_r1_k1 = coeff(1:Ni);
    C_r1_k1 = coeff(1+Ni:Ni+Nl);
    G_r1_k1 = coeff(1+Ni+Nl:Ni+2*Nl);
	H_r1_k1 = coeff(1+Ni+2*Nl:end);
   
    B_r1_k1 = K_ltau*A_r1_k1;
	D_r1_k1 = M_ntau*C_r1_k1;
    E_r1_k1 = M_ntau*G_r1_k1;
    F_r1_k1 = L_jtau*A_r1_k1;

% 	if w<=length(Omega)
		Zr(1,1,w) = -pi*rho*(Rbo*I1'*C_r1_k1 - Rbi*I1'*G_r1_k1);
		Zr(1,4,w) = -pi*rho*(Rp*I2*A_r1_k1 + Rc*I10'*H_r1_k1);

		tmp1 = T_1n_int(1)*D_r1_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r1_k1(2:end));
		tmp2 = T_1n_tild_int(1)*E_r1_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r1_k1(2:end));
		Zr(1,3,w) = - pi*rho*(Rbo*I3'*C_r1_k1 - Rbi*I3'*G_r1_k1 + tmp1 + tmp2);

		tmp1 = .25*F_r1_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k1(2:end);
		tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r1_k1) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r1_k1);
		tmp3 = T_1n_int(1)*D_r1_k1(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r1_k1(2:end));
		tmp4 = T_1n_tild_int(1)*E_r1_k1(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r1_k1(2:end));
		tmp5 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r1_k1) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r1_k1);
		Zr(1,6,w) = - pi*rho*( Rp*I4*A_r1_k1 + Rc*I6'*H_r1_k1 + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5) );
% 	else
% 		A_inf(1,1) = -pi*rho*(Rbo*I1'*C_r1_k1 - Rbi*I1'*G_r1_k1);
% 		A_inf(1,4) = -pi*rho*(Rp*I2*A_r1_k1 + Rc*I10'*H_r1_k1);
% 
% 		tmp1 = T_1n_int(1)*D_r1_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r1_k1(2:end));
% 		tmp2 = T_1n_tild_int(1)*E_r1_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r1_k1(2:end));
% 		A_inf(1,3) = - pi*rho*(Rbo*I3'*C_r1_k1 - Rbi*I3'*G_r1_k1 + tmp1 + tmp2);
% 
% 		tmp1 = .25*F_r1_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k1(2:end);
% 		tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r1_k1) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r1_k1);
% 		tmp3 = T_1n_int(1)*D_r1_k1(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r1_k1(2:end));
% 		tmp4 = T_1n_tild_int(1)*E_r1_k1(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r1_k1(2:end));
% 		tmp5 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r1_k1) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r1_k1);
% 		A_inf(1,6) = - pi*rho*( Rp*I4*A_r1_k1 + Rc*I6'*H_r1_k1 + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5) );
% 	end

%%% ***********************************************************************
%%%     k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************

	H_q1_k2 = zeros(Ni+3*Nl,1);
    H_q1_k2(1:Ni) = (1/h).*I2';
    H_q1_k2(1+Ni:Ni+Nl) = zeros(Nl,1);
    H_q1_k2(1+Ni+Nl:Ni+2*Nl) = zeros(Nl,1);
	H_q1_k2(1+Ni+2*Nl:Ni+3*Nl) = (1/e1).*I10;
 
    coeff = D_m1\H_q1_k2;

    A_r1_k2 = coeff(1:Ni);
    C_r1_k2 = coeff(1+Ni:Ni+Nl);
    G_r1_k2 = coeff(1+Ni+Nl:Ni+2*Nl);
	H_r1_k2 = coeff(1+Ni+2*Nl:end);

    B_r1_k2 = K_ltau*A_r1_k2;
	D_r1_k2 = M_ntau*C_r1_k2;
    E_r1_k2 = M_ntau*G_r1_k2;
    F_r1_k2 = L_jtau*A_r1_k2; 
    
% 	if w<=length(Omega)
		Zr(4,1,w) = -pi*rho*(Rbo*I1'*C_r1_k2 - Rbi*I1'*G_r1_k2);

		Zr(4,4,w) = -pi*rho*(Rp*I2*A_r1_k2 + Rc*I10'*H_r1_k2);

		tmp1 = T_1n_int(1)*D_r1_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r1_k2(2:end));
		tmp2 = T_1n_tild_int(1)*E_r1_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r1_k2(2:end));
		Zr(4,3,w) = - pi*rho*(Rbo*I3'*C_r1_k2 - Rbi*I3'*G_r1_k2 + tmp1 + tmp2);

		tmp1 = .25*F_r1_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k2(2:end);
		tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r1_k2) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r1_k2);
		tmp3 = T_1n_int(1)*D_r1_k2(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r1_k2(2:end));
		tmp4 = T_1n_tild_int(1)*E_r1_k2(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r1_k2(2:end));
		tmp5 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r1_k2) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r1_k2);
		Zr(4,6,w) = - pi*rho*( Rp*I4*A_r1_k2 + Rc*I6'*H_r1_k2 + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5));

		Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*A_r1_k1(1)) / (cosh(k_e*h)*besselh(1,k_e*Rp));
		Fe_Haskind(w,4) = (4*rho*g*h*sqrt(N_lambda_i(1))*A_r1_k2(1)) / (cosh(k_e*h)*besselh(1,k_e*Rp));
% 	else
% 		A_inf(4,1) = -pi*rho*(Rbo*I1'*C_r1_k2 - Rbi*I1'*G_r1_k2);
% 		A_inf(4,4) = -pi*rho*(Rp*I2*A_r1_k2 + Rc*I10'*H_r1_k2);
% 
% 		tmp1 = T_1n_int(1)*D_r1_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r1_k2(2:end));
% 		tmp2 = T_1n_tild_int(1)*E_r1_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r1_k2(2:end));
% 		A_inf(4,3) = - pi*rho*(Rbo*I3'*C_r1_k2 - Rbi*I3'*G_r1_k2 + tmp1 + tmp2);
% 
% 		tmp1 = .25*F_r1_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k2(2:end);
% 		tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r1_k2) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r1_k2);
% 		tmp3 = T_1n_int(1)*D_r1_k2(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r1_k2(2:end));
% 		tmp4 = T_1n_tild_int(1)*E_r1_k2(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r1_k2(2:end));
% 		tmp5 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r1_k2) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r1_k2);
% 		A_inf(4,6) = - pi*rho*( Rp*I4*A_r1_k2 + Rc*I6'*H_r1_k2 + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5));
% 	end
% end%%% END OF SURGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         HEAVE MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if strcmp(options.dof,'all') || strcmp(options.dof,'heave')
%%% ***********************************************************************
%%%		k = 1	-->	Buoy oscillating, column fixed
%%% ***********************************************************************
    % % solution particulière pour r = R
    P_3n_k1(1,1) = ( (e1-b) / 12 ) * ( 2 - 3*(Rbo^2/(e1-b)^2) );
    P_3n_k1(2:Nn,1) = (sqrt(2)*(e1-b).*(-1).^n)./(n.*pi).^2;
		
    % % solution particulière pour r = R
    Q_3n_k1(1,1) = ( (e1-b) / 12 ) * ( 2 - 3*(Rbi^2/(e1-b)^2) );
    Q_3n_k1(2:Nn,1) = (sqrt(2)*(e1-b).*(-1).^n) ./ (n.*pi).^2;
	
	H_q3_k1 = zeros(Ni+3*Nl,1);
    H_q3_k1(1:Ni) = zeros(Ni,1);
    H_q3_k1(1+Ni:Ni+Nl) = -((e1-b)/e1).*M_ntau'*(T_0n_prim_Rbo.*P_3n_k1 + T_0n_tild_prim_Rbo.*Q_3n_k1) - (.5*Rbo/e1).*M_ntau(1,:)';
    H_q3_k1(1+Ni+Nl:Ni+2*Nl) = -((e1-b)/e1).*M_ntau'*(T_0n_prim_Rbi.*P_3n_k1 + T_0n_tild_prim_Rbi.*Q_3n_k1) - (.5*Rbi/e1).*M_ntau(1,:)';
	H_q3_k1(1+Ni+2*Nl:Ni+3*Nl) = zeros(Nl,1);

    coeff = D_m0\H_q3_k1;

    A_r3_k1 = coeff(1:Ni);
    C_r3_k1 = coeff(1+Ni:Ni+Nl);
    G_r3_k1 = coeff(1+Ni+Nl:Ni+2*Nl);
	H_r3_k1 = coeff(1+Ni+2*Nl:end);
   
    B_r3_k1 = K_ltau*A_r3_k1;
	D_r3_k1 = M_ntau*C_r3_k1 - P_3n_k1;
    E_r3_k1 = M_ntau*G_r3_k1 - Q_3n_k1;
    F_r3_k1 = L_jtau*A_r3_k1;

% 	if w<=length(Omega)
		tmp1 = .25*(e1-b)*( (Rbo^2-Rbi^2) - (1/4)*( (Rbo^4-Rbi^4)/(e1-b)^2 ));                
		tmp2 = T_0n_int(1)*D_r3_k1(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*D_r3_k1(2:end));
		tmp3 = T_0n_tild_int(1)*E_r3_k1(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*E_r3_k1(2:end));
		Zr(2,2,w) = 2*pi*rho*( tmp1 + tmp2 + tmp3);

		tmp1 = .5*Rp^2*F_r3_k1(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k1(2:end);
		tmp2 = sum(N_gamma_l.^(-.5).*S_0l_int.*B_r3_k1) + sum(N_gamma_l.^(-.5).*S_0l_tild_int.*C_r3_k1);
		tmp3 = T_0n_int(1)*D_r3_k1(1) + sqrt(2)*sum(T_0n_int(2:end).*D_r3_k1(2:end));
		tmp4 = T_0n_tild_int(1)*E_r3_k1(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*E_r3_k1(2:end));
		tmp5 = sum(N_gamma_l.^(-.5).*N_0l_int.*G_r3_k1) + sum(N_gamma_l.^(-.5).*N_0l_tild_int.*H_r3_k1);
		Zr(2,5,w) = 2*pi*rho*( tmp1 - (tmp2 + tmp3 + tmp4 + tmp5 - (Rbo^4 - Rbi^4)/(16*(e1-b)) ) );
% 	else
% 		tmp1 = .25*(e1-b)*( (Rbo^2-Rbi^2) - (1/4)*( (Rbo^4-Rbi^4)/(e1-b)^2 ));                
% 		tmp2 = T_0n_int(1)*D_r3_k1(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*D_r3_k1(2:end));
% 		tmp3 = T_0n_tild_int(1)*E_r3_k1(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*E_r3_k1(2:end));
% 		A_inf(2,2) = 2*pi*rho*( tmp1 + tmp2 + tmp3);
% 
% 		tmp1 = .5*Rp^2*F_r3_k1(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k1(2:end);
% 		tmp2 = sum(N_gamma_l.^(-.5).*S_0l_int.*B_r3_k1) + sum(N_gamma_l.^(-.5).*S_0l_tild_int.*C_r3_k1);
% 		tmp3 = T_0n_int(1)*D_r3_k1(1) + sqrt(2)*sum(T_0n_int(2:end).*D_r3_k1(2:end));
% 		tmp4 = T_0n_tild_int(1)*E_r3_k1(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*E_r3_k1(2:end));
% 		tmp5 = sum(N_gamma_l.^(-.5).*N_0l_int.*G_r3_k1) + sum(N_gamma_l.^(-.5).*N_0l_tild_int.*H_r3_k1);
% 		A_inf(2,5) = 2*pi*rho*( tmp1 - (tmp2 + tmp3 + tmp4 + tmp5 - (Rbo^4 - Rbi^4)/(16*(e1-b)) ) );
% 	end
%%% ***********************************************************************    
%%%		k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************
% if w<=length(Omega)
    %%% solution particulière pour r = Rc
    P_3n_k2(1,1) = - ((e1-b)/12)*(2-3*( Rbo^2/(e1-b)^2)) - (.5*(b^2-e1^2)/(e1-b) + (g/omega^2));
    P_3n_k2(2:Nn,1) = -(sqrt(2)*(e1-b))./(n.*pi).^2 - sqrt(2).*f2(alpha_n,-e1,-b,0,e1)./(e1-b);

    Q_3n_k2(1,1) = - ((e1-b)/12)*(2-3*( Rbi^2/(e1-b)^2)) - (.5*(b^2-e1^2)/(e1-b) + (g/omega^2));
    Q_3n_k2(2:Nn,1) = -(sqrt(2)*(e1-b))./(n.*pi).^2 - sqrt(2).*f2(alpha_n,-e1,-b,0,e1)./(e1-b);

    P_3j_k2(1,1) = ( (h-e2) / 12 ) * ( 2 - 3*(Rp/(h-e2))^2);
    P_3j_k2(2:Nj,1) = (sqrt(2)*(h-e2).*(-1).^j) ./ (j.*pi).^2;

    P_3l_k2 = ( I5 + (g/omega^2).*N_gamma_l.^(-.5).*sin(gamma_l.*e1)./gamma_l ) ./ e1;    

	H_q3_k2 = zeros(Ni+3*Nl,1);
    H_q3_k2(1:Ni) = -((h-e2)/h).*(L_jtau'*(Gamma_0j.*P_3j_k2)) - (e1/h).*(K_ltau'*(S_0l_prim_Rp.*P_3l_k2)) - (.5*Rp/h).*L_jtau(1,:)';
    H_q3_k2(1+Ni:Ni+Nl) = S_0l_prim_Rbo.*P_3l_k2 - ((e1-b)/e1).*M_ntau'*(T_0n_prim_Rbo.*P_3n_k2 + T_0n_tild_prim_Rbo.*Q_3n_k2) + (.5*Rbo/e1).*M_ntau(1,:)';
    H_q3_k2(1+Ni+Nl:Ni+2*Nl) = -((e1-b)/e1).*M_ntau'*(T_0n_prim_Rbi.*P_3n_k2 + T_0n_tild_prim_Rbi.*Q_3n_k2) + (.5*Rbi/e1).*M_ntau(1,:)';
	H_q3_k2(1+Ni+2*Nl:Ni+3*Nl) = zeros(Nl,1);
% else
%     %%% solution particulière pour r = Rc
%     P_3n_k2(1,1) = - ((e1-b)/12)*(2-3*( Rbo^2/(e1-b)^2)) - (.5*(b^2-e1^2)/(e1-b));
%     P_3n_k2(2:Nn,1) = -(sqrt(2)*(e1-b))./(n.*pi).^2 - sqrt(2).*f2(alpha_n,-e1,-b,0,e1)./(e1-b);
% 
%     Q_3n_k2(1,1) = - ((e1-b)/12)*(2-3*( Rbi^2/(e1-b)^2)) - (.5*(b^2-e1^2)/(e1-b));
%     Q_3n_k2(2:Nn,1) = -(sqrt(2)*(e1-b))./(n.*pi).^2 - sqrt(2).*f2(alpha_n,-e1,-b,0,e1)./(e1-b);
% 
%     P_3j_k2(1,1) = ( (h-e2) / 12 ) * ( 2 - 3*(Rp/(h-e2))^2);
%     P_3j_k2(2:Nj,1) = (sqrt(2)*(h-e2).*(-1).^j) ./ (j.*pi).^2;
% 
%     P_3l_k2 = I5 ./ e1;    
% 
% 	H_q3_k2 = zeros(Ni+3*Nl,1);
%     H_q3_k2(1:Ni) = -((h-e2)/h).*(L_jtau'*(Gamma_0j.*P_3j_k2)) - (e1/h).*(K_ltau'*(S_0l_prim_Rp.*P_3l_k2)) - (.5*Rp/h).*L_jtau(1,:)';
%     H_q3_k2(1+Ni:Ni+Nl) = S_0l_prim_Rbo.*P_3l_k2 - ((e1-b)/e1).*M_ntau'*(T_0n_prim_Rbo.*P_3n_k2 + T_0n_tild_prim_Rbo.*Q_3n_k2) + (.5*Rbo/e1).*M_ntau(1,:)';
%     H_q3_k2(1+Ni+Nl:Ni+2*Nl) = -((e1-b)/e1).*M_ntau'*(T_0n_prim_Rbi.*P_3n_k2 + T_0n_tild_prim_Rbi.*Q_3n_k2) + (.5*Rbi/e1).*M_ntau(1,:)';
% 	H_q3_k2(1+Ni+2*Nl:Ni+3*Nl) = zeros(Nl,1);
% end

    %%% Calcul des coefficients A_mi, C_mn & E_mn
    coeff = D_m0\H_q3_k2;

    A_r3_k2 = coeff(1:Ni);
    C_r3_k2 = coeff(1+Ni:Ni+Nl);
    G_r3_k2 = coeff(1+Ni+Nl:Ni+2*Nl);
	H_r3_k2 = coeff(1+Ni+2*Nl:end);
	
    B_r3_k2 = K_ltau*A_r3_k2 - P_3l_k2;
	D_r3_k2 = M_ntau*C_r3_k2 - P_3n_k2;
    E_r3_k2 = M_ntau*G_r3_k2 - Q_3n_k2;
    F_r3_k2 = L_jtau*A_r3_k2 - P_3j_k2;

% 	if w<=length(Omega)
		tmp2 = T_0n_int(1)*D_r3_k2(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*D_r3_k2(2:end));
		tmp3 = T_0n_tild_int(1)*E_r3_k2(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*E_r3_k2(2:end));
		Zr(5,2,w) = 2*pi*rho*(tmp2 + tmp3 + (Rbo^4-Rbi^4)/(16*(e1-b)) );

		tmp1 = .5*Rp^2*F_r3_k2(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k2(2:end);
		tmp2 = .5*(h-e2)*( .5*Rp^2 - .125*(Rp^4/(h-e2)^2) );
		tmp3 = sum(N_gamma_l.^(-.5).*S_0l_int.*B_r3_k2) + sum(N_gamma_l.^(-.5).*S_0l_tild_int.*C_r3_k2);
		tmp4 = .5*(g/omega^2 - e1)*(Rp^2-Rbo^2);
		tmp5 = T_0n_int(1)*D_r3_k2(1) + sqrt(2)*sum(T_0n_int(2:end).*D_r3_k2(2:end));
		tmp6 = T_0n_tild_int(1)*E_r3_k2(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*E_r3_k2(2:end));
		tmp7 = - .5*(e1-b)*( .5*(Rbo^2-Rbi^2) - .125*((Rbo^4-Rbi^4)/(e1-b)^2) );
		tmp8 = sum(N_gamma_l.^(-.5).*N_0l_int.*G_r3_k2) + sum(N_gamma_l.^(-.5).*N_0l_tild_int.*H_r3_k2);
		tmp9 = .5*(g/omega^2 - e1)*(Rbi^2-Rc^2);
		Zr(5,5,w) = 2*pi*rho*(tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9));

		Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*A_r3_k1(1)) / (cosh(k_e*h)*besselh(0,k_e*Rp));
		Fe_Haskind(w,5) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*A_r3_k2(1)) / (cosh(k_e*h)*besselh(0,k_e*Rp));
% 	else
% 		tmp2 = T_0n_int(1)*D_r3_k2(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*D_r3_k2(2:end));
% 		tmp3 = T_0n_tild_int(1)*E_r3_k2(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*E_r3_k2(2:end));
% 		A_inf(5,2) = 2*pi*rho*(tmp2 + tmp3 + (Rbo^4-Rbi^4)/(16*(e1-b)) );
% 
% 		tmp1 = .5*Rp^2*F_r3_k2(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k2(2:end);
% 		tmp2 = .5*(h-e2)*( .5*Rp^2 - .125*(Rp^4/(h-e2)^2) );
% 		tmp3 = sum(N_gamma_l.^(-.5).*S_0l_int.*B_r3_k2) + sum(N_gamma_l.^(-.5).*S_0l_tild_int.*C_r3_k2);
% 		tmp4 = .5*(-e1)*(Rp^2-Rbo^2);
% 		tmp5 = T_0n_int(1)*D_r3_k2(1) + sqrt(2)*sum(T_0n_int(2:end).*D_r3_k2(2:end));
% 		tmp6 = T_0n_tild_int(1)*E_r3_k2(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*E_r3_k2(2:end));
% 		tmp7 = - .5*(e1-b)*( .5*(Rbo^2-Rbi^2) - .125*((Rbo^4-Rbi^4)/(e1-b)^2) );
% 		tmp8 = sum(N_gamma_l.^(-.5).*N_0l_int.*G_r3_k2) + sum(N_gamma_l.^(-.5).*N_0l_tild_int.*H_r3_k2);
% 		tmp9 = .5*(-e1)*(Rbi^2-Rc^2);
% 		A_inf(5,5) = 2*pi*rho*(tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9));
% 	end
% end% *** END OF HEAVE

%%% ***********************************************************************
%%%                         PITCH MODE
% if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
%%% ***********************************************************************
%%% ***********************************************************************
%%%     k = 1       Buoy oscillating, plate fixed
%%% *********************************************************************** 

    P_5n_k1(1,1) = (1/(24*(e1-b)))*(3*Rbo^3 - 4*Rbo*(e1-b)^2);
    P_5n_k1(2:Nn,1) = -sqrt(2)*Rbo*(e1-b).*(-1).^n./(n.*pi).^2;

    Q_5n_k1(1,1) = (1/(24*(e1-b)))*(3*Rbi^3 - 4*Rbi*(e1-b)^2);
    Q_5n_k1(2:Nn,1) = -sqrt(2)*Rbi*(e1-b).*(-1).^n./(n.*pi).^2;
	
	H_q5_k1 = zeros(Ni+3*Nl,1);
    H_q5_k1(1:Ni) = zeros(Ni,1);
    H_q5_k1(1+Ni:Ni+Nl) = -((e1-b)/e1).*M_ntau'*(T_1n_prim_Rbo.*P_5n_k1 + T_1n_tild_prim_Rbo.*Q_5n_k1) - (1/(2*e1*(e1-b))).*I7_k1 + ((3*Rbo^2)/(8*e1)).*M_ntau(1,:)' + (1/e1).*I3;
    H_q5_k1(1+Ni+Nl:Ni+2*Nl) = -((e1-b)/e1).*M_ntau'*(T_1n_prim_Rbi.*P_5n_k1 + T_1n_tild_prim_Rbi.*Q_5n_k1) - (1/(2*e1*(e1-b))).*I7_k1 + ((3*Rbi^2)/(8*e1)).*M_ntau(1,:)' + (1/e1).*I3;
	H_q5_k1(1+Ni+2*Nl:Ni+3*Nl) = zeros(Nl,1);

    coeff = D_m1\H_q5_k1;

    A_r5_k1 = coeff(1:Ni);
    C_r5_k1 = coeff(1+Ni:Ni+Nl);
    G_r5_k1 = coeff(1+Ni+Nl:Ni+2*Nl);
	H_r5_k1 = coeff(1+Ni+2*Nl:end);

    B_r5_k1 = K_ltau*A_r5_k1;
	D_r5_k1 = M_ntau*C_r5_k1 - P_5n_k1;
    E_r5_k1 = M_ntau*G_r5_k1 - Q_5n_k1;
    F_r5_k1 = L_jtau*A_r5_k1;

% 	if w<=length(Omega)
		tmp1 =  ( .125/(e1-b) ) * ( (Rbo^6-Rbi^6)/6 - (Rbo^4-Rbi^4)*(e1-b)^2);
		tmp2 = T_1n_int(1)*D_r5_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r5_k1(2:end));
		tmp3 = T_1n_tild_int(1)*E_r5_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r5_k1(2:end));
		Zr(3,3,w) = - pi*rho*( Rbo*I3'*C_r5_k1 - Rbi*I3'*G_r5_k1 + tmp1 + tmp2 + tmp3);

		tmp1 = .25*F_r5_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k1(2:end);
		tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r5_k1) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r5_k1);
		tmp3 = T_1n_int(1)*D_r5_k1(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r5_k1(2:end));
		tmp4 = T_1n_tild_int(1)*E_r5_k1(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r5_k1(2:end));
		tmp5 = (Rbo^6-Rbi^6)/(48*(e1-b));
		tmp6 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r5_k1) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r5_k1);

		Zr(3,6,w) = - pi*rho*(Rp*I4*A_r5_k1 + Rc*I6'*H_r5_k1 + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5 + tmp6));

		Zr(3,1,w) = -pi*rho*(Rbo*I1'*C_r5_k1 - Rbi*I1'*G_r5_k1);
		Zr(3,4,w) = -pi*rho*( Rp*I2*A_r5_k1 + Rc*I10'*H_r5_k1 );
% 	else
% 		tmp1 =  ( .125/(e1-b) ) * ( (Rbo^6-Rbi^6)/6 - (Rbo^4-Rbi^4)*(e1-b)^2);
% 		tmp2 = T_1n_int(1)*D_r5_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r5_k1(2:end));
% 		tmp3 = T_1n_tild_int(1)*E_r5_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r5_k1(2:end));
% 		A_inf(3,3) = - pi*rho*( Rbo*I3'*C_r5_k1 - Rbi*I3'*G_r5_k1 + tmp1 + tmp2 + tmp3);
% 
% 		tmp1 = .25*F_r5_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k1(2:end);
% 		tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r5_k1) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r5_k1);
% 		tmp3 = T_1n_int(1)*D_r5_k1(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r5_k1(2:end));
% 		tmp4 = T_1n_tild_int(1)*E_r5_k1(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r5_k1(2:end));
% 		tmp5 = (Rbo^6-Rbi^6)/(48*(e1-b));
% 		tmp6 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r5_k1) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r5_k1);
% 
% 		A_inf(3,6) = - pi*rho*(Rp*I4*A_r5_k1 + Rc*I6'*H_r5_k1 + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5 + tmp6));
% 
% 		A_inf(3,1) = -pi*rho*(Rbo*I1'*C_r5_k1 - Rbi*I1'*G_r5_k1);
% 		A_inf(3,4) = -pi*rho*( Rp*I2*A_r5_k1 + Rc*I10'*H_r5_k1 );
% 	end

%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed, plate oscillating 
%%% ***********************************************************************

%     P_5n_k2(1,1) = -(1/(24*(e1-b))) * (3*Rbo^3 - 4*Rbo*(e1-b)^2) + (Rbo/(e1-b))*(.5*(b^2-e1^2) + (g/omega^2)*(e1-b));
%     P_5n_k2(2:Nn,1) = 2*sqrt(2)*Rbo*(e1-b).*(-1).^n./(n.*pi).^2;% + (sqrt(2)*Rb/(e1-b))*f2(alpha_n,-e1,-b,0,e1);
% 
%     Q_5n_k2(1,1) = -(1/(24*(e1-b)))*(3*Rbi^3 - 4*Rbi*(e1-b)^2) + (Rbi/(e1-b))*(.5*(b^2-e1^2) + (g/omega^2)*(e1-b));
%     Q_5n_k2(2:Nn,1) = 2*sqrt(2)*Rbi*(e1-b).*(-1).^n./(n.*pi).^2;
% if w<=length(Omega)
    p_05 = -(Rbo^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rbo^2));
    p_n5 = ((sqrt(2)*Rbo)/(2*(e1-b)^2)).*f1(alpha_n,-e1,-b,b,e1);
    P_n5 = [p_05;p_n5];
    
    s_05 = - (Rbo/(e1-b))*(.5*(b^2-e1^2) + (g/omega^2)*(e1-b));
    s_n5 = - (sqrt(2)*Rbo/(e1-b))*f2(alpha_n,-e1,-b,0,e1);
    S_n5 = [s_05;s_n5];

	P_5n_k2 = [P_n5-S_n5];
	
    p2_05 = -(Rbi^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rbi^2));
    p2_n5 = ((sqrt(2)*Rbi)/(2*(e1-b)^2)).*f1(alpha_n,-e1,-b,b,e1);
    P2_n5 = [p2_05;p2_n5];
    
    s2_05 = - (Rbi/(e1-b))*(.5*(b^2-e1^2) + (g/omega^2)*(e1-b));
    s2_n5 = - (sqrt(2)*Rbi/(e1-b))*f2(alpha_n,-e1,-b,0,e1);
    S2_n5 = [s2_05;s2_n5];

	Q_5n_k2 = [P2_n5-S2_n5];
	

    P_5j_k2(1,1) = (1/(24*(h-e2))) * ( 3*Rp^3 - 4*Rp*(h-e2)^2);
    P_5j_k2(2:Nj,1) = -(sqrt(2)*Rp*(h-e2).*(-1).^j) ./ (j.*pi).^2;
    
    P_5l_k2 =  - (Rp/e1)*N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,g/omega^2,e1);

	H_q5_k2 = zeros(Ni+3*Nl,1);
	%%% Interface r = Rp
	tmp1 = N_lambda_i.^(-.5).*f2(lambda_i,-e1,0,g/omega^2,h)';
    H_q5_k2(1:Ni) =  -((h-e2)/h).*(L_jtau'*(Gamma_1j.*P_5j_k2)) - (e1/h).*(K_ltau'*(S_1l_prim_Rp.*P_5l_k2))...
				+ ((3*Rp^2)/(8*h)).*L_jtau(1,:)' - (.5/(h*(h-e2))).*I9' - (1/h).*tmp1' + (1/h).*I4';
    
	tmp1 = N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,g/omega^2,e1);
    H_q5_k2(1+Ni:Ni+Nl) = S_1l_prim_Rbo.*P_5l_k2 - ((e1-b)/e1).*M_ntau'*(T_1n_prim_Rbo.*P_5n_k2 + T_1n_tild_prim_Rbo.*Q_5n_k2) + (1/e1).*tmp1...
							 - ((3*Rbo^2)/(8*e1)).*M_ntau(1,:)' + (.5/(e1*(e1-b))).*I7_k2;
    H_q5_k2(1+Ni+Nl:Ni+2*Nl) = - ((e1-b)/e1).*M_ntau'*(T_1n_prim_Rbi.*P_5n_k2 + T_1n_tild_prim_Rbi.*Q_5n_k2) + (1/e1).*tmp1...
							 - ((3*Rbi^2)/(8*e1)).*M_ntau(1,:)' + (.5/(e1*(e1-b))).*I7_k2;
	H_q5_k2(1+Ni+2*Nl:Ni+3*Nl) = (1/e1)*I6 + (1/e1)*N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,g/omega^2,e1);
% else
%     p_05 = -(Rbo^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rbo^2));
%     p_n5 = ((sqrt(2)*Rbo)/(2*(e1-b)^2)).*f1(alpha_n,-e1,-b,b,e1);
%     P_n5 = [p_05;p_n5];
%     
%     s_05 = - (Rbo/(e1-b))*(.5*(b^2-e1^2));
%     s_n5 = - (sqrt(2)*Rbo/(e1-b))*f2(alpha_n,-e1,-b,0,e1);
%     S_n5 = [s_05;s_n5];
% 
% 	P_5n_k2 = [P_n5-S_n5];
% 	
%     p2_05 = -(Rbi^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rbi^2));
%     p2_n5 = ((sqrt(2)*Rbi)/(2*(e1-b)^2)).*f1(alpha_n,-e1,-b,b,e1);
%     P2_n5 = [p2_05;p2_n5];
%     
%     s2_05 = - (Rbi/(e1-b))*(.5*(b^2-e1^2));
%     s2_n5 = - (sqrt(2)*Rbi/(e1-b))*f2(alpha_n,-e1,-b,0,e1);
%     S2_n5 = [s2_05;s2_n5];
% 
% 	Q_5n_k2 = [P2_n5-S2_n5];
% 	
% 
%     P_5j_k2(1,1) = (1/(24*(h-e2))) * ( 3*Rp^3 - 4*Rp*(h-e2)^2);
%     P_5j_k2(2:Nj,1) = -(sqrt(2)*Rp*(h-e2).*(-1).^j) ./ (j.*pi).^2;
%     
%     P_5l_k2 =  - (Rp/e1)*N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,0,e1);
% 
% 	H_q5_k2 = zeros(Ni+3*Nl,1);
% 	%%% Interface r = Rp
% 	tmp1 = N_lambda_i.^(-.5).*f2(lambda_i,-e1,0,0,h)';
%     H_q5_k2(1:Ni) =  -((h-e2)/h).*(L_jtau'*(Gamma_1j.*P_5j_k2)) - (e1/h).*(K_ltau'*(S_1l_prim_Rp.*P_5l_k2))...
% 				+ ((3*Rp^2)/(8*h)).*L_jtau(1,:)' - (.5/(h*(h-e2))).*I9' - (1/h).*tmp1' + (1/h).*I4';
%     
% 	tmp1 = N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,0,e1);
%     H_q5_k2(1+Ni:Ni+Nl) = S_1l_prim_Rbo.*P_5l_k2 - ((e1-b)/e1).*M_ntau'*(T_1n_prim_Rbo.*P_5n_k2 + T_1n_tild_prim_Rbo.*Q_5n_k2) + (1/e1).*tmp1...
% 							 - ((3*Rbo^2)/(8*e1)).*M_ntau(1,:)' + (.5/(e1*(e1-b))).*I7_k2;
%     H_q5_k2(1+Ni+Nl:Ni+2*Nl) = - ((e1-b)/e1).*M_ntau'*(T_1n_prim_Rbi.*P_5n_k2 + T_1n_tild_prim_Rbi.*Q_5n_k2) + (1/e1).*tmp1...
% 							 - ((3*Rbi^2)/(8*e1)).*M_ntau(1,:)' + (.5/(e1*(e1-b))).*I7_k2;
% 	H_q5_k2(1+Ni+2*Nl:Ni+3*Nl) = (1/e1)*I6 + (1/e1)*N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,0,e1);
% end
    coeff = D_m1\H_q5_k2;

    A_r5_k2 = coeff(1:Ni);
    C_r5_k2 = coeff(1+Ni:Ni+Nl);
    G_r5_k2 = coeff(1+Ni+Nl:Ni+2*Nl);
	H_r5_k2 = coeff(1+Ni+2*Nl:end);

    B_r5_k2 = K_ltau*A_r5_k2 - P_5l_k2;
	D_r5_k2 = M_ntau*C_r5_k2 - P_5n_k2;
    E_r5_k2 = M_ntau*G_r5_k2 - Q_5n_k2;
    F_r5_k2 = L_jtau*A_r5_k2 - P_5j_k2;

% 	if w<=length(Omega)
		tmp1 = (-b^3/3 + .5*(g/omega^2 - zG(1))*b^2 + b*zG(1)*(g/omega^2));
		tmp2 = T_1n_int(1)*D_r5_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r5_k2(2:end));
		tmp3 = T_1n_tild_int(1)*E_r5_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r5_k2(2:end));
		tmp4 = - (Rbo^6-Rbi^6) / (48*(e1-b));
		Zr(6,3,w) = - pi*rho*( Rbo*(I3'*C_r5_k2 + Rbo*tmp1) - Rbi*(I3'*G_r5_k2 + Rbi*tmp1) + tmp2 + tmp3 + tmp4);

		%%%
		tmp0 = Rc^2*e1*(-(1/3)*e1^2 + .5*e1*(g/omega^2-zG(2)) + zG(2)*g/omega^2);
		tmp1 = .25*F_r5_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k2(2:end);
		tmp2 = ( .125/(h-e2) ) * ( Rp^6/6 - Rp^4*(h-e2)^2);
		tmp3 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r5_k2) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r5_k2);
		tmp4 = -.25*(g/omega^2 - e1)*(Rp^4-Rbo^4);

		tmp5 = T_1n_int(1)*D_r5_k2(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r5_k2(2:end));
		tmp6 = T_1n_tild_int(1)*E_r5_k2(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r5_k2(2:end));
		tmp7 = - ( .125/(e1-b) ) * ( (Rbo^6-Rbi^6)/6 - (Rbo^4-Rbi^4)*(e1-b)^2 );
		tmp8 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r5_k2) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r5_k2);
		tmp9 = -.25*(g/omega^2 - e1)*(Rbi^4-Rc^4);
		Zr(6,6,w) = - pi*rho*( Rp*I4*A_r5_k2 + Rc*I6'*H_r5_k2 + tmp0 + tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9));

		tmp1 = (.5*b^2 - (g/omega^2)*b);
		Zr(6,1,w) = - pi*rho*(Rbo*(I1'*C_r5_k2 + Rbo*tmp1) - Rbi*(I1'*G_r5_k2 + Rbi*tmp1));

		tmp1 = Rc^2*(.5*e1^2 - (g/omega^2)*e1);
		Zr(6,4,w) = - pi*rho*(Rp*I2*A_r5_k2 + Rc*I10'*H_r5_k2 + tmp1);

		Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*A_r5_k1(1)) / (cosh(k_e*h)*besselh(1,k_e*Rp));
		Fe_Haskind(w,6) = (4*rho*g*h*sqrt(N_lambda_i(1))*A_r5_k2(1)) / (cosh(k_e*h)*besselh(1,k_e*Rp));
% 	else
% 		tmp1 = (-b^3/3 + .5*(-zG(1))*b^2);
% 		tmp2 = T_1n_int(1)*D_r5_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r5_k2(2:end));
% 		tmp3 = T_1n_tild_int(1)*E_r5_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r5_k2(2:end));
% 		tmp4 = - (Rbo^6-Rbi^6) / (48*(e1-b));
% 		A_inf(6,3) = - pi*rho*( Rbo*(I3'*C_r5_k2 + Rbo*tmp1) - Rbi*(I3'*G_r5_k2 + Rbi*tmp1) + tmp2 + tmp3 + tmp4);
% 
% 		%%%
% 		tmp0 = Rc^2*e1*(-(1/3)*e1^2 + .5*e1*(-zG(2)));
% 		tmp1 = .25*F_r5_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k2(2:end);
% 		tmp2 = ( .125/(h-e2) ) * ( Rp^6/6 - Rp^4*(h-e2)^2);
% 		tmp3 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r5_k2) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r5_k2);
% 		tmp4 = -.25*(-e1)*(Rp^4-Rbo^4);
% 
% 		tmp5 = T_1n_int(1)*D_r5_k2(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r5_k2(2:end));
% 		tmp6 = T_1n_tild_int(1)*E_r5_k2(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r5_k2(2:end));
% 		tmp7 = - ( .125/(e1-b) ) * ( (Rbo^6-Rbi^6)/6 - (Rbo^4-Rbi^4)*(e1-b)^2 );
% 		tmp8 = sum(N_gamma_l.^(-.5).*N_1l_int.*G_r5_k2) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*H_r5_k2);
% 		tmp9 = -.25*(-e1)*(Rbi^4-Rc^4);
% 		A_inf(6,6) = - pi*rho*( Rp*I4*A_r5_k2 + Rc*I6'*H_r5_k2 + tmp0 + tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8 + tmp9));
% 
% 		tmp1 = (.5*b^2);
% 		A_inf(6,1) = - pi*rho*(Rbo*(I1'*C_r5_k2 + Rbo*tmp1) - Rbi*(I1'*G_r5_k2 + Rbi*tmp1));
% 
% 		tmp1 = Rc^2*(.5*e1^2);
% 		A_inf(6,4) = - pi*rho*(Rp*I2*A_r5_k2 + Rc*I10'*H_r5_k2 + tmp1);
% 	end
% end% *** END OF PITCH

%%% ***********************************************************************
%%% ***********************************************************************
%%%             Problème de DIFFRACTION
%%% ***********************************************************************
%%% ***********************************************************************
% if w<=length(Omega)
a0 = 1i*omega;

B0 = (-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k_e*h));
B1 = 2*1i*B0;

% if strcmp(options.dof,'all') || strcmp(options.dof,'heave')
        
	H_q7_m0 = zeros(Ni+3*Nl,1);
	H_q7_m0(1:Ni) = B0.*( besselj(0,k_e*Rp).*( ((h-e2)/h).*(L_jtau'*(Gamma_0j.*L_jtau(:,1)))...
			+ (e1/h).*(K_ltau'*(S_0l_prim_Rp.*K_ltau(:,1)))) + k_e*besselj(1,k_e*Rp)*[1;zeros(Ni-1,1)]);
	H_q7_m0(1+Ni:Ni+Nl) = -B0*besselj(0,k_e*Rp).*(S_0l_prim_Rbo.*K_ltau(:,1));
	H_q7_m0(1+Ni+Nl:Ni+2*Nl) = zeros(Nl,1);
	H_q7_m0(1+Ni+2*Nl:end) = zeros(Nl,1);

	coeff = D_m0\H_q7_m0;

    a_0i = coeff(1:Ni);
    c_0l = coeff(1+Ni:Ni+Nl);
    g_0l = coeff(1+Ni+Nl:Ni+2*Nl);
	h_0l = coeff(1+Ni+2*Nl:end);
   
    b_0l = K_ltau*a_0i + B0.*besselj(0,k_e*Rp).*K_ltau(:,1);
	d_0n = M_ntau*c_0l;
    e_0n = M_ntau*g_0l;
    f_0j = L_jtau*a_0i + B0.*besselj(0,k_e*Rp).*L_jtau(:,1);

	tmp1 = T_0n_int(1)*d_0n(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*d_0n(2:end));
	tmp2 = T_0n_tild_int(1)*e_0n(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*e_0n(2:end));
	Fe(w,2) = 2*pi*a0*rho*( tmp1 + tmp2 );

	tmp1 = .5*Rp^2*f_0j(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*f_0j(2:end);
	tmp2 = sum(N_gamma_l.^(-.5).*S_0l_int.*b_0l) + sum(N_gamma_l.^(-.5).*S_0l_tild_int.*c_0l);
	tmp3 = T_0n_int(1)*d_0n(1) + sqrt(2)*sum(T_0n_int(2:end).*d_0n(2:end));
	tmp4 = T_0n_tild_int(1)*e_0n(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*e_0n(2:end));
	tmp5 = sum(N_gamma_l.^(-.5).*N_0l_int.*g_0l) + sum(N_gamma_l.^(-.5).*N_0l_tild_int.*h_0l);
	Fe(w,5) = 2*pi*a0*rho*(tmp1 - (tmp2 + tmp3 + tmp4 + tmp5));

	%%% FROUDE-KRILOV
	Fe_FK(w,2) = 2*pi*a0*rho*B0*N_lambda_i(1)^(-.5)*cosh(k_e*(h-b))*(Rbo*besselj(1,k_e*Rbo) - Rbi*besselj(1,k_e*Rbi))/k_e;
	Fe_FK(w,5) = 2*pi*a0*rho*B0*N_lambda_i(1)^(-.5)*(cosh(k_e*(h-e2))*(Rp*besselj(1,k_e*Rp)) - cosh(k_e*(h-e1))*(Rp*besselj(1,k_e*Rp) - Rc*besselj(1,k_e*Rc)))/k_e;


% %         
% end
% %     
% if ~strcmp(options.dof,'heave')
	H_q7_m1 = zeros(Ni+3*Nl,1);
	H_q7_m1(1:Ni) = B1.*( besselj(1,k_e*Rp).*( ((h-e2)/h).*(L_jtau'*(Gamma_1j.*L_jtau(:,1)))...
			+ (e1/h).*(K_ltau'*(S_1l_prim_Rp.*K_ltau(:,1)))) - (k_e*besselj(0,k_e*Rp)-besselj(1,k_e*Rp)/Rp)*[1;zeros(Ni-1,1)]);
	H_q7_m1(1+Ni:Ni+Nl) = -B1*besselj(1,k_e*Rp).*(S_1l_prim_Rbo.*K_ltau(:,1));
	H_q7_m1(1+Ni+Nl:Ni+2*Nl) = zeros(Nl,1);
	H_q7_m1(1+Ni+2*Nl:end) = zeros(Nl,1);

	coeff = D_m1\H_q7_m1;

    a_1i = coeff(1:Ni);
    c_1l = coeff(1+Ni:Ni+Nl);
    g_1l = coeff(1+Ni+Nl:Ni+2*Nl);
	h_1l = coeff(1+Ni+2*Nl:end);
   
    b_1l = K_ltau*a_1i + B1.*besselj(1,k_e*Rp).*K_ltau(:,1);
	d_1n = M_ntau*c_1l;
    e_1n = M_ntau*g_1l;
    f_1j = L_jtau*a_1i + B1.*besselj(1,k_e*Rp).*L_jtau(:,1);
% end
% %     
% if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
	
    Fe(w,1) = - pi*a0*rho*(Rbo*I1'*c_1l - Rbi*I1'*g_1l);
	Fe(w,4) = - pi*a0*rho*(Rp*(B1*besselj(1,k_e*Rp)*I2(1) + I2*a_1i) + Rc*I10'*h_1l);
		
	tmp = N_lambda_i(1)^(-.5)*(sinh(k_e*h) - sinh(k_e*(h-b)))/k_e;			%%% Int_(-b)^(0){ Z_i(z) dz }
	Fe_FK(w,1) = -pi*a0*rho*B1*(Rbo*besselj(1,k_e*Rbo)*tmp - Rbi*besselj(1,k_e*Rbi)*tmp);
	Fe_FK(w,4) = -pi*a0*rho*B1*(Rp*besselj(1,k_e*Rp)*I2(1) + Rc*besselj(1,k_e*Rc)*N_lambda_i(1)^(-.5)*(sinh(k_e*h) - sinh(k_e*(h-e1)))/k_e);
% end
    
% if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')

    tmp1 = T_1n_int(1)*d_1n(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*d_1n(2:end));
    tmp2 = T_1n_tild_int(1)*e_1n(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*e_1n(2:end));
	Fe(w,3) = -pi*a0*rho*( Rbo*I3'*c_1l - Rbi*I3'*g_1l + tmp1 + tmp2);

	tmp1 = .25*f_1j(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*f_1j(2:end);
	tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*b_1l) + sum(N_gamma_l.^(-.5).*S_1l_tild_int.*c_1l);
	tmp3 = T_1n_int(1)*d_1n(1) + sqrt(2)*sum(T_1n_int(2:end).*d_1n(2:end));
	tmp4 = T_1n_tild_int(1)*e_1n(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*e_1n(2:end));
	tmp5 = sum(N_gamma_l.^(-.5).*N_1l_int.*g_1l) + sum(N_gamma_l.^(-.5).*N_1l_tild_int.*h_1l);
	Fe(w,6)  = - pi*a0*rho*( Rp*(B1*besselj(1,k_e*Rp)*I4(1) + I4*a_1i) + Rc*I6'*h_1l + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5));

	tmp = N_lambda_i(1)^(-.5)*f2(lambda_i(1),-b,0,-zG(1),h);%% int_(-b)^(0){ (z-Zc)Z_i(z) dz}
	tmp1 = Rbo*besselj(1,k_e*Rbo)*tmp - Rbi*besselj(1,k_e*Rbi)*tmp;
	tmp2 = N_lambda_i(1)^(-.5)*cosh(k_e*(h-b))*(Rbo^2*besselj(2,k_e*Rbo) - Rbi^2*besselj(2,k_e*Rbi))/k_e;
	Fe_FK(w,3) = -pi*a0*rho*B1*(tmp1 + tmp2);
	
	tmp1 = Rp*besselj(1,k_e*Rp)*N_lambda_i(1)^(-.5)*f2(lambda_i(1),-e2,-e1,-zG(2),h);
	tmp2 = Rc*besselj(1,k_e*Rc)*N_lambda_i(1)^(-.5)*f2(lambda_i(1),-e1,0,-zG(2),h);
	tmp3 = N_lambda_i(1)^(-.5)*cos(lambda_i(1)*(h-e2))*(Rp^2*besselj(2,k_e*Rp))/k_e;
	tmp4 = N_lambda_i(1)^(-.5)*cos(lambda_i(1)*(h-e1))*(Rp^2*besselj(2,k_e*Rp)-Rc^2*besselj(2,k_e*Rc))/k_e;
	Fe_FK(w,6) = -pi*a0*rho*B1*(tmp1 + tmp2 + tmp3 - tmp4);
% end
% end

end% END OF FOR() LOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%								EXIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** ATTENTION -> ici on reprend une dépendance en temps positive exp(-iwt) --> exp(iwt)
for i=1:6
	for j=1:6
		A(i,j,:)=real(Zr(i,j,1:end-1));
		B(i,j,:)=squeeze(imag(Zr(i,j,1:end-1))).*Omega;
	end
end
A_inf=real(Zr(:,:,end));
Fe=conj(Fe(1:end-1,:));
Fe_Haskind=conj(Fe_Haskind(1:end-1,:));
Fe_FK=conj(Fe_FK(1:end-1,:));

% save([DIR,filesep,'DATA',filesep,'hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');
% % *** ATTENTION -> ici on reprend une dépendance en temps positive exp(-iwt) --> exp(iwt)
% for i=1:6
% 	for j=1:6
% 		A(i,j,:) = real(Zr(i,j,:));
% 		B(i,j,:) = squeeze(imag(Zr(i,j,:))).*Omega;
% 	end
% end
% A_inf = real(A_inf);
% Fe = conj(Fe);
% Fe_Haskind = conj(Fe_Haskind);
% Fe_FK = conj(Fe_FK);

% save([DIR,filesep,'DATA',filesep,'hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');

if isfield(options,'IRF')
	if options.IRF==1	
		[IRF,t]=IRFfcn(Omega,B,25,.1);
		save([DIR,filesep,'DATA',filesep,'IRF.mat'],'IRF','t');
	end
end

end% END OF FUNCTION

% function out = f1( alpha, a, b, c, d )
% %%% Int_(a)^(b){(z+c)^2*cos(alpha*(z+d))}
% epsilon = 1e-8;
% for i=1:length(alpha)
% 	if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
%         out(i,1) = ((b+c)^3 - (a+c)^3) / 3;
% 	else
% 		out(i,1) = (((b+c)^2 - 2./alpha(i).^2).*sin(alpha(i).*(b+d)) - ((a+c)^2 - 2./alpha(i).^2).*sin((i).*(a+d)))./alpha(i)...
% 					+ 2.*((b+c).*cos(alpha(i).*(b+d)) - (a+c).*cos(alpha(i).*(a+d)))./alpha(i).^2;
% 	end
% end
% end
%% Additional functions
function out = f1(alpha, a, b, c, d )
%%% Int_(a)^(b){(z+c)^2*cos(alpha*(z+d))}
[n,m]=size(alpha);
out=zeros(n,m);
out(1:end)=((b+c)^2-2./alpha.^2).*sin(alpha.*(b+d))./alpha + 2*(b+c).*cos(alpha.*(b+d))./alpha.^2-...
	      (((a+c)^2-2./alpha.^2).*sin(alpha.*(a+d))./alpha + 2*(a+c).*cos(alpha.*(a+d))./alpha.^2);
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

function [S_prim,S_tild_prim] = Sprim_func(m,r,a,b,alpha)
	cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);

	S_prim = (alpha./cst).*( besselk(m,alpha.*b).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r))+...
							  besseli(m,alpha.*b).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));

	S_tild_prim = -(alpha./cst).*( besselk(m,alpha.*a).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r))+...
							  besseli(m,alpha.*a).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));
end

function [S_int,S_tild_int] = Sint_func(m,a,b,alpha)
        cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);
        
        S_int = ( besselk(m,alpha.*b).*(a^(m+1).*besseli(m+1,alpha.*a)-b^(m+1).*besseli(m+1,alpha.*b)) +...
				  besseli(m,alpha.*b).*(a^(m+1).*besselk(m+1,alpha.*a)-b^(m+1).*besselk(m+1,alpha.*b)) )./(alpha.*cst);
		  
        S_tild_int = -( besselk(m,alpha.*a).*(a^(m+1).*besseli(m+1,alpha.*a)-b^(m+1).*besseli(m+1,alpha.*b)) +...
					    besseli(m,alpha.*a).*(a^(m+1).*besselk(m+1,alpha.*a)-b^(m+1).*besselk(m+1,alpha.*b)) )./(alpha.*cst);    
end

function [T_prim,T_tild_prim] = Tprim_func(m,r,a,b,alpha)
	T_prim = zeros(length(alpha)+1,1);
	T_tild_prim = zeros(length(alpha)+1,1);

	switch m
		case 0
			T_prim(1,1) = (r*log(a/b))^(-1);
			T_tild_prim(1,1) = -(r*log(a/b))^(-1);

		case 1
			T_prim(1,1) = ((1/b + b/r^2)/(a/b-b/a));
			T_tild_prim(1,1) = ((-a/r^2 - 1/a)/(a/b-b/a));
	end
	cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);
	T_prim(2:end,1) = (alpha./cst).*( besselk(m,alpha.*b).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r)) +...
								besseli(m,alpha.*b).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));

	T_tild_prim(2:end,1) = -(alpha./cst).*( besselk(m,alpha.*a).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r)) +...
								besseli(m,alpha.*a).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));
end
