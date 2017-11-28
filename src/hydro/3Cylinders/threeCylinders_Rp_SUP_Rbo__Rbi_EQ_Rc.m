function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = threeCylinders_Rp_SUP_Rbo__Rbi_EQ_Rc(Omega, depth, WECstructure, options)
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
%           |-----------------------> (Rbo)
%           |-----> (Rbi=Rc)
%           |-----> (Rc)
%           .                                                              
% --------  +-----++----------------+ --------------------------- (z=0)
%           .ooooo||cccccccccccccccc|  ^      ^      ^        ^
%           |ooooo||cccccccccccccccc|  |      |      |        |
%           .ooooo||cccccccccccccccc|  |      |      |        |
%           |ooooo|+----------------+  v (b)  |      |        |
%           .ooooo|                           |      |        |
%           |ooooo|                           |      |        |
%           .ooooo|                           |      |        |
%           |ooooo|                           |      |        |
%           .ooooo|                           |      |        |
%           |ooooo|                           |      |        |
%           .ooooo|                           |      |        |
%           |ooooo|                           |      |        |
%           .ooooo|                           |      |        |
%           |ooooo|                           |      |        |
%           .ooooo|                           |      |        |
%           |ooooo|                           |      |        |
%           .ooooo|                           v (e1) |        |
%           |ooooo+--------------------------+       |        |
%           .ooooo|oooooooooooooooooooooooooo|       |        |
%           |ooooo|oooooooooooooooooooooooooo|       |        |
%           +-----+--------------------------+       v (e2)   .
%           .                                                 .
%           |--------------------------------> (Rp)           .
%                                                             |
%                                                             |
%                                                             |
%															  |
%                                                             v (h)
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

Rbo = WECstructure.Rbo;
Rc = WECstructure.Rc;
Rp = WECstructure.Rp;

if Rp<=Rbo
	error('Error! Rp must be strictly greater than Rbo');
	return;
end

b = WECstructure.b;
e1 = WECstructure.e1;
e2 = WECstructure.e2;

if ~isfield(options,'Trunc')
	Ni=80; Nl=80; Nn=100; Nj=100;
else
	Ni = options.Trunc.Ni;
	Nn = options.Trunc.Nn;
	Nj = options.Trunc.Nj;
	Nl = options.Trunc.Nl;
end


if ~isfield(options,'Zc')
	Zc=[0;0];
else
	Zc=options.Zc;
end

epsilon = 1e-7;

%%% Initialisation des sorties
Zr = zeros(6,6,length(Omega)); %%% matrice d'impédance --> A + (1i/omega)*B
A = zeros(6,6,length(Omega));
B = zeros(6,6,length(Omega));
Fe = zeros(length(Omega),6); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),6); %%% effort d'excitation
Fe_FK = zeros(length(Omega),6);%%% Froude-Krilov	

Hq1=zeros(Ni+Nl+Nn,2);
Hq3=zeros(Ni+Nl+Nn,2);
Hq5=zeros(Ni+Nl+Nn,2);
Hq7 = zeros(Ni+Nl+Nn,2);

Dm0=zeros(Ni+Nl+Nn,Ni+Nl+Nn); Dm1=zeros(Ni+Nl+Nn,Ni+Nl+Nn);

n = (1:1:Nn-1)'; lambda_n = (pi/(e1-b)).*n;
j = (1:1:Nj-1)'; lambda_j = (pi/(h-e2)).*j;

for w=1:length(Omega)+1
	clc;
	disp([num2str(round(w*100/(length(Omega)+1))),'%']);
	
	if w<=length(Omega)
		omega = Omega(w);
		lambda_i=fcn_waveDispersion(Omega(w),h,Ni);
		gamma_l=fcn_waveDispersion(Omega(w),e1,Nl); gamma_l=conj(gamma_l)';
		
		k0 = -imag(lambda_i(1));
	else
		omega = Inf;% omega >> 1
		i=1:1:Ni;
		l=(1:1:Nl)';
		lambda_i=.5*(2*i-1)*pi/h;
		gamma_l=.5*(2*l-1)*pi/e1;
	end
	
    N_lambda_i = .5.*( 1 + sin(2.*lambda_i.*h)./(2.*lambda_i.*h));
    N_gamma_l = .5.*( 1 + sin(2.*gamma_l.*e1)./(2.*gamma_l.*e1));
	
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
		    if (lambda_n(p)/gamma_l(i) < 1+epsilon)&&(lambda_n(p)/gamma_l(i) > 1-epsilon)
		        M_ntau(p+1,i) = .5*sqrt(2).*N_gamma_l(i)^(-.5);
		    else
		        M_ntau(p+1,i) = .5*sqrt(2)*N_gamma_l(i)^(-.5)*...
					( sin((lambda_n(p)-gamma_l(i))*(e1-b))/(lambda_n(p)-gamma_l(i))...
					+ sin((lambda_n(p)+gamma_l(i))*(e1-b))/(lambda_n(p)+gamma_l(i)) ) / (e1-b);
		    end
		end
	end

	Gamma_0j(1,1) = 0;
	Gamma_0j(2:Nj,1) = lambda_j.*besseli(1,lambda_j.*Rp)./besseli(0,lambda_j.*Rp);

	Gamma_1j(1,1) = 1/Rp;
	Gamma_1j(2:Nj,1) = lambda_j.*besseli(2,lambda_j.*Rp)./besseli(1,lambda_j.*Rp) + 1/Rp;

	Delta_0i(1,1:Ni) = -lambda_i.*(besselk(1,lambda_i.*Rp) ./ besselk(0,lambda_i.*Rp));
	Delta_1i(1,1:Ni) = -lambda_i.*(besselk(2,lambda_i.*Rp)./besselk(1,lambda_i.*Rp)) + 1/Rp;

    [Sp_0l_Rp,S_0l_tild_prim_Rp] = Sprim_func(0,Rp,Rp,Rbo,gamma_l);
    [S_0l_prim_Rb,S_0l_tild_prim_Rb] = Sprim_func(0,Rbo,Rp,Rbo,gamma_l);
    
    [S_1l_prim_Rp,S_1l_tild_prim_Rp] = Sprim_func(1,Rp,Rp,Rbo,gamma_l);
    [S_1l_prim_Rb,S_1l_tild_prim_Rb] = Sprim_func(1,Rbo,Rp,Rbo,gamma_l);
    
    [S_0l_int,S_0l_tild_int] = Si_func(0,Rp,Rbo,gamma_l);
    [S_1l_int,S_1l_tild_int] = Si_func(1,Rp,Rbo,gamma_l);

    [T_0n_prim_Rb,T_0n_tild_prim_Rb] = Tp_func(0,Rbo,Rc,Rbo,[0;lambda_n]);	
    [T_0n_prim_Rc,T_0n_tild_prim_Rc] = Tp_func(0,Rc,Rc,Rbo,[0;lambda_n]);
  
    [T_1n_prim_Rb,T_1n_tild_prim_Rb] = Tp_func(1,Rbo,Rc,Rbo,[0;lambda_n]);
    [T_1n_prim_Rc,T_1n_tild_prim_Rc] = Tp_func(1,Rc,Rc,Rbo,[0;lambda_n]);
	
    [Ti1_0n,Ti2_0n]=Ti_func(0,Rc,Rbo,[0;lambda_n]);
    [Ti1_1n,Ti2_1n]=Ti_func(1,Rc,Rbo,[0;lambda_n]);
	
    I1 = N_gamma_l.^(-.5).*( sin(gamma_l.*e1)-sin(gamma_l.*(e1-b)) ) ./ gamma_l;					%% Int_(-b)^(0){ Z_l^II }
    I2 = N_lambda_i.^(-.5).*( sin(lambda_i.*(h-e1)) - sin(lambda_i.*(h-e2)) ) ./ lambda_i;			%% Int_(-e2)^(-e1){ Z_i^I }
    I3 = N_gamma_l.^(-.5).*f2(gamma_l,-b,0,-Zc(1),e1);												%% Int_(-b)^(0){ (z-Zc)*Z_l^II }
    I4 = N_lambda_i.^(-.5).*f2(lambda_i,-e2,-e1,-Zc(2),h)';											%% Int_(-e2)^(-e1){ (z-Zc)*Z_i^I }
    I5 = N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,0,e1);													%% Int_(-e1)^(0){ z*Z_l^II }
    I6 = [.5*((b+Zc(2))^2-(e1+Zc(2))^2);sqrt(2).*f2(lambda_n,-e1,-b,-Zc(2),e1)];					%% Int_(-e1)^(-b){ (z-Zc)*Z_n^III }
    I7_k1 = N_gamma_l.^(-.5).*f1(gamma_l,-e1,-b,e1,e1);												%% Int_(-e1)^(-b){ (z+e1)^2*Z_l^II }
    I7_k2 = N_gamma_l.^(-.5).*f1(gamma_l,-e1,-b,b,e1);												%% Int_(-e1)^(-b){ (z+b)^2*Z_l^II }
    I8_k1 = [(e1-b)^3/3;sqrt(2).*f1(lambda_n,-e1,-b,e1,e1)];										%% Int_(-e1)^(-b){ (z+e1)^2*Z_n^III }
    I8_k2 = [(e1-b)^3/3;sqrt(2).*f1(lambda_n,-e1,-b,b,e1)];											%% Int_(-e1)^(-b){ (z+b)^2*Z_n^III }
    I9 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-e2,h,h);												%% Int_(-h)^(-e2){ (z+h)^2*Z_i^I }

	%%% *******************************************************************
    %%% Définition des matrices intrinsèques à la structure
	%%% *******************************************************************
    d_r3_tau_11 = diag(Delta_0i) - ((h-e2)/h).*(L_jtau'*((Gamma_0j*ones(1,Ni)).*L_jtau)) - (e1/h).*(K_ltau'*(diag(Sp_0l_Rp)*K_ltau));
    d_r3_tau_12 = -(e1/h).*( diag(S_0l_tild_prim_Rp)*K_ltau )';
    d_r3_tau_13 = zeros(Ni,Nn);

    d_r3_tau_21 = diag(S_0l_prim_Rb)*K_ltau;
    d_r3_tau_22 = diag(S_0l_tild_prim_Rb) - ((e1-b)/e1).*(M_ntau'*(diag(T_0n_prim_Rb)*M_ntau));
    d_r3_tau_23 = -((e1-b)/e1).*(diag(T_0n_tild_prim_Rb)*M_ntau)';

    d_r3_tau_31 = zeros(Nn,Ni);
    d_r3_tau_32 = diag(T_0n_prim_Rc)*M_ntau;
    d_r3_tau_33 = diag(T_0n_tild_prim_Rc);
    
    Dm0 = [d_r3_tau_11 d_r3_tau_12 d_r3_tau_13;...
                d_r3_tau_21 d_r3_tau_22 d_r3_tau_23;...
                d_r3_tau_31 d_r3_tau_32 d_r3_tau_33];
    % ***
    d_r1_tau_11 = diag(Delta_1i) - ((h-e2)/h).*(L_jtau'*(diag(Gamma_1j)*L_jtau)) - (e1/h).*(K_ltau'*(diag(S_1l_prim_Rp)*K_ltau));
    d_r1_tau_12 = -(e1/h).*( diag(S_1l_tild_prim_Rp)*K_ltau )';
    d_r1_tau_13 = zeros(Ni,Nn);

    d_r1_tau_21 = diag(S_1l_prim_Rb)*K_ltau;
    d_r1_tau_22 = diag(S_1l_tild_prim_Rb) - ((e1-b)/e1) .*(M_ntau'*(diag(T_1n_prim_Rb)*M_ntau));
    d_r1_tau_23 = -((e1-b)/e1).*(diag(T_1n_tild_prim_Rb)*M_ntau)';

    d_r1_tau_31 = zeros(Nn,Ni);
    d_r1_tau_32 = diag(T_1n_prim_Rc)*M_ntau;
    d_r1_tau_33 = diag(T_1n_tild_prim_Rc);
    
    Dm1 = [d_r1_tau_11 d_r1_tau_12 d_r1_tau_13;...
                d_r1_tau_21 d_r1_tau_22 d_r1_tau_23;...
                d_r1_tau_31 d_r1_tau_32 d_r1_tau_33];
			

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%                  Problème de RADIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%                         SURGE MODE
%%% ***********************************************************************                  
%%%     k = 1       Buoy oscillating, column fixed
%%% ***********************************************************************
	Hq1=zeros(Ni+Nl+Nn,2);
    Hq1(1:Ni,1) = zeros(Ni,1);
    Hq1(Ni+1:Ni+Nl,1) = (1/e1).*I1;
    Hq1(Ni+Nl+1:end,1) = zeros(Nn,1);
    
    coeff = Dm1\Hq1(:,1);

    Aq1(:,1) = coeff(1:Ni);
    Cq1(:,1) = coeff(Ni+1:Ni+Nl);
    Eq1(:,1) = coeff(Ni+Nl+1:end);
   
    Bq1(:,1) = K_ltau*Aq1(:,1);
    Dq1(:,1) = M_ntau*Cq1(:,1);
    Fq1(:,1) = L_jtau*Aq1(:,1);


    Zr(1,1,w) = -pi*rho*Rbo*(I1'*Cq1(:,1));

    Zr(4,1,w) = -pi*rho*(Rp*I2*Aq1(:,1) + Rc*Eq1(1,1)*(e1-b));

    tmp1 = Ti1_1n(1)*Dq1(1,1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq1(2:end,1));
    tmp2 = Ti2_1n(1)*Eq1(1,1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq1(2:end,1));
    Zr(3,1,w) = - pi*rho*(Rbo*(I3'*Cq1(:,1)) + tmp1 + tmp2);

    tmp1 = .25*Fq1(1,1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq1(2:end,1);
    tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*Bq1(:,1));
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*Cq1(:,1));
    tmp4 = Ti1_1n(1)*Dq1(1,1) + sqrt(2)*sum(Ti1_1n(2:end).*Dq1(2:end,1));
    tmp5 = Ti2_1n(1)*Eq1(1,1) + sqrt(2)*sum(Ti2_1n(2:end).*Eq1(2:end,1));
    Zr(6,1,w) = - pi*rho*( Rp*I4*Aq1(:,1) + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5) + Rc*I6'*Eq1(:,1) );

	
%%% ***********************************************************************
%%%     k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************

    Hq1(1:Ni,2) = (1/h).*I2';
    Hq1(Ni+1:Ni+Nl,2) = zeros(Nl,1);
    Hq1(Ni+Nl+1:end,2) = [1;zeros(Nn-1,1)];


    coeff = Dm1\Hq1(:,2);

    Aq1(:,2) = coeff(1:Ni);
    Cq1(:,2) = coeff(Ni+1:Ni+Nl);
    Eq1(:,2) = coeff(Ni+Nl+1:end);
   
    Bq1(:,2) = K_ltau*Aq1(:,2);
    Dq1(:,2) = M_ntau*Cq1(:,2);
    Fq1(:,2) = L_jtau*Aq1(:,2);   
    
    Zr(1,4,w) = -pi*rho*Rbo*(I1'*Cq1(:,2));
	
    Zr(4,4,w) = -pi*rho*(Rp*I2*Aq1(:,2) + Rc*Eq1(1,2)*(e1-b));

    tmp1 = Ti1_1n(1)*Dq1(1,2) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq1(2:end,2));
    tmp2 = Ti2_1n(1)*Eq1(1,2) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq1(2:end,2));
    Zr(3,4,w) = - pi*rho*(Rbo*(I3'*Cq1(:,2)) + tmp1 + tmp2);

    tmp1 = .25*Fq1(1,2)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq1(2:end,2);
    tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*Bq1(:,2));
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*Cq1(:,2));
    tmp4 = Ti1_1n(1)*Dq1(1,2) + sqrt(2)*sum(Ti1_1n(2:end).*Dq1(2:end,2));
    tmp5 = Ti2_1n(1)*Eq1(1,2) + sqrt(2)*sum(Ti2_1n(2:end).*Eq1(2:end,2));
    Zr(6,4,w) = - pi*rho*( Rp*I4*Aq1(:,2) + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5) + Rc*I6'*Eq1(:,2));

	Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1,1)) / (cosh(k0*h)*besselh(1,k0*Rp));
	Fe_Haskind(w,4) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1,2)) / (cosh(k0*h)*besselh(1,k0*Rp));


%%                         HEAVE MODE
%%% ***********************************************************************
%%%		k = 1	-->	Buoy oscillating, column fixed
%%% ***********************************************************************
    % % solution particulière pour r = R
    P_3n_k1(1,1) = ( (e1-b) / 12 ) * ( 2 - 3*(Rbo^2 / (e1-b)^2) );
    P_3n_k1(2:Nn,1) = (sqrt(2)*(e1-b).*(-1).^n) ./ (n.*pi).^2;
 
	Hq3 = zeros(Ni+Nl+Nn,2);
    Hq3(1:Ni,1) = zeros(Ni,1);
    Hq3(1+Ni:Ni+Nl,1) = -((e1-b)/e1) .* M_ntau'*(T_0n_prim_Rb.*P_3n_k1) - (.5*Rbo/e1).*M_ntau(1,:)';
    Hq3(1+Ni+Nl:end,1) = T_0n_prim_Rc.*P_3n_k1 + [(.5*Rc/(e1-b));zeros(Nn-1,1)];

    coeff = Dm0\Hq3(:,1);

    Aq3(:,1) = coeff(1:Ni);
    Cq3(:,1) = coeff(Ni+1:Ni+Nl);
    Eq3(:,1) = coeff(Ni+Nl+1:end);
   
    Bq3(:,1) = K_ltau*Aq3(:,1);
    Dq3(:,1) = M_ntau*Cq3(:,1) - P_3n_k1;
    Fq3(:,1) = L_jtau*Aq3(:,1);

    tmp1 = .25*(e1-b)*( (Rbo^2-Rc^2) - (1/4)*( (Rbo^4-Rc^4)/(e1-b)^2 ));                
    tmp2 = Ti1_0n(1)*Dq3(1,1) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Dq3(2:end,1));
    tmp3 = Ti2_0n(1)*Eq3(1,1) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Eq3(2:end,1));
    Zr(2,2,w) = 2*pi*rho*( tmp1 + tmp2 + tmp3);

    tmp1 = .5*Rp^2*Fq3(1,1) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*Fq3(2:end,1);
    tmp2 = sum(N_gamma_l.^(-.5).*S_0l_int.*Bq3(:,1));
    tmp3 = sum(N_gamma_l.^(-.5).*S_0l_tild_int.*Cq3(:,1));
    tmp4 = Ti1_0n(1)*Dq3(1,1) + sqrt(2)*sum(Ti1_0n(2:end).*Dq3(2:end,1));
    tmp5 = Ti2_0n(1)*Eq3(1,1) + sqrt(2)*sum(Ti2_0n(2:end).*Eq3(2:end,1));
    Zr(5,2,w) = 2*pi*rho*( tmp1 - (tmp2 + tmp3 + tmp4 + tmp5 - (Rbo^4 - Rc^4)/(16*(e1-b)) ) );
	
%%% ***********************************************************************    
%%%		k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************    
    %%% solution particulière pour r = Rc
    P_3n_k2(1,1) = - ((e1-b)/12)*(2-3*( Rbo^2/(e1-b)^2)) - (.5*(b^2-e1^2)/(e1-b) + (g/omega^2));
    P_3n_k2(2:Nn,1) = -(sqrt(2)*(e1-b))./(n.*pi).^2 - sqrt(2).*f2(lambda_n,-e1,-b,0,e1)./(e1-b);
    
    P_3j_k2(1,1) = ( (h-e2) / 12 ) * ( 2 - 3*(Rp/(h-e2))^2);
    P_3j_k2(2:Nj,1) = (sqrt(2)*(h-e2).*(-1).^j) ./ (j.^2 .* pi^2);

    P_3l_k2 = ( I5 + (g/omega^2).*N_gamma_l.^(-.5).*sin(gamma_l.*e1)./gamma_l ) ./ e1;    

    Hq3(1:Ni,2) = -((h-e2)/h).*(L_jtau'*(Gamma_0j.*P_3j_k2)) - (e1/h).*(K_ltau'*(Sp_0l_Rp.*P_3l_k2)) - (.5*Rp/h).*L_jtau(1,:)';
    Hq3(1+Ni:Ni+Nl,2) = S_0l_prim_Rb.*P_3l_k2 - ((e1-b)/e1).*M_ntau'*(T_0n_prim_Rb.*P_3n_k2) + (.5*Rbo/e1).*M_ntau(1,:)';
    Hq3(1+Ni+Nl:end,2) = T_0n_prim_Rc.*P_3n_k2 - [(.5*Rc/(e1-b));zeros(Nn-1,1)];
    
    %%% Calcul des coefficients A_mi, C_mn & E_mn
    coeff = Dm0\Hq3(:,2);

    Aq3(:,2) = coeff(1:Ni);
    Cq3(:,2) = coeff(Ni+1:Ni+Nl);
    Eq3(:,2) = coeff(Ni+Nl+1:end);

    Bq3(:,2) = K_ltau*Aq3(:,2) - P_3l_k2;
    Dq3(:,2) = M_ntau*Cq3(:,2) - P_3n_k2;
    Fq3(:,2) = L_jtau*Aq3(:,2) - P_3j_k2;

    tmp2 = Ti1_0n(1)*Dq3(1,2) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Dq3(2:end,2));
    tmp3 = Ti2_0n(1)*Eq3(1,2) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Eq3(2:end,2));
    Zr(2,5,w) = 2*pi*rho*(tmp2 + tmp3 + (Rbo^4-Rc^4)/(16*(e1-b)) );

    tmp1 = .5*Rp^2*Fq3(1,2) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*Fq3(2:end,2);
    tmp2 = .5*(h-e2)*( .5*Rp^2 - .125*(Rp^4/(h-e2)^2) );
    tmp3 = sum(N_gamma_l.^(-.5).*S_0l_int.*Bq3(:,2));
    tmp4 = sum(N_gamma_l.^(-.5).*S_0l_tild_int.*Cq3(:,2));
    tmp5 = .5*(g/omega^2 - e1)*(Rp^2-Rbo^2);
    tmp6 = Ti1_0n(1)*Dq3(1,2) + sqrt(2)*sum(Ti1_0n(2:end).*Dq3(2:end,2));
    tmp7 = Ti2_0n(1)*Eq3(1,2) + sqrt(2)*sum(Ti2_0n(2:end).*Eq3(2:end,2));
    tmp8 = - .5*(e1-b)*( .5*(Rbo^2-Rc^2) - .125*((Rbo^4-Rc^4)/(e1-b)^2) );
    Zr(5,5,w) = 2*pi*rho*( tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8) );
	
	Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1,1)) / (cosh(k0*h)*besselh(0,k0*Rp));
	Fe_Haskind(w,5) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1,2)) / (cosh(k0*h)*besselh(0,k0*Rp));


%%                         PITCH MODE
%%% ***********************************************************************
%%%     k = 1       Buoy oscillating, plate fixed
%%% *********************************************************************** 

    P_5n_k1(1,1) = (1/(24*(e1-b)))*(3*Rbo^3 - 4*Rbo*(e1-b)^2);
    P_5n_k1(2:Nn,1) = -sqrt(2)*Rbo*(e1-b).*(-1).^n./(n.*pi).^2;
    
	Hq5 = zeros(Ni+Nl+Nn,2);
    Hq5(1:Ni,1) = zeros(Ni,1);
    Hq5(Ni+1:Ni+Nl,1) = -((e1-b)/e1).*M_ntau'*(T_1n_prim_Rb.*P_5n_k1) - (1/(2*e1*(e1-b))).*I7_k1 + ((3*Rbo^2)/(8*e1)).*M_ntau(1,:)' + (1/e1).*I3;
    Hq5(Ni+Nl+1:end,1) = (T_1n_prim_Rc.*P_5n_k1) - [(3*Rc^2)/(8*(e1-b));zeros(Nn-1,1)] + (1/(2*(e1-b)^2)).*I8_k1;
    
    coeff = Dm1\Hq5(:,1);

    Aq5(:,1) = coeff(1:Ni);
    Cq5(:,1) = coeff(Ni+1:Ni+Nl);
    Eq5(:,1) = coeff(Ni+Nl+1:end);
   
    Bq5(:,1) = K_ltau*Aq5(:,1);
    Dq5(:,1) = M_ntau*Cq5(:,1) - P_5n_k1;
    Fq5(:,1) = L_jtau*Aq5(:,1);
    
    tmp1 =  ( .125/(e1-b) ) * ( (Rbo^6-Rc^6)/6 - (Rbo^4-Rc^4)*(e1-b)^2);
    tmp2 = Ti1_1n(1)*Dq5(1,1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq5(2:end,1));
    tmp3 = Ti2_1n(1)*Eq5(1,1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq5(2:end,1));
    Zr(3,3,w) = - pi*rho*(Rbo*I3'*Cq5(:,1) + tmp1 + tmp2 + tmp3);

    tmp1 = .25*Fq5(1,1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq5(2:end,1);
    tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*Bq5(:,1));
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*Cq5(:,1));
    tmp4 = Ti1_1n(1)*Dq5(1,1) + sqrt(2)*sum(Ti1_1n(2:end).*Dq5(2:end,1));
    tmp5 = Ti2_1n(1)*Eq5(1,1) + sqrt(2)*sum(Ti2_1n(2:end).*Eq5(2:end,1));
    tmp6 = (Rbo^6-Rc^6)/(48*(e1-b));
    tmp7 = .25*(b^4-e1^4) - (1/3)*(b^3-e1^3)*(2*e1-Zc(2)) + .5*(b^2-e1^2)*(e1^2-2*e1*Zc(2)) + e1^2*Zc(2)*(b-e1);
    tmp8 = Rc*(I6'*Eq5(:,1) + ( .125/(e1-b) )*(.5*Rc^3*((b+Zc(2))^2-(e1+Zc(2))^2) - 4*Rc*tmp7 ) );

    Zr(6,3,w) = - pi*rho*( Rp*I4*Aq5(:,1) + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5 + tmp6) + tmp8);
	
    Zr(1,3,w) = -pi*rho*Rbo*(I1'*Cq5(:,1));
    
    tmp1 = Rc^3/8 - Rc*(e1-b)^2/6;
    Zr(4,3,w) = -pi*rho*( Rp*I2*Aq5(:,1) + Rc*(Eq5(1,1)*(e1-b) + tmp1) );

%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed, plate oscillating 
%%% ***********************************************************************
   
    p_05 = -(Rbo^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rbo^2));
    p_n5 = ((sqrt(2)*Rbo)/(2*(e1-b)^2)).*f1(lambda_n,-e1,-b,b,e1);
    P_n5 = [p_05;p_n5];
    
    s_05 = - (Rbo/(e1-b))*(.5*(b^2-e1^2) + (g/omega^2)*(e1-b));
    s_n5 = - (sqrt(2)*Rbo/(e1-b))*f2(lambda_n,-e1,-b,0,e1);
    S_n5 = [s_05;s_n5];

	P_5n_k2 = [P_n5-S_n5];

    P_5j_k2(1,1) = (1/(24*(h-e2))) * ( 3*Rp^3 - 4*Rp*(h-e2)^2);
    P_5j_k2(2:Nj,1) = -(sqrt(2)*Rp*(h-e2).*(-1).^j) ./ (j.*pi).^2;
    
    P_5l_k2 =  - (Rp/e1)*N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,g/omega^2,e1);

	
	%%% Interface r = Rp
	tmp1 = N_lambda_i.^(-.5).*f2(lambda_i,-e1,0,g/omega^2,h)';
    Hq5(1:Ni,2) =  - ((h-e2)/h).*(L_jtau'*(Gamma_1j.*P_5j_k2)) - (e1/h).*(K_ltau'*(S_1l_prim_Rp.*P_5l_k2))...
				+ ((3*Rp^2)/(8*h)).*L_jtau(1,:)' - (.5/(h*(h-e2))).*I9' - (1/h).*tmp1' + (1/h).*I4';

	tmp1 = N_gamma_l.^(-.5).*f2(gamma_l,-e1,0,g/omega^2,e1);
    Hq5(Ni+1:Ni+Nl,2) = S_1l_prim_Rb.*P_5l_k2 - ((e1-b)/e1).*M_ntau'*(T_1n_prim_Rb.*P_5n_k2) + (1/e1).*tmp1...
							 - ((3*Rbo^2)/(8*e1)).*M_ntau(1,:)' + (.5/(e1*(e1-b))).*I7_k2;
                          
	Hq5(Ni+Nl+1:end,2) = T_1n_prim_Rc.*P_5n_k2 + [3*Rc^2/(8*(e1-b));zeros(Nn-1,1)] - (1/(2*(e1-b)^2)).*I8_k2  + (1/(e1-b)).*I6;
    
    coeff = Dm1\Hq5(:,2);

    Aq5(:,2) = coeff(1:Ni);
    Cq5(:,2) = coeff(Ni+1:Ni+Nl);
    Eq5(:,2) = coeff(Ni+Nl+1:end);
   
    Bq5(:,2) = K_ltau*Aq5(:,2) - P_5l_k2;
    Dq5(:,2) = M_ntau*Cq5(:,2) - P_5n_k2;
    Fq5(:,2) = L_jtau*Aq5(:,2) - P_5j_k2;


    tmp1 = Rbo*(-b^3/3 + .5*(g/omega^2 - Zc(1))*b^2 + b*Zc(1)*(g/omega^2));
    tmp2 = Ti1_1n(1)*Dq5(1,2) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq5(2:end,2));
    tmp3 = Ti2_1n(1)*Eq5(1,2) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq5(2:end,2));
    tmp4 = - (Rbo^6-Rc^6) / (48*(e1-b));
    Zr(3,6,w) = - pi*rho*( Rbo*(I3'*Cq5(:,2) + tmp1) + tmp2 + tmp3 + tmp4);
	
    tmp1 = .25*Fq5(1,2)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq5(2:end,2);
    tmp2 = ( .125/(h-e2) ) * ( Rp^6/6 - Rp^4*(h-e2)^2);
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_int.*Bq5(:,2));
    tmp4 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*Cq5(:,2));
    tmp5 = -.25*(g/omega^2 - e1)*(Rp^4-Rbo^4);
    tmp6 = Ti1_1n(1)*Dq5(1,2) + sqrt(2)*sum(Ti1_1n(2:end).*Dq5(2:end,2));
    tmp7 = Ti2_1n(1)*Eq5(1,2) + sqrt(2)*sum(Ti2_1n(2:end).*Eq5(2:end,2));
    tmp8 = - ( .125/(e1-b) ) * ( (Rbo^6-Rc^6)/6 - (Rbo^4-Rc^4)*(e1-b)^2 );
    
    tmp9 = -(1/3)*(-Zc(2) - e1)*(b-e1)^3 + (1/12)*(b-e1)^4;
	
    tmp10 = Rc*(I6'*Eq5(:,2) - ( .125/(e1-b) )*(.5*Rc^3*((b+Zc(2))^2-(e1+Zc(2))^2) - 4*Rc*tmp9 ) );

    Zr(6,6,w) = - pi*rho*( Rp*I4*Aq5(:,2) + tmp10 + tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7 + tmp8));
	
    tmp1 = Rbo*(.5*b^2 - (g/omega^2)*b);
    Zr(1,6,w) = - pi*rho*Rbo*(I1'*Cq5(:,2) + tmp1);
    
    tmp1 = - Rc^3/8 + Rc*(e1-b)^2/6;
    Zr(4,6,w) = - pi*rho*(Rp*I2*Aq5(:,2) + Rc*(Eq5(1,2)*(e1-b) + tmp1) );

	Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1,1)) / (cosh(k0*h)*besselh(1,k0*Rp));
	Fe_Haskind(w,6) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1,2)) / (cosh(k0*h)*besselh(1,k0*Rp));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%             Problème de DIFFRACTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a0 = 1i*omega;
%% m=0
B0 = (-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));

	Hq7 = zeros(Ni+Nl+Nn,2);
	Hq7(1:Ni,1) = B0.*( besselj(0,k0*Rp).*( ((h-e2)/h).*(L_jtau'*(Gamma_0j.*L_jtau(:,1)))...
			+ (e1/h).*(K_ltau'*(Sp_0l_Rp.*K_ltau(:,1)))) + k0*besselj(1,k0*Rp)*[1;zeros(Ni-1,1)]);
	Hq7(1+Ni:Ni+Nl,1) = -B0*besselj(0,k0*Rp).*(S_0l_prim_Rb.*K_ltau(:,1));
	Hq7(1+Ni+Nl:end,1) = zeros(Nn,1);

	coeff = Dm0\Hq7(:,1);

	Aq7(:,1) = coeff(1:Ni);
	Cq7(:,1) = coeff(Ni+1:Ni+Nl);
	Eq7(:,1) = coeff(Ni+Nl+1:end);

	Bq7(:,1) = K_ltau*Aq7(:,1) + B0.*besselj(0,k0*Rp).*K_ltau(:,1);
	Dq7(:,1) = M_ntau*Cq7(:,1);
	Fq7(:,1) = L_jtau*Aq7(:,1) + B0.*besselj(0,k0*Rp).*L_jtau(:,1);    

	tmp1 = Ti1_0n(1)*Dq7(1,1) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Dq7(2:end,1));
	tmp2 = Ti2_0n(1)*Eq7(1,1) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Eq7(2:end,1));
	Fe(w,2) = 2*pi*a0*rho*( tmp1 + tmp2 );

	tmp1 = .5*Rp^2*Fq7(1,1) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*Fq7(2:end,1);
	tmp2 = sum(N_gamma_l.^(-.5).*S_0l_int.*Bq7(:,1));
	tmp3 = sum(N_gamma_l.^(-.5).*S_0l_tild_int.*Cq7(:,1));
	tmp4 = Ti1_0n(1)*Dq7(1,1) + sqrt(2)*sum(Ti1_0n(2:end).*Dq7(2:end,1));
	tmp5 = Ti2_0n(1)*Eq7(1,1) + sqrt(2)*sum(Ti2_0n(2:end).*Eq7(2:end,1));
	Fe(w,5) = 2*pi*a0*rho*(tmp1 - (tmp2 + tmp3 + tmp4 + tmp5));

	%%% FROUDE-KRILOV
	Fe_FK(w,2) = 2*pi*a0*rho*B0*N_lambda_i(1)^(-.5)*cosh(k0*(h-b))*(Rbo*besselj(1,k0*Rbo) - Rc*besselj(1,k0*Rc))/k0;
	Fe_FK(w,5) = 2*pi*a0*rho*B0*N_lambda_i(1)^(-.5)*(cosh(k0*(h-e2))*(Rp*besselj(1,k0*Rp)) - cosh(k0*(h-e1))*(Rp*besselj(1,k0*Rp) - Rc*besselj(1,k0*Rc)))/k0;

%% m=1
B1 = 2*1i*B0;

	Hq7(1:Ni,2) = B1.*( besselj(1,k0*Rp).*( ((h-e2)/h).*(L_jtau'*(Gamma_1j.*L_jtau(:,1)))...
			+ (e1/h).*(K_ltau'*(S_1l_prim_Rp.*K_ltau(:,1)))) - (k0*besselj(0,k0*Rp)-besselj(1,k0*Rp)/Rp).*[1;zeros(Ni-1,1)]);
	Hq7(1+Ni:Ni+Nl,2) = -B1*besselj(1,k0*Rp).*(S_1l_prim_Rb.*K_ltau(:,1));
	Hq7(1+Ni+Nl:end,2) = zeros(Nn,1);

	coeff = Dm1\Hq7(:,2);

	Aq7(:,2) = coeff(1:Ni);
	Cq7(:,2) = coeff(Ni+1:Ni+Nl);
	Eq7(:,2) = coeff(Ni+Nl+1:end);

	Bq7(:,2) = K_ltau*Aq7(:,2) + B1*besselj(1,k0*Rp).*K_ltau(:,1);
	Dq7(:,2) = M_ntau*Cq7(:,2);
	Fq7(:,2) = L_jtau*Aq7(:,2) + B1*besselj(1,k0*Rp).*L_jtau(:,1);

	Fe(w,1) = -pi*a0*rho*Rbo*( I1'*Cq7(:,2) );
	Fe(w,4) = - pi*a0*rho*( Rp*(B1*besselj(1,k0*Rp)*I2(1) + I2*Aq7(:,2)) + Rc*Eq7(1,2)*(e1-b));

	tmp1 = Ti1_1n(1)*Dq7(1,2) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq7(2:end,2));
	tmp2 = Ti2_1n(1)*Eq7(1,2) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq7(2:end,2));
	Fe(w,3) = - pi*a0*rho*( Rbo*I3'*Cq7(:,2) + tmp1 + tmp2);

	tmp1 = .25*Fq7(1,2)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq7(2:end,2);
	tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*Bq7(:,2));
	tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*Cq7(:,2));
	tmp4 = Ti1_1n(1)*Dq7(1,2) + sqrt(2)*sum(Ti1_1n(2:end).*Dq7(2:end,2));
	tmp5 = Ti2_1n(1)*Eq7(1,2) + sqrt(2)*sum(Ti2_1n(2:end).*Eq7(2:end,2));

	Fe(w,6)  = - pi*a0*rho*( Rp*(B1*besselj(1,k0*Rp)*I4(1) + I4*Aq7(:,2)) + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5)  + Rc*(I6'*Eq7(:,2)));

	%%% FROUDE-KRILOV
	tmp = N_lambda_i(1)^(-.5)*(sinh(k0*h) - sinh(k0*(h-b)))/k0;			%%% Int_(-b)^(0){ Z_i(z) dz }
	Fe_FK(w,1) = -pi*a0*rho*B1*(Rbo*besselj(1,k0*Rbo)*tmp);% - Rc*besselj(1,k_e*Rc)*tmp);
	Fe_FK(w,4) = -pi*a0*rho*B1*(Rp*besselj(1,k0*Rp)*I2(1) + Rc*besselj(1,k0*Rc)*N_lambda_i(1)^(-.5)*(sinh(k0*(h-b)) - sinh(k0*(h-e1)))/k0);
	
	tmp = N_lambda_i(1)^(-.5)*f2(lambda_i(1),-b,0,-Zc(1),h);
	tmp1 = Rbo*besselj(1,k0*Rbo)*tmp;
	tmp2 = N_lambda_i(1)^(-.5)*cosh(k0*(h-b))*(Rbo^2*besselj(2,k0*Rbo) - Rc^2*besselj(2,k0*Rc))/k0;
	Fe_FK(w,3) = -pi*a0*rho*B1*(tmp1 + tmp2);

	tmp1 = Rp*besselj(1,k0*Rp)*N_lambda_i(1)^(-.5)*f2(lambda_i(1),-e2,-e1,-Zc(2),h);
	tmp2 = Rc*besselj(1,k0*Rc)*N_lambda_i(1)^(-.5)*f2(lambda_i(1),-e1,0,-Zc(2),h);
	tmp3 = N_lambda_i(1)^(-.5)*cos(lambda_i(1)*(h-e2))*(Rp^2*besselj(2,k0*Rp))/k0;
	tmp4 = N_lambda_i(1)^(-.5)*cos(lambda_i(1)*(h-e1))*(Rp^2*besselj(2,k0*Rp)-Rc^2*besselj(2,k0*Rc))/k0;
	Fe_FK(w,6) = -pi*a0*rho*B1*(tmp1 + tmp2 + tmp3 - tmp4);

end% END OF FOR() LOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%								EXIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** ATTENTION -> ici on reprend une dépendance en temps positive exp(-iwt) --> exp(iwt)
for i=1:6
	for j=1:6
		A(i,j,:) = real(Zr(i,j,1:end-1));
		B(i,j,:) = squeeze(imag(Zr(i,j,1:end-1))).*Omega;
	end
end
A_inf = real(Zr(:,:,end));
Fe = conj(Fe(1:end-1,:));
Fe_Haskind = conj(Fe_Haskind(1:end-1,:));
Fe_FK = conj(Fe_FK(1:end-1,:));

save([DIR,filesep,'DATA',filesep,'hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');

end% END OF FUNCTION

%% Additional functions
function out = f1(alpha, a, b, c, d )
% Int_(a)^(b){(z+c)^2*cos(alpha*(z+d))}
[n,m]=size(alpha);
out=zeros(n,m);
out(1:end)=((b+c)^2-2./alpha.^2).*sin(alpha.*(b+d))./alpha + 2*(b+c).*cos(alpha.*(b+d))./alpha.^2-...
	      (((a+c)^2-2./alpha.^2).*sin(alpha.*(a+d))./alpha + 2*(a+c).*cos(alpha.*(a+d))./alpha.^2);
end

function out = f2(alpha, a, b, c, d )
% Int_(a)^(b){(z+c)*cos(alpha*(z+d))}
epsilon = 1e-8;
for i=1:length(alpha)
	if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
		out(i,1) = .5*((b+c)^2 - (a+c)^2);
	else
		out(i,1) = ((b+c).*sin(alpha(i).*(b+d)) - (a+c).*sin(alpha(i).*(a+d)))./alpha(i)  + (cos(alpha(i).*(b+d)) - cos(alpha(i).*(a+d)))./alpha(i).^2;
	end
end
end

function [S1,S2] = S_func(m,r,a,b,alpha)
	cst=besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);
	S1=(besseli(m,alpha.*r).*besselk(m,alpha.*b)-besseli(m,alpha.*b).*besselk(m,alpha.*r))./cst;
	S2=(besseli(m,alpha.*a).*besselk(m,alpha.*r)-besseli(m,alpha.*r).*besselk(m,alpha.*a))./cst;
end

function [S_prim,S_tild_prim] = Sprim_func(m,r,a,b,alpha)
	cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);

	S_prim = (alpha./cst).*( besselk(m,alpha.*b).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r))+...
							  besseli(m,alpha.*b).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));

	S_tild_prim = -(alpha./cst).*( besselk(m,alpha.*a).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r))+...
							  besseli(m,alpha.*a).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));
end

function [S_int,S_tild_int] = Si_func(m,a,b,alpha)
        cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);
        
        S_int = ( besselk(m,alpha.*b).*(a^(m+1).*besseli(m+1,alpha.*a)-b^(m+1).*besseli(m+1,alpha.*b)) +...
				  besseli(m,alpha.*b).*(a^(m+1).*besselk(m+1,alpha.*a)-b^(m+1).*besselk(m+1,alpha.*b)) )./(alpha.*cst);
		  
        S_tild_int = -( besselk(m,alpha.*a).*(a^(m+1).*besseli(m+1,alpha.*a)-b^(m+1).*besseli(m+1,alpha.*b)) +...
					    besseli(m,alpha.*a).*(a^(m+1).*besselk(m+1,alpha.*a)-b^(m+1).*besselk(m+1,alpha.*b)) )./(alpha.*cst);    
end


function [T1,T2] = T_func(m,r,a,b,alpha)
[p,q]=size(alpha);
T1=zeros(p,q);
T2=zeros(p,q);

cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);
T1(:)=(besseli(m,alpha.*a).*besselk(m,alpha.*r)-besseli(m,alpha.*r).*besselk(m,alpha.*a))./cst;
T2(:)=(besseli(m,alpha.*r).*besselk(m,alpha.*b)-besseli(m,alpha.*b).*besselk(m,alpha.*r))./cst;

if abs(alpha(1))==0	
	switch m
		case 0
			T1(1) = log(r/a)./log(b/a);
			T2(1) = log(b/r)./log(b/a);
		case 1
			T1(1) = (a/r-r/a)./(a/b-b/a);
			T2(1) = (r/b-b/r)./(a/b-b/a);
	end
end
end

function [Tp1,Tp2] = Tp_func(m,r,a,b,alpha)

[p,q]=size(alpha);
Tp1=zeros(p,q);
Tp2=zeros(p,q);

cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);

Tp1(:)=-(alpha./cst).*( besseli(m,alpha.*a).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)) +...
						besselk(m,alpha.*a).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r)));

Tp2(:)=(alpha./cst).*( besselk(m,alpha.*b).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r)) +...
					   besseli(m,alpha.*b).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));				
							
if abs(alpha(1))==0	
	switch m
		case 0
			Tp1(1,1) = (r*log(b/a))^(-1);
			Tp2(1,1) = -(r*log(b/a))^(-1);
		case 1
			Tp1(1,1) = -(m/a)*((r/a)^(-m-1)+(r/a)^(m-1))/(a/b-b/a);
			Tp2(1,1) = (m/b)*((r/b)^(-m-1)+(r/b)^(m-1))/(a/b-b/a);
	end
end

end

function [Ti1,Ti2] = Ti_func(m,a,b,alpha)
[p,q]=size(alpha);
Ti1=zeros(p,q);
Ti2=zeros(p,q);

cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);

Ti1(:) = (-1./ (alpha.*cst)) .*...
		( b^(m+1).*(besseli(m,alpha.*a).*besselk(m+1,alpha.*b) + besseli(m+1,alpha.*b).*besselk(m,alpha.*a))-...
		  a^(m+1).*(besseli(m,alpha.*a).*besselk(m+1,alpha.*a) + besseli(m+1,alpha.*a).*besselk(m,alpha.*a)));

Ti2(:) = (1./ (alpha.*cst)) .*...
		( b^(m+1).*(besseli(m+1,alpha.*b).*besselk(m,alpha.*b)+besseli(m,alpha.*b).*besselk(m+1,alpha.*b))-...
		  a^(m+1).*(besseli(m+1,alpha.*a).*besselk(m,alpha.*b)+besseli(m,alpha.*b).*besselk(m+1,alpha.*a)));

if abs(alpha(1))==0	
	switch m
		case 0
			Ti1(1)=(.5*b^2*(log(b/a)-.5) - .5*a^2*(log(a/a)-.5))/log(b/a);
			Ti2(1)=-(.5*b^2*(log(b/b)-.5) - .5*a^2*(log(a/b)-.5))/log(b/a);
		case 1
			Ti1(1)=-(.25*(b^4-a^4)/a-.5*(b^2-a^2)*a)/(a/b-b/a);
			Ti2(1)= (.25*(b^4-a^4)/b-.5*(b^2-a^2)*b)/(a/b-b/a);
	end
end

end

