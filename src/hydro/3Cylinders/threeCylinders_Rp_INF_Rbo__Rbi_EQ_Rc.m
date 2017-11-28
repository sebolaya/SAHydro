function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = threeCylinders_Rp_INF_Rbo__Rbi_EQ_Rc(Omega, depth, WECstructure, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the hydrodynamic coefficients for the three cylinder
% structure depicted below (Rp<Rbo).
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
%           |ooooo+----------+   			         |        |
%           .ooooo|oooooooooo|  			         |        |
%           |ooooo|oooooooooo|  			         |        |
%           +-----+----------+  					 v (e2)   .
%           .                                                 .
%           |-------------------> (Rp)			              .
%                                                             |
%                                                             |
%                                                             |
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



g = 9.81;
rho = 1025;
h = depth; %%% water depth
    
Rbo = WECstructure.Rbo;
Rc = WECstructure.Rc;
Rp = WECstructure.Rp;

if Rp>=Rbo
	error('Error! Rp must be strictly lower than Rbo');
	return;
end

b = WECstructure.b;
e1 = WECstructure.e1;
e2 = WECstructure.e2;

if ~isfield(options,'Trunc')
	Ni=80; Nl=100; Nn=100; Nj=100;
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
A = zeros(6,6,length(Omega)); A_inf = zeros(6,6);
B = zeros(6,6,length(Omega));
Fe = zeros(length(Omega),6); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),6); %%% effort d'excitation
Fe_FK = zeros(length(Omega),6);%%% Froude-Krilov

Hq1=zeros(Ni+Nl+Nn,2);
Hq3=zeros(Ni+Nl+Nn,2);
Hq5=zeros(Ni+Nl+Nn,2);
Hq7=zeros(Ni+Nl+Nn,2);

Dm0=zeros(Ni+Nl+Nn,Ni+Nl+Nn); Dm1=zeros(Ni+Nl+Nn,Ni+Nl+Nn);

% Initialisation des coefficients de Fourier
Aq1=zeros(Ni,2); Bq1=zeros(Nl,2); Cq1=zeros(Nl,2); Dq1=zeros(Nn,2); Eq1=zeros(Nn,2); Fq1=zeros(Nj,2);
Aq3=zeros(Ni,2); Bq3=zeros(Nl,2); Cq3=zeros(Nl,2); Dq3=zeros(Nn,2); Eq3=zeros(Nn,2); Fq3=zeros(Nj,2); 
Aq5=zeros(Ni,2); Bq5=zeros(Nl,2); Cq5=zeros(Nl,2); Dq5=zeros(Nn,2); Eq5=zeros(Nn,2); Fq5=zeros(Nj,2);
Aq7=zeros(Ni,2); Bq7=zeros(Nl,2); Cq7=zeros(Nl,2); Dq7=zeros(Nn,2); Eq7=zeros(Nn,2); Fq7=zeros(Nj,2);

n = (1:1:Nn-1)'; lambda_n = (pi/(e1-b)).*n;
j = (1:1:Nj-1)'; lambda_j = (pi/(h-e2)).*j;
l = (1:1:Nl-1)'; lambda_l = (pi/(h-b)).*l;

for w=1:length(Omega)+1
clc;
disp([num2str(round(w*100/(length(Omega)+1))),'%']);

if w<=length(Omega)
	omega = Omega(w);
	lambda_i=fcn_waveDispersion(Omega(w),h,Ni);
	k0 = -imag(lambda_i(1));
else
	omega = Inf;% omega >> 1
	i=1:1:Ni;
	lambda_i=.5*(2*i-1)*pi/h;
end
    
    N_lambda_i = .5.*( 1 + sin(2.*lambda_i.*h)./(2.*lambda_i.*h));

	K_ltau = zeros(Nl,Ni);
	K_ltau(1,:)=N_lambda_i.^(-.5).*sin(lambda_i.*(h-b))./( (h-b).*lambda_i );
	for i=1:Ni
		for p=1:Nl-1
		    if (lambda_l(p)/lambda_i(i) < 1+epsilon)&&(lambda_l(p)/lambda_i(i) > 1-epsilon)
		        K_ltau(p+1,i) = .5*sqrt(2).*N_lambda_i(i)^(-.5);
		    else
		        K_ltau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*...
					( sin((lambda_l(p)-lambda_i(i))*(h-b))/(lambda_l(p)-lambda_i(i))...
					+ sin((lambda_l(p)+lambda_i(i))*(h-b))/(lambda_l(p)+lambda_i(i)))/(h-b);
		    end
		end
	end

	L_jtau=zeros(Nj,Nl);
	L_jtau(1,1)=1;
	L_jtau(1,2:end)=sqrt(2).*sin(lambda_l.*(h-e2))./((h-e2).*lambda_l);
	for i=1:Nl-1
		for p=1:Nj-1
		    if (lambda_j(p)/lambda_l(i) < 1+epsilon)&&(lambda_j(p)/lambda_l(i) > 1-epsilon)
		        L_jtau(p+1,i+1)=1+sin(2*lambda_j(i).*(h-e2))./(2*lambda_j(i).*(h-e2));
		    else
		        L_jtau(p+1,i+1)= ...
					( sin((lambda_j(p)-lambda_l(i))*(h-e2))/(lambda_j(p)-lambda_l(i))...
					+ sin((lambda_j(p)+lambda_l(i))*(h-e2))/(lambda_j(p)+lambda_l(i)))/(h-e2);
		    end
		end
	end
	
	
	M_ntau=zeros(Nn,Nl);
	M_ntau(1,1)=1;
	M_ntau(1,2:end)=sqrt(2).*sin(lambda_l.*(e1-h))./((e1-b).*lambda_l);
	for i=1:Nl-1
		for p=1:Nn-1
		    if (lambda_l(i)/lambda_n(p) < 1+epsilon)&&(lambda_l(i)/lambda_n(p) > 1-epsilon)
				M_ntau(p+1,i+1) = 2*(.25*(sin(lambda_n(p)*(e1+h-2*b))-sin(lambda_n(p)*(h-e1)))/lambda_n(p)+.5*(e1-b)*cos(lambda_n(p)*(e1-h)))/(e1-b);
			else
				M_ntau(p+1,i+1) = ...
					( (sin(-(lambda_n(p)-lambda_l(i))*b+lambda_n(p)*e1-lambda_l(i)*h)-sin(-(lambda_n(p)-lambda_l(i))*e1+lambda_n(p)*e1-lambda_l(i)*h))/(lambda_n(p)-lambda_l(i))...
					+ (sin(-(lambda_n(p)+lambda_l(i))*b+lambda_n(p)*e1+lambda_l(i)*h)-sin(-(lambda_n(p)+lambda_l(i))*e1+lambda_n(p)*e1+lambda_l(i)*h))/(lambda_n(p)+lambda_l(i)) ) / (e1-b);
		    end
		end
	end
	
	
    Gamma_0j(1,1) = 0;
    Gamma_0j(2:Nj,1) = lambda_j.*(besseli(1,lambda_j.*Rp)./besseli(0,lambda_j.*Rp));
	
	Gamma_1j(1,1) = 1/Rp;
	Gamma_1j(2:Nj,1) = lambda_j.*besseli(2,lambda_j.*Rp)./besseli(1,lambda_j.*Rp) + 1/Rp;
    
	Delta_0i(1,1:Ni) = -lambda_i.*(besselk(1,lambda_i.*Rbo)./besselk(0,lambda_i.*Rbo));
	Delta_1i(1,1:Ni) = -lambda_i.*(besselk(2,lambda_i.*Rbo)./besselk(1,lambda_i.*Rbo)) + 1/Rbo;

    [Tp1_0l_Rbo,Tp2_0l_Rbo] = Tp_func(0,Rbo,Rp,Rbo,[0;lambda_l]);
    [Tp1_0l_Rp,Tp2_0l_Rp] = Tp_func(0,Rp,Rp,Rbo,[0;lambda_l]);
    
    [Tp1_1l_Rbo,Tp2_1l_Rbo] = Tp_func(1,Rbo,Rp,Rbo,[0;lambda_l]);
    [Tp1_1l_Rp,Tp2_1l_Rp] = Tp_func(1,Rp,Rp,Rbo,[0;lambda_l]);
    
    [Ti1_0l,Ti2_0l] = Ti_func(0,Rp,Rbo,[0;lambda_l]);
    [Ti1_1l,Ti2_1l] = Ti_func(1,Rp,Rbo,[0;lambda_l]);
	
    [Tp1_0n_Rp,Tp2_0n_Rp] = Tp_func(0,Rp,Rc,Rp,[0;lambda_n]);
    [Tp1_0n_Rc,Tp2_0n_Rc] = Tp_func(0,Rc,Rc,Rp,[0;lambda_n]);

    [Tp1_1n_Rp,Tp2_1n_Rp] = Tp_func(1,Rp,Rc,Rp,[0;lambda_n]);
    [Tp1_1n_Rc,Tp2_1n_Rc] = Tp_func(1,Rc,Rc,Rp,[0;lambda_n]);   
    
    [Ti1_0n,Ti2_0n] = Ti_func(0,Rc,Rp,[0;lambda_n]);
    [Ti1_1n,Ti2_1n] = Ti_func(1,Rc,Rp,[0;lambda_n]);

    I1 = N_lambda_i.^(-.5).*( sin(lambda_i.*h)-sin(lambda_i.*(h-b)) ) ./ lambda_i;              % % Int_(-b)^(0){ Z_i^I }
    I2 = [(e2-e1);sqrt(2).*(sin(lambda_l.*(h-e1)) - sin(lambda_l.*(h-e2)))./lambda_l];             % % Int_(-e2)^(-e1){ Z_l^II }
    I3 = N_lambda_i.^(-.5).*f2(lambda_i,-b,0,0,h);                                                 % % Int_(-b)^(0){ z*Z_i^I }
    I4 = [.5*(e1^2-e2^2);sqrt(2).*f2(lambda_l,-e2,-e1,0,h)];                                         % % Int_(-e2)^(-e1){ z*Z_l^II }
    I5_1 = [(e1-b)^3/3;sqrt(2).*f1(lambda_l,-e1,-b,e1,h)];                                         % % Int_{-e1}^{-b}{(z+e1)^2*Z_l^II}
    I5_2 = [(e1-b)^3/3;sqrt(2).*f1(lambda_l,-e1,-b,b,h)];                                          % % Int_{-e1}^{-b}{(z+e1)^2*Z_l^II}
    I6 = [(h-b)^3/3;sqrt(2).*f1(lambda_l,-h,-b,h,h)];                                              % % Int_{-h}^{-b}{(z+h)^2*Z_l^II}
    I7 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-b,h,h);                                               % % Int_(-h)^(-b){ (z+h)^2*Z_i^I }
    I8 = [(h-e2)^3/3;sqrt(2).*f1(lambda_l,-h,-e2,h,h)];                                            % % Int_{-h}^{-e2}{(z+h)^2*Z_l^II}
    I9 = [.5*(b^2-e1^2);sqrt(2).*f2(lambda_n,-e1,-b,0,e1)];                                          % % Int_(-e1)^(-b){ z*Z_n^III }
    I10_1 = [(e1-b)^3/3;sqrt(2).*f1(lambda_n,-e1,-b,e1,e1)];                                        % % Int_(-e1)^(-b){ (z+e1)^2*Z_n^III }
    I10_2 = [(e1-b)^3/3;sqrt(2).*f1(lambda_n,-e1,-b,b,e1)];                                        % % Int_(-e1)^(-b){ (z+b)^2*Z_n^III }
    

%% Définition des matrices intrinsèques à la structure
    % Interface r = Rb
    Dm0(1:Ni,1:Ni) = diag(Delta_0i) - ((h-b)/h).*(K_ltau'*(diag(Tp1_0l_Rbo)*K_ltau));
    Dm0(1:Ni,Ni+1:Ni+Nl) = -((h-b)/h).*(diag(Tp2_0l_Rbo)*K_ltau)';
    Dm0(1:Ni,Ni+Nl+1:end) = zeros(Ni,Nn);
    % Interface r = Rp
    Dm0(Ni+1:Ni+Nl,1:Ni)=diag(Tp1_0l_Rp)*K_ltau ;
    Dm0(Ni+1:Ni+Nl,Ni+1:Ni+Nl)=diag(Tp2_0l_Rp) - ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_0j)*L_jtau)) - ((e1-b)/(h-b)).*(M_ntau'*(diag(Tp1_0n_Rp)*M_ntau));
    Dm0(Ni+1:Ni+Nl,Ni+Nl+1:end)=-((e1-b)/(h-b)).*(diag(Tp2_0n_Rp)*M_ntau)';
    % Interface r = Rc
    Dm0(Ni+Nl+1:end,1:Ni) = zeros(Nn,Ni);
    Dm0(Ni+Nl+1:end,Ni+1:Ni+Nl) = diag(Tp1_0n_Rc)*M_ntau;
    Dm0(Ni+Nl+1:end,Ni+Nl+1:end) = diag(Tp2_0n_Rc);
    
	% Interface r = Rb
    Dm1(1:Ni,1:Ni)=diag(Delta_1i) - ((h-b)/h).*(K_ltau'*(diag(Tp1_1l_Rbo)*K_ltau));
    Dm1(1:Ni,Ni+1:Ni+Nl)=-((h-b)/h).*(diag(Tp2_1l_Rbo)*K_ltau)';
    Dm1(1:Ni,Ni+Nl+1:end)=zeros(Ni,Nn);
    % % Interface r = Rp
    Dm1(Ni+1:Ni+Nl,1:Ni)=diag(Tp1_1l_Rp)*K_ltau ;
    Dm1(Ni+1:Ni+Nl,Ni+1:Ni+Nl)=diag(Tp2_1l_Rp) - ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_1j)*L_jtau)) - ((e1-b)/(h-b)).*(M_ntau'*(diag(Tp1_1n_Rp)*M_ntau));
    Dm1(Ni+1:Ni+Nl,Ni+Nl+1:end)=-((e1-b)/(h-b)).*(diag(Tp2_1n_Rp)*M_ntau)';
    % % Interface r = Rc
    Dm1(Ni+Nl+1:end,1:Ni) = zeros(Nn,Ni);
    Dm1(Ni+Nl+1:end,Ni+1:Ni+Nl) = diag(Tp1_1n_Rc)*M_ntau;
    Dm1(Ni+Nl+1:end,Ni+Nl+1:end) = diag(Tp2_1n_Rc);

%% Problème de RADIATION
%% SURGE MODE

%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, column fixed
%%% ***********************************************************************
Hq1(1:Ni,1)=(1/h).*I1';					% Interface r = Rb
Hq1(Ni+1:Ni+Nl,1)=zeros(Nl,1);			% Interface r = Rp
Hq1(Ni+Nl+1:end,1)=zeros(Nn,1);			% Interface r = Rc

coeff = Dm1\Hq1(:,1);

Aq1(:,1) = coeff(1:Ni);
Cq1(:,1) = coeff(Ni+1:Ni+Nl);
Eq1(:,1) = coeff(Ni+Nl+1:end);

Bq1(:,1) = K_ltau*Aq1(:,1);
Dq1(:,1) = M_ntau*Cq1(:,1);
Fq1(:,1) = L_jtau*Cq1(:,1);

Zr(1,1,w) = - pi*rho*Rbo*(I1*Aq1(:,1));
Zr(4,1,w) = - pi*rho*(Rp*I2'*Cq1(:,1) + Rc*Eq1(1,1)*(e1-b));

tmp1 = Ti1_1l(1)*Bq1(1,1) + sqrt(2)*sum((-1).^l.*Ti1_1l(2:end).*Bq1(2:end,1));
tmp2 = Ti2_1l(1)*Cq1(1,1) + sqrt(2)*sum((-1).^l.*Ti2_1l(2:end).*Cq1(2:end,1));
tmp3 = Ti1_1n(1)*Dq1(1,1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq1(2:end,1));
tmp4 = Ti2_1n(1)*Eq1(1,1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq1(2:end,1));
Zr(3,1,w) = - pi*rho*(Rbo*I3*Aq1(:,1) + tmp1 + tmp2 + tmp3 + tmp4);

tmp1 = .25*Fq1(1,1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq1(2:end,1);
tmp2 = Ti1_1n(1)*Dq1(1,1) + sqrt(2)*sum(Ti1_1n(2:end).*Dq1(2:end,1));
tmp3 = Ti2_1n(1)*Eq1(1,1) + sqrt(2)*sum(Ti2_1n(2:end).*Eq1(2:end,1));
Zr(6,1,w) = - pi*rho*( Rp*I4'*Cq1(:,1) + tmp1 - (tmp2 + tmp3) + Rc*(I9'*Eq1(:,1)) );

%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************
Hq1(1:Ni,2) = zeros(Ni,1);          % Interface r = Rb
Hq1(Ni+1:Ni+Nl,2) = (1/(h-b)).*I2;        % Interface r = Rp
Hq1(Ni+Nl+1:end,2) = [1;zeros(Nn-1,1)];       % Interface r = Rc

coeff = Dm1\Hq1(:,2);

Aq1(:,2) = coeff(1:Ni);
Cq1(:,2) = coeff(Ni+1:Ni+Nl);
Eq1(:,2) = coeff(Ni+Nl+1:end);

Bq1(:,2) = K_ltau*Aq1(:,2);
Dq1(:,2) = M_ntau*Cq1(:,2);
Fq1(:,2) = L_jtau*Cq1(:,2);

Zr(1,4,w) = - pi*rho*Rbo*(I1*Aq1(:,2));
Zr(4,4,w) = - pi*rho*(Rp*I2'*Cq1(:,2) + Rc*Eq1(1,2)*(e1-b));

tmp1 = Ti1_1l(1)*Bq1(1,2) + sqrt(2)*sum((-1).^l.*Ti1_1l(2:end).*Bq1(2:end,2));
tmp2 = Ti2_1l(1)*Cq1(1,2) + sqrt(2)*sum((-1).^l.*Ti2_1l(2:end).*Cq1(2:end,2));
tmp3 = Ti1_1n(1)*Dq1(1,2) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq1(2:end,2));
tmp4 = Ti2_1n(1)*Eq1(1,2) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq1(2:end,2));
Zr(3,4,w) = - pi*rho*(Rbo*I3*Aq1(:,2) + tmp1 + tmp2 + tmp3 + tmp4);

tmp1 = .25*Fq1(1,2)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq1(2:end,2);
tmp2 = Ti1_1n(1)*Dq1(1,2) + sqrt(2)*sum(Ti1_1n(2:end).*Dq1(2:end,2));
tmp3 = Ti2_1n(1)*Eq1(1,2) + sqrt(2)*sum(Ti2_1n(2:end).*Eq1(2:end,2));
Zr(6,4,w) = - pi*rho*( Rp*I4'*Cq1(:,2) + tmp1 - (tmp2+tmp3) + Rc*(I9'*Eq1(:,2)) );

Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1,1)) / (cosh(k0*h)*besselh(1,k0*Rbo));
Fe_Haskind(w,4) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1,2)) / (cosh(k0*h)*besselh(1,k0*Rbo));
	
%%                         HEAVE MODE
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, column fixed
%%% *********************************************************************** 
    p_03 = ( (e1-b) / 12 ) * ( 2 - 3*(Rp^2 / (e1-b)^2) );
    p_n3 = sqrt(2).*f1(lambda_n,-e1,-b,e1,e1)./(2*(e1-b)^2);
    P_n3 = [p_03;p_n3];
    
    r1_03 = ( 2*((h-b)^3-(h-e1)^3) - 3*Rp^2*(e1-b) ) / (12*(h-b)*(e1-b)) ;
    r_n3 = sqrt(2).*f1(lambda_n,-e1,-b,h,e1)./(2*(h-b)*(e1-b));
    R_n3 = [r1_03;r_n3];
    
	P3n(:,1)=R_n3-P_n3;
	
    r2_03 = ( (h-b) / 12 ) * ( 2 - 3*(Rbo^2 / (h-b)^2) );
    r_l3 = sqrt(2).*f1(lambda_l,-h,-b,h,h)./(2*(h-b)^2);
    P3l(:,1) = -[r2_03;r_l3];

    r3_03 = ( (h-e2)^2 / (12*(h-b)) ) * ( 2 - 3*(Rp^2 / (h-e2)^2 ));
    r_j3 = sqrt(2).*f1(lambda_j,-h,-e2,h,h)./(2*(h-b)*(h-e2));
    P3j(:,1) = [r3_03;r_j3];

    Hq3(1:Ni,1)=((h-b)/h).*K_ltau'*(Tp1_0l_Rbo.*P3l(:,1)) - (.5*Rbo/h).*K_ltau(1,:)';
    Hq3(Ni+1:Ni+Nl,1)=((e1-b)/(h-b)).*(M_ntau'*(Tp1_0n_Rp.*P3n(:,1))) + ((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*P3j(:,1))) - Tp1_0l_Rp.*P3l(:,1) - (.5*Rp/(h-b)).*M_ntau(1,:)' + [.5*Rp/(h-b);zeros(Nl-1,1)];
    Hq3(Ni+Nl+1:end,1) = -Tp1_0n_Rc.*P3n(:,1) + [(.5*Rc/(e1-b));zeros(Nn-1,1)];

    coeff = Dm0\Hq3(:,1);

    Aq3(:,1) = coeff(1:Ni);
    Cq3(:,1) = coeff(Ni+1:Ni+Nl);
    Eq3(:,1) = coeff(Ni+Nl+1:end);
    
    Bq3(:,1) = K_ltau*Aq3(:,1)+P3l(:,1);
    Dq3(:,1) = M_ntau*Cq3(:,1)+P3n(:,1);
    Fq3(:,1) = L_jtau*Cq3(:,1)+P3j(:,1);

    tmp1 = Ti1_0l(1)*Bq3(1,1) + sqrt(2)*sum((-1).^l.*Ti1_0l(2:end).*Bq3(2:end,1));
    tmp2 = Ti2_0l(1)*Cq3(1,1) + sqrt(2)*sum((-1).^l.*Ti2_0l(2:end).*Cq3(2:end,1));
	tmp3 = .5*(h-b)*(.5*(Rbo^2 - Rp^2) - (1/8)*( (Rbo^4-Rp^4)/(h-b)^2 ));
    tmp4 = Ti1_0n(1)*Dq3(1,1) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Dq3(2:end,1));
    tmp5 = Ti2_0n(1)*Eq3(1,1) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Eq3(2:end,1));
    tmp6 = .5*(e1-b)*(.5*(Rp^2-Rc^2) - (1/8)*( (Rp^4-Rc^4)/(e1-b)^2 ));
    Zr(2,2,w) = 2*pi*rho*(tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6);

    tmp1 = .5*Rp^2*Fq3(1,1) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*Fq3(2:end,1);
    tmp2 = Ti1_0n(1)*Dq3(1,1) + sqrt(2)*sum(Ti1_0n(2:end).*Dq3(2:end,1));
    tmp3 = Ti2_0n(1)*Eq3(1,1) + sqrt(2)*sum(Ti2_0n(2:end).*Eq3(2:end,1));
    tmp4 = - (Rp^4-Rc^4)/(16*(e1-b));
    Zr(5,2,w) = 2*pi*rho*( tmp1 - (tmp2 + tmp3 + tmp4) );
            
%%% ***********************************************************************    
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************     
    % % solution particulière pour r = Rc
    P3n(1,2)=((e1-b)/12)*(2 - 3*Rp^2/(e1-b)^2);
    P3n(2:Nn,2)=.5*sqrt(2).*f1(lambda_n,-e1,-b,b,e1)./(e1-b)^2;%(sqrt(2)*(e1-b)) ./ (n.^2 .* pi^2);
    
    P3j(1,2)= -((h-e2)/12)*(2 - 3*Rp/(h-e2)^2);
    P3j(2:Nj,2)= -.5*sqrt(2).*f1(lambda_j,-h,-e2,h,h)./(h-e2)^2;%-(sqrt(2)*(h-e2).*(-1).^j) ./ (j.^2 .* pi^2);
    
    Hq3(1:Ni,2) = zeros(Ni,1);
    Hq3(Ni+1:Ni+Nl,2) = ((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*P3j(:,2))) + ((e1-b)/(h-b)).*(M_ntau'*(Tp1_0n_Rp.*P3n(:,2))) + (.5*Rp/(h-b)).*(M_ntau(1,:)'-L_jtau(1,:)');
    Hq3(Ni+Nl+1:end,2) = -Tp1_0n_Rc.*P3n(:,2) - [.5*Rc/(e1-b);zeros(Nn-1,1)];

    coeff = Dm0\Hq3(:,2);

    Aq3(:,2) = coeff(1:Ni);
    Cq3(:,2) = coeff(Ni+1:Ni+Nl);
    Eq3(:,2) = coeff(Ni+Nl+1:end);
   
    Bq3(:,2) = K_ltau*Aq3(:,2);
    Dq3(:,2) = M_ntau*Cq3(:,2) + P3n(:,2);
    Fq3(:,2) = L_jtau*Cq3(:,2) + P3j(:,2);
    
    tmp1 = Ti1_0l(1)*Bq3(1,2) + sqrt(2)*sum((-1).^l.*Ti1_0l(2:end).*Bq3(2:end,2));
    tmp2 = Ti2_0l(1)*Cq3(1,2) + sqrt(2)*sum((-1).^l.*Ti2_0l(2:end).*Cq3(2:end,2));
    tmp3 = Ti1_0n(1)*Dq3(1,2) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Dq3(2:end,2));
    tmp4 = Ti2_0n(1)*Eq3(1,2) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Eq3(2:end,2));
    tmp5 = (Rp^4-Rc^4)/(16*(e1-b));
    Zr(2,5,w) = 2*pi*rho*( tmp1 + tmp2 + tmp3 + tmp4 + tmp5);

    tmp1 = .5*Rp^2*Fq3(1,2) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*Fq3(2:end,2);
    tmp2 = .5*(h-e2)*(.5*Rp^2 - .125*Rp^4/(h-e2)^2);
    tmp3 = Ti1_0n(1)*Dq3(1,2) + sqrt(2)*sum(Ti1_0n(2:end).*Dq3(2:end,2));
    tmp4 = Ti2_0n(1)*Eq3(1,2) + sqrt(2)*sum(Ti2_0n(2:end).*Eq3(2:end,2));
    tmp5 = - .5*(e1-b)*(.5*(Rp^2-Rc^2) - .125*(Rp^4-Rc^4)/(e1-b)^2);
    Zr(5,5,w) = 2*pi*rho*(tmp1 + tmp2 - (tmp3 + tmp4 + tmp5));

	Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1,1)) / (cosh(k0*h)*besselh(0,k0*Rbo));
	Fe_Haskind(w,5) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1,2)) / (cosh(k0*h)*besselh(0,k0*Rbo));

%% PITCH MODE
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, plate fixed
%%% *********************************************************************** 

    % % solution particulière pour r = R
    p_05 = (Rp^3/(24*(e1-b))) * ( 3 - 4*((e1-b)/Rp)^2);% Phi3
    p_n5 = -((sqrt(2)*Rp)/(2*(e1-b)^2)).*f1(lambda_n,-e1,-b,e1,e1);
    P_n5 = [p_05;p_n5];
    
    s_05 =  (Rp^3/(24*(e1-b)*(h-b)))*( 3*(e1-b) - 4*((h-b)^3-(h-e1)^3)/Rp^2 ) ;%Phi2
    s_n5 = -(sqrt(2)*Rp)/(2*(e1-b)*(h-b)).*f1(lambda_n,-e1,-b,h,e1);
    S_n5 = [s_05;s_n5];
    
    P5l(1,1) = (Rbo^3/(24*(h-b))) * ( 3 - 4*((h-b)/Rbo)^2 );
    P5l(2:Nl,1) = -((sqrt(2)*Rbo)/(2*(h-b)^2)).*f1(lambda_l,-h,-b,h,h);
	P5l=-P5l;

    P5j(1,1) = (Rp^3/(24*(h-b))) * ( 3 - 4*((h-e2)/Rp)^2);
    P5j(2:Nj,1) = -((sqrt(2)*Rp)/(2*(h-b)*(h-e2))).*f1(lambda_j,-h,-e2,h,h);

    Hq5(1:Ni,1) = ((h-b)/h).*K_ltau'*(Tp1_1l_Rbo.*P5l(:,1)) + ((3*Rbo^2)/(8*h)).*K_ltau(1,:)' - (.5/(h*(h-b))).*I7' + (1/h).*I3';
    Hq5(Ni+1:Ni+Nl,1) = ((e1-b)/(h-b)).*(M_ntau'*(Tp1_1n_Rp.*(S_n5-P_n5))) + ((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*P5j(:,1))) - Tp1_1l_Rp.*P5l(:,1)...
							- [3*Rp^2/(8*(h-b));zeros(Nl-1,1)] + (.5/(h-b)^2).*I6 + ((3*Rp^2)/(8*(h-b))).*M_ntau(1,:)' - (.5/((e1-b)*(h-b))).*I5_1;
    Hq5(Ni+Nl+1:end,1) = -Tp1_1n_Rc.*(S_n5-P_n5) - [(3*Rc^2)/(8*(e1-b));zeros(Nn-1,1)] + (.5/(e1-b)^2).*I10_1;
	
    coeff = Dm1\Hq5(:,1);

    Aq5(:,1) = coeff(1:Ni);
    Cq5(:,1) = coeff(Ni+1:Ni+Nl);
    Eq5(:,1) = coeff(Ni+Nl+1:end);
    
    Bq5(:,1) = K_ltau*Aq5(:,1) + P5l(:,1);
    Dq5(:,1) = M_ntau*Cq5(:,1) + (S_n5-P_n5);
    Fq5(:,1) = L_jtau*Cq5(:,1) + P5j(:,1);
    
    tmp1 = Ti1_1l(1)*Bq5(1,1) + sqrt(2)*sum((-1).^l.*Ti1_1l(2:end).*Bq5(2:end,1));
    tmp2 = Ti2_1l(1)*Cq5(1,1) + sqrt(2)*sum((-1).^l.*Ti2_1l(2:end).*Cq5(2:end,1));
    tmp3 = ((Rbo^6-Rp^6)/6 - (Rbo^4-Rp^4)*(h-b)^2) / (8*(h-b));
    tmp4 = Ti1_1n(1)*Dq5(1,1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq5(2:end,1));
    tmp5 = Ti2_1n(1)*Eq5(1,1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq5(2:end,1));
    tmp6 = ((Rp^6-Rc^6)/6 - (Rp^4-Rc^4)*(e1-b)^2) / (8*(e1-b));
    Zr(3,3,w) = - pi*rho*(Rbo*I3*Aq5(:,1) + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6);
    
    tmp1 = (.5*Rp^3*(e1^2-e2^2) - 4*Rp*( e2*(h-e2)^3 - e1*(h-e1)^3 - ((h-e1)^4-(h-e2)^4)/4 )/3 ) / (8*(h-b));
    tmp2 = .25*Fq5(1,1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq5(2:end,1);
    tmp3 = Ti1_1n(1)*Dq5(1,1) + sqrt(2)*sum(Ti1_1n(2:end).*Dq5(2:end,1));
    tmp4 = Ti2_1n(1)*Eq5(1,1) + sqrt(2)*sum(Ti2_1n(2:end).*Eq5(2:end,1));
    tmp5 = (Rp^6-Rc^6)/(48*(e1-b));
    tmp6 = ( .125/(e1-b) )*(.5*Rc^3*(b^2-e1^2) - 4*Rc*( -(1/12)*(e1-b)^4 - (1/3)*b*(e1-b)^3 ));
    Zr(6,3,w) = - pi*rho*( Rp*(I4'*Cq5(:,1) + tmp1) + tmp2 - (tmp3 + tmp4 + tmp5) + Rc*(I9'*Eq5(:,1) + tmp6) );
    
    Zr(1,3,w) =  - pi*rho*Rbo*(I1*Aq5(:,1));
    
    tmp1 = (Rp^3*(e2-e1) - 4*Rp*((h-e1)^3-(h-e2)^3)/3)/(8*(h-b));
    tmp2 = Rc^3/8 - Rc*(e1-b)^2/6;
    Zr(4,3,w) = - pi*rho*( Rp*(I2'*Cq5(:,1) + tmp1) + Rc*(Eq5(1,1)*(e1-b) + tmp2) );

%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% *********************************************************************** 

    p_05 = -(Rp^3/(24*(e1-b))) * ( 3 - 4*((e1-b)/Rp)^2);
    p_n5 = ((sqrt(2)*Rp)/(2*(e1-b)^2)).*f1(lambda_n,-e1,-b,b,e1);
    P_n5 = -[p_05;p_n5];

    q_05 = (Rp^3/(24*(h-e2))) * ( 3 - 4*((h-e2)/Rp)^2);
    q_j5 = -((sqrt(2)*Rp)/(2*(h-e2)^2)).*f1(lambda_j,-h,-e2,h,h);
    P5j(:,2) = -[q_05;q_j5];

    Hq5(1:Ni,2) = zeros(Ni,1);
    Hq5(Ni+1:Ni+Nl,2) = ((e1-b)/(h-b)).*(M_ntau'*(Tp1_1n_Rp.*P_n5)) - ((3*Rp^2)/(8*(h-b))).*M_ntau(1,:)' + (1/(2*(e1-b)*(h-b))).*I5_2...
               + ((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*P5j(:,2))) + ((3*Rp^2)/(8*(h-b))).*L_jtau(1,:)' - (1/(2*(h-e2)*(h-b))).*I8 + (1/(h-b)).*I4;
    Hq5(Ni+Nl+1:end,2) = -Tp1_1n_Rc.*P_n5 + [(3*Rc^2)/(8*(e1-b));zeros(Nn-1,1)] - (1/(2*(e1-b)^2)).*I10_2  + (1/(e1-b)).*I9;

    coeff = Dm1\Hq5(:,2);
    
    Aq5(:,2) = coeff(1:Ni);
    Cq5(:,2) = coeff(Ni+1:Ni+Nl);
    Eq5(:,2) = coeff(Ni+Nl+1:end);
    
    Bq5(:,2) = K_ltau*Aq5(:,2);
    Dq5(:,2) = M_ntau*Cq5(:,2) + P_n5 ;
    Fq5(:,2) = L_jtau*Cq5(:,2) + P5j(:,2);
    
    
    tmp1 = Ti1_1l(1)*Bq5(1,2) + sqrt(2)*sum((-1).^l.*Ti1_1l(2:end).*Bq5(2:end,2));
    tmp2 = Ti2_1l(1)*Cq5(1,2) + sqrt(2)*sum((-1).^l.*Ti2_1l(2:end).*Cq5(2:end,2));
    tmp3 = Ti1_1n(1)*Dq5(1,2) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq5(2:end,2));
    tmp4 = Ti2_1n(1)*Eq5(1,2) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq5(2:end,2));
    tmp5 = -((Rp^6-Rc^6)/6) / (8*(e1-b));
    Zr(3,6,w) = - pi*rho*(Rbo*I3*Aq5(:,2) + tmp1 + tmp2 + tmp3 + tmp4 + tmp5);
    
    tmp1 = .25*Fq5(1,2)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq5(2:end,2);
    tmp2 = (Rp^6/6 - Rp^4*(h-e2)^2) / (8*(h-e2));
    tmp3 = Ti1_1n(1)*Dq5(1,2) + sqrt(2)*sum(Ti1_1n(2:end).*Dq5(2:end,2));
    tmp4 = Ti2_1n(1)*Eq5(1,2) + sqrt(2)*sum(Ti2_1n(2:end).*Eq5(2:end,2));
    tmp5 = -( (Rp^6-Rc^6)/6 - (Rp^4-Rc^4)*(e1-b)^2) / (8*(e1-b));
    tmp6 = - ( .125/(e1-b) )*(.5*Rc^3*(b^2-e1^2) - 4*Rc*( (1/12)*(e1-b)^4 - (1/3)*e1*(e1-b)^3 ));
    Zr(6,6,w) = - pi*rho*( Rp*I4'*Cq5(:,2) + tmp1 + tmp2 - (tmp3 + tmp4 + tmp5) + Rc*(I9'*Eq5(:,2) + tmp6) );

    Zr(1,6,w) = - pi*rho*Rbo*(I1*Aq5(:,2));
    
    tmp1 = - Rc^3/8 + Rc*(e1-b)^2/6;
    Zr(4,6,w) = - pi*rho*(Rp*I2'*Cq5(:,2) + Rc*((e1-b)*Eq5(1,2) + tmp1));
	
	Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1,1)) / (cosh(k0*h)*besselh(1,k0*Rbo));
	Fe_Haskind(w,6) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1,2)) / (cosh(k0*h)*besselh(1,k0*Rbo));

%% Problème de DIFFRACTION
a0 = 1i*omega;
%% m=0
B0=(-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));
        
Hq7(1:Ni,1) = B0.*(besselj(0,k0*Rbo)*((h-b)/h).*(K_ltau'*(Tp1_0l_Rbo.*K_ltau(:,1))) + k0*besselj(1,k0*Rbo).*[1;zeros(Ni-1,1)]);
Hq7(Ni+1:Ni+Nl,1) = -B0*besselj(0,k0*Rbo).*(Tp1_0l_Rp.*K_ltau(:,1));
Hq7(Ni+Nl+1:end,1) = zeros(Nn,1);

coeff = Dm0\Hq7(:,1);

Aq7(:,1) = coeff(1:Ni);
Cq7(:,1) = coeff(Ni+1:Ni+Nl);
Eq7(:,1) = coeff(Ni+Nl+1:end);

Bq7(:,1) = K_ltau*Aq7(:,1)+B0.*besselj(0,k0*Rbo).*K_ltau(:,1);
Dq7(:,1) = M_ntau*Cq7(:,1);
Fq7(:,1) = L_jtau*Cq7(:,1);    

tmp1 = Ti1_0l(1)*Bq7(1,1) + sqrt(2)*sum((-1).^l.*Ti1_0l(2:end).*Bq7(2:end,1));
tmp2 = Ti2_0l(1)*Cq7(1,1) + sqrt(2)*sum((-1).^l.*Ti2_0l(2:end).*Cq7(2:end,1));
tmp3 = Ti1_0n(1)*Dq7(1,1) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Dq7(2:end,1));
tmp4 = Ti2_0n(1)*Eq7(1,1) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Eq7(2:end,1));
Fe(w,2) = 2*pi*a0*rho*(tmp1 + tmp2 + tmp3 + tmp4);

tmp1 = .5*Rp^2*Fq7(1,1) + sqrt(2)*Rp*((-1).^j.*besseli(1,lambda_j.*Rp)./(lambda_j.*besseli(0,lambda_j.*Rp)))'*Fq7(2:end,1);
tmp2 = Ti1_0n(1)*Dq7(1,1) + sqrt(2)*sum(Ti1_0n(2:end).*Dq7(2:end,1));
tmp3 = Ti2_0n(1)*Eq7(1,1) + sqrt(2)*sum(Ti2_0n(2:end).*Eq7(2:end,1));
Fe(w,5) = 2*pi*a0*rho*( tmp1 - (tmp2 + tmp3) );

%% m=1  
B1 = 2*1i*B0;

Hq7(1:Ni,2)=B1.*(besselj(1,k0*Rbo)*((h-b)/h).*(K_ltau'*(Tp1_1l_Rbo.*K_ltau(:,1))) - (k0*besselj(0,k0*Rbo)-besselj(1,k0*Rbo)/Rbo).*[1;zeros(Ni-1,1)]);
Hq7(Ni+1:Ni+Nl,2)=-B1*besselj(1,k0*Rbo).*(Tp1_1l_Rp.*K_ltau(:,1));
Hq7(Ni+Nl+1:end,2)=zeros(Nn,1);

coeff = Dm1\Hq7(:,2);

Aq7(:,2) = coeff(1:Ni);
Cq7(:,2) = coeff(Ni+1:Ni+Nl);
Eq7(:,2) = coeff(Ni+Nl+1:end);

Bq7(:,2) = K_ltau*Aq7(:,2)+B1*besselj(1,k0*Rbo).*K_ltau(:,1);
Dq7(:,2) = M_ntau*Cq7(:,2);
Fq7(:,2) = L_jtau*Cq7(:,2);

Fe(w,1) = -pi*a0*rho*Rbo*(B1*besselj(1,k0*Rbo)*I1(1) + I1*Aq7(:,2));
Fe(w,4) = -pi*a0*rho*( Rp*I2'*Cq7(:,2) + Rc*Eq7(1,2)*(e1-b) );

tmp1 = Ti1_1l(1)*Bq7(1,2) + sqrt(2)*sum((-1).^l.*Ti1_1l(2:end).*Bq7(2:end,2));
tmp2 = Ti2_1l(1)*Cq7(1,2) + sqrt(2)*sum((-1).^l.*Ti2_1l(2:end).*Cq7(2:end,2));
tmp3 = Ti1_1n(1)*Dq7(1,2) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Dq7(2:end,2));
tmp4 = Ti2_1n(1)*Eq7(1,2) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Eq7(2:end,2));
Fe(w,3) = - pi*a0*rho*( Rbo*(B1*besselj(1,k0*Rbo)*I3(1) + I3*Aq7(:,2)) + tmp1 + tmp2 + tmp3 + tmp4 );

tmp1 = .25*Fq7(1,2)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,lambda_j.*Rp))./(lambda_j.*besseli(1,lambda_j.*Rp)))'*Fq7(2:end,2);
tmp2 = Ti1_1n(1)*Dq7(1,2) + sqrt(2)*sum(Ti1_1n(2:end).*Dq7(2:end,2));
tmp3 = Ti2_1n(1)*Eq7(1,2) + sqrt(2)*sum(Ti2_1n(2:end).*Eq7(2:end,2));  
Fe(w,6) = - pi*a0*rho*( Rp*I4'*Cq7(:,2) + tmp1 - (tmp2 + tmp3) + Rc*(I9'*Eq7(:,2)) );
 
end% END OF FOR() LOOP

%% EXIT
% *** ATTENTION -> ici on reprend une dépendance en temps positive exp(-iwt) --> exp(iwt)
for i=1:6
	for j=1:6
		A(i,j,:) = real(Zr(i,j,1:end-1));
		B(i,j,:) = squeeze(imag(Zr(i,j,1:end-1))).*Omega;
	end
end
A_inf(:,:) = real(Zr(:,:,end));
Fe = conj(Fe(1:end-1,:));
Fe_Haskind = conj(Fe_Haskind(1:end-1,:));
Fe_FK = conj(Fe_FK(1:end-1,:));

save([DIR,filesep,'DATA',filesep,'hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');
end% END OF FUNCTION

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
[n,m]=size(alpha);
out=zeros(n,m);
for i=1:length(alpha)
	if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
		out(i) = .5*((b+c)^2 - (a+c)^2);
	else
		out(i) = ((b+c).*sin(alpha(i).*(b+d)) - (a+c).*sin(alpha(i).*(a+d)))./alpha(i)  + (cos(alpha(i).*(b+d)) - cos(alpha(i).*(a+d)))./alpha(i).^2;
	end
end
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
			Ti1(1)=(.5*b^2*(log(b/a)-.5) - .5*a^2*(log(a/a)-.5))/log(b/a);%(.5*b^2*log(b/a) - .25*(b^2-a^2))/log(b/a);
			Ti2(1)=-(.5*b^2*(log(b/b)-.5) - .5*a^2*(log(a/b)-.5))/log(b/a);%(-.5*a^2*log(b/a) + .25*(b^2-a^2))/log(b/a);
		case 1
			Ti1(1)=-(.25*(b^4-a^4)/a-.5*(b^2-a^2)*a)/(a/b-b/a);
			Ti2(1)= (.25*(b^4-a^4)/b-.5*(b^2-a^2)*b)/(a/b-b/a);
	end
end

end
