function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = twoCylinder_Rp_EQ_Rbo(Omega, depth, WECstructure, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the hydrodynamic coefficients for the three cylinder
% structure depicted below (Rp=Rbo).
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
%           |-----------------------> (Rbo)
%           .                                                              
% --------  +-----------------------+ --------------------------- (z=0)
%           |ccccccccccccccccccccccc|  ^      ^      ^        ^
%           |ccccccccccccccccccccccc|  |      |      |        |
%           |ccccccccccccccccccccccc|  |      |      |        |
%           +-----------------------+  v (b)  |      |        |
%			.								  |      |        |
%			|								  |      |        |
%			.								  |      |        |
%           +-----------------------+         v (e1) |        |
%           |ooooooooooooooooooooooo|                |        |
%           |ooooooooooooooooooooooo|                |        |
%           +-----------------------+                v (e2)   .
%           .                                                 .
%           |-----------------------> (Rp=Rbo)                .
%           .                                                 |
%           |                                                 |
%                                                             |
%                                                             v (h)
% 
%
% Written by S. OLAYA, IRDL/ENIB*, olaya@enib.fr
% * ENIB - École Nationale d'Ingénieurs de Brest, France
% Date : 2016/15/03
% Please include the following citation when you use this code :
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
b = WECstructure.b;
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
	Zc=[0;0];
else
	Zc=options.Zc;
end

epsilon = 1e-7;

% Initialisation des sorties
Zr = zeros(6,6,length(Omega)); %%% matrice d'impédance --> A + (1i/omega)*B
A = zeros(6,6,length(Omega));
B = zeros(6,6,length(Omega));
Fe = zeros(length(Omega),6); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),6); %%% effort d'excitation
Fe_FK = zeros(length(Omega),6);%%% Froude-Krilov	
	
Hq1=zeros(Ni,2); Hq3=zeros(Ni,2); Hq5=zeros(Ni,2); Hq7=zeros(Ni,2);
Dm0=zeros(Ni,Ni); Dm1=zeros(Ni,Ni);

% Initialisation des coefficients de Fourier
Aq1=zeros(Ni,2); Bq1=zeros(Nn,2); Cq1=zeros(Nj,2);
Aq3=zeros(Ni,2); Bq3=zeros(Nn,2); Cq3=zeros(Nj,2); 
Aq5=zeros(Ni,2); Bq5=zeros(Nn,2); Cq5=zeros(Nj,2);
Aq7=zeros(Ni,2); Bq7=zeros(Nn,2); Cq7=zeros(Nj,2);

L_jtau = zeros(Nj,Ni); M_ntau = zeros(Nn,Ni);

n=(1:1:Nn-1)'; lambda_n=(pi/(e1-b)).*n;
j=(1:1:Nj-1)'; lambda_j=(pi/(h-e2)).*j;


P3n=zeros(Nn,2);
% k=1
P3n(1,1) =-((e1-b) / 12)*(2 - 3*(Rbo^2 / (e1-b)^2));
P3n(2:Nn,1) =-(sqrt(2)*(e1-b).*(-1).^n)./(n.^2 .* pi^2);
% k=2
P3n(1,2)=((e1-b) / 12)*(2 - 3*( Rbo/(e1-b))^2);
P3n(2:Nn,2)=sqrt(2)*(e1-b)./(n.^2.*pi^2);

P3j=zeros(Nj,2);
P3j(1,2)=-((h-e2)/12)*(2 - 3*(Rbo/(h-e2))^2);
P3j(2:Nj,2)=-(sqrt(2)*(h-e2).*(-1).^j)./(j.^2 .* pi^2);



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


N_lambda_i = 0.5.*( 1 + sin(2.*lambda_i.*h)./(2.*lambda_i.*h));

L_jtau(1,:) = N_lambda_i.^(-.5).*sin(lambda_i.*(h-e2))./( (h-e2).*lambda_i );
for i=1:Ni
	for p=1:Nj-1
		if (lambda_j(p)/lambda_i(i) < 1+epsilon)&&(lambda_j(p)/lambda_i(i) > 1-epsilon)
			L_jtau(p+1,i) = .5*sqrt(2).*N_lambda_i(i)^(-.5);
		else
			L_jtau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*...
				( sin((lambda_j(p)-lambda_i(i))*(h-e2))/(lambda_j(p)-lambda_i(i))...
				+ sin((lambda_j(p)+lambda_i(i))*(h-e2))/(lambda_j(p)+lambda_i(i)))/(h-e2);
		end
	end
end

M_ntau(1,:) = N_lambda_i.^(-.5).*(sin(lambda_i.*(h-b))-sin(lambda_i.*(h-e1)))./((e1-b).*lambda_i);
for i=1:Ni
	for p=1:Nn-1
		if (lambda_n(p)/lambda_i(i) < 1+epsilon)&&(lambda_n(p)/lambda_i(i) > 1-epsilon)
			M_ntau(p+1,i) = sqrt(2)*N_lambda_i(i)^(-.5)*(.25*(sin(lambda_n(p)*(e1+h-2*b))-sin(lambda_n(p)*(h-e1)))/lambda_n(p)+.5*(e1-b)*cos(lambda_n(p)*(e1-h)))/(e1-b);
		else
			M_ntau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*...
				( (sin(-(lambda_n(p)-lambda_i(i))*b+lambda_n(p)*e1-lambda_i(i)*h)-sin(-(lambda_n(p)-lambda_i(i))*e1+lambda_n(p)*e1-lambda_i(i)*h))/(lambda_n(p)-lambda_i(i))...
				+ (sin(-(lambda_n(p)+lambda_i(i))*b+lambda_n(p)*e1+lambda_i(i)*h)-sin(-(lambda_n(p)+lambda_i(i))*e1+lambda_n(p)*e1+lambda_i(i)*h))/(lambda_n(p)+lambda_i(i)) ) / (e1-b);
		end
	end
end

Gamma_0j(1:Nj,1) = [0;lambda_j].*besseli(1,[0;lambda_j].*Rbo)./besseli(0,[0;lambda_j].*Rbo);
Gamma_1j(1:Nj,1) = [0;lambda_j].*besseli(2,[0;lambda_j].*Rbo)./besseli(1,[0;lambda_j].*Rbo) + 1/Rbo; Gamma_1j(1,1)=1/Rbo;

Pi_0n(1:Nn,1) = [0;lambda_n].*besseli(1,[0;lambda_n].*Rbo)./besseli(0,[0;lambda_n].*Rbo);
Pi_1n(1:Nn,1) = [0;lambda_n].*besseli(2,[0;lambda_n].*Rbo)./besseli(1,[0;lambda_n].*Rbo) + 1/Rbo; Pi_1n(1,1)=1/Rbo;

Delta_0i(1,1:Ni) = -lambda_i.*(besselk(1,lambda_i.*Rbo)./besselk(0,lambda_i.*Rbo));
Delta_1i(1,1:Ni) = -lambda_i.*(besselk(2,lambda_i.*Rbo)./besselk(1,lambda_i.*Rbo)) + 1/Rbo;

I1 = N_lambda_i.^(-.5).*( sin(lambda_i.*h)-sin(lambda_i.*(h-b)) ) ./ lambda_i;            % % Int_(-b)^(0){ Z_i^I }
I2 = N_lambda_i.^(-.5).*( sin(lambda_i.*(h-e1)) - sin(lambda_i.*(h-e2)) ) ./ lambda_i;    % % Int_(-e2)^(-e1){ Z_i^I }
I3_1 = N_lambda_i.^(-.5).*f1(lambda_i,-e1,-b,e1,h);                                       % % Int_(-e1)^(-b){ (z+e1)^2*Z_i^I }
I3_2 = N_lambda_i.^(-.5).*f1(lambda_i,-e1,-b,b,h);                                        % % Int_(-e1)^(-b){ (z+b)^2*Z_i^I }
I4 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-e2,h,h);                                          % % Int_(-h)^(-e2){ (z+h)^2*Z_i^I }
I5 = N_lambda_i.^(-.5).*f2(lambda_i,-b,0,0,h);                                            % % Int_(-b)^(0){ z*Z_i^I }
I6 = N_lambda_i.^(-.5).*f2(lambda_i,-e2,-e1,0,h);                                         % % Int_(-e2)^(-e1){ z*Z_i^I }
    
%% ***         Définition des matrices intrinsèques à la structure
Dm0(:,:) = diag(Delta_0i) - ((h-e2)/h).*(L_jtau'*(diag(Gamma_0j)*L_jtau)) - ((e1-b)/h).*(M_ntau'*(diag(Pi_0n)*M_ntau));
Dm1(:,:) = diag(Delta_1i) - ((h-e2)/h).*(L_jtau'*(diag(Gamma_1j)*L_jtau)) - ((e1-b)/h).*(M_ntau'*(diag(Pi_1n)*M_ntau));

%% Problème de RADIATION
%% SURGE MODE
% ***********************************************************************
%     Cas k = 1       Buoy oscillating, column fixed
% ***********************************************************************
Hq1(:,1) = (1/h).*I1';

Aq1(:,1) = Dm1\Hq1(:,1);
Bq1(:,1) = M_ntau*Aq1(:,1);
Cq1(:,1) = L_jtau*Aq1(:,1);

Zr(1,1,w) = -pi*rho*Rbo*(I1*Aq1(:,1));
Zr(4,1,w) = -pi*rho*Rbo*(I2*Aq1(:,1));

tmp2 = .25*Bq1(1,1)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^n.*besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq1(2:end,1);
Zr(3,1,w) = -pi*rho*(Rbo*I5*Aq1(:,1) + tmp2);

tmp1 = .25*Cq1(1,1)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Cq1(2:end,1);
tmp2 = .25*Bq1(1,1)*Rbo^3 + sqrt(2)*Rbo^2*((besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq1(2:end,1);
Zr(6,1,w) = - pi*rho*( Rbo*I6*Aq1(:,1) + tmp1 - tmp2 );
    
% ***********************************************************************
%     Cas k = 2       Buoy fixed , column oscillating 
% ***********************************************************************
Hq1(:,2) = (1/h).*I2';

Aq1(:,2) = Dm1\Hq1(:,2);
Bq1(:,2) = M_ntau*Aq1(:,2);
Cq1(:,2) = L_jtau*Aq1(:,2);

Zr(1,4,w) = -pi*rho*Rbo*(I1*Aq1(:,2));
Zr(4,4,w) = -pi*rho*Rbo*(I2*Aq1(:,2));

tmp2 = .25*Bq1(1,2)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^n.*besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq1(2:end,2);
Zr(3,4,w) = -pi*rho*(Rbo*I5*Aq1(:,2) + tmp2);

tmp1 = .25*Cq1(1,2)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Cq1(2:end,2);
tmp2 = .25*Bq1(1,2)*Rbo^3 + sqrt(2)*Rbo^2*((besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq1(2:end,2);
Zr(6,4,w) = - pi*rho*( Rbo*I6*Aq1(:,2) + tmp1 - tmp2 );


Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1,1)) / (cosh(k0*h)*besselh(1,k0*Rbo));
Fe_Haskind(w,4) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1,2)) / (cosh(k0*h)*besselh(1,k0*Rbo));

%% HEAVE MODE
% ***********************************************************************
%     Cas k = 1       Buoy oscillating, column fixed
% ***********************************************************************   
Hq3(:,1) = ((e1-b)/h).*M_ntau'*(Pi_0n.*P3n(:,1)) - (Rbo/ (2*h)).*M_ntau(1,:)';

Aq3(:,1) = Dm0\Hq3(:,1);
Bq3(:,1) = M_ntau*Aq3(:,1) + P3n(:,1);
Cq3(:,1) = L_jtau*Aq3(:,1);

tmp1 = ( ( (-1).^n.*Rbo./(lambda_n) ).*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo)) )';
Zr(2,2,w) = 2*pi*rho*( .5*Rbo^2*Bq3(1,1) + sqrt(2)*tmp1*Bq3(2:end,1) + .5*(e1-b)*(.5*Rbo^2 - (1/8)*( Rbo^4/(e1-b)^2 )) );

tmp1 = ( ((-1).^j.*Rbo./(lambda_j)).*(besseli(1,lambda_j.*Rbo)./besseli(0,lambda_j.*Rbo)) )'*Cq3(2:end,1);
tmp2 = ( (Rbo./(lambda_n)).*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo)) )'*Bq3(2:end,1);
Zr(5,2,w) = 2*pi*rho*( .5*Rbo^2*(Cq3(1,1)-Bq3(1,1)) + sqrt(2)*(tmp1-tmp2) + (1/16)*(Rbo^4/(e1-b)) );


% ***********************************************************************    
%     Cas k = 2       Buoy fixed, plate oscillating 
% ***********************************************************************
Hq3(:,2) = ((h-e2)/h).*(L_jtau'*(Gamma_0j.*P3j(:,2))) + ((e1-b)/h).*M_ntau'*(Pi_0n.*P3n(:,2)) - (.5*Rbo/h).*L_jtau(1,:)' + (.5*Rbo/h).*M_ntau(1,:)';

Aq3(:,2) = Dm0\Hq3(:,2);
Bq3(:,2) = M_ntau*Aq3(:,2) + P3n(:,2);
Cq3(:,2) = L_jtau*Aq3(:,2) + P3j(:,2); 

tmp1 = ( ((-1).^n.*Rbo./lambda_n).*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo)))';
Zr(2,5,w) =2*pi*rho*( .5*Rbo^2*Bq3(1,2) + sqrt(2)*tmp1*Bq3(2:end,2) + (1/16)*(Rbo^4/(e1-b)) );

tmp1 = ( ((-1).^j.*Rbo./lambda_j).*(besseli(1,lambda_j.*Rbo)./besseli(0,lambda_j.*Rbo)) )'*Cq3(2:end,2);
tmp2 = ( (Rbo./lambda_n).*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo)) )'*Bq3(2:end,2);
tmp3 = .5*(h-e2) * ( .5*Rbo^2 - (1/8)*( Rbo^4/(h-e2)^2) );
tmp4 = .5*(e1-b) * ( .5*Rbo^2 - (1/8)*( Rbo^4/(e1-b)^2) );
Zr(5,5,w) = 2*pi*rho*( .5*Rbo^2*(Cq3(1,2)-Bq3(1,2)) + sqrt(2)*(tmp1-tmp2) + tmp3 + tmp4);

Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1,1)) / (cosh(k0*h)*besselh(0,k0*Rbo));
Fe_Haskind(w,5) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1,2)) / (cosh(k0*h)*besselh(0,k0*Rbo));	

		
%%% ***********************************************************************

%%%                             PITCH MODE

%%% ***********************************************************************
%%% ***********************************************************************
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, plate fixed
%%% *********************************************************************** 
    P5n(1,1) = -(Rbo^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rbo^2));
    P5n(2:Nn,1) = ((sqrt(2)*Rbo)/(2*(e1-b)^2)).*f1(lambda_n,-e1,-b,e1,e1);

    Hq5(:,1) = ((e1-b)/h).*M_ntau'*(Pi_1n.*P5n(:,1)) + ((3*Rbo^2)/(8*h)).*M_ntau(1,:)' - (1/(2*h*(e1-b))).*I3_1' + (1/h).*I5';
              
    coeff = Dm1\Hq5(:,1);

    Aq5(:,1) = coeff(1:end,1);
    Bq5(:,1) = M_ntau*Aq5(:,1) + P5n(:,1);
    Cq5(:,1) = L_jtau*Aq5(:,1);

    tmp1 =  ( 1/(8*(e1-b)) ) * ( Rbo^6/6 - Rbo^4*(e1-b)^2);
    tmp2 = .25*Bq5(1,1)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^n.*besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq5(2:end,1);
    Zr(3,3,w) = - pi*rho*(Rbo*(I5*Aq5(:,1)) + tmp1 + tmp2);

    tmp1 = .25*Cq5(1,1)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Cq5(2:end,1);
    tmp2 = .25*Bq5(1,1)*Rbo^3 + sqrt(2)*Rbo^2*((besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq5(2:end,1);
    Zr(6,3,w) = - pi*rho*( Rbo*I6*Aq5(:,1) + tmp1 - (tmp2 + Rbo^6/(48*(e1-b))) );    
    
    Zr(1,3,w) = -pi*rho*(Rbo*I1*Aq5(:,1));
    Zr(4,3,w) = -pi*rho*(Rbo*I2*Aq5(:,1));

%%% ***********************************************************************    
% %     Cas k = 2       Buoy fixed, plate oscillating 
%%% ***********************************************************************

    P5n(1,2) = (Rbo^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rbo^2));
    P5n(2:Nn,2) =- ((sqrt(2)*Rbo)/(2*(e1-b)^2)).*f1(lambda_n,-e1,-b,b,e1);

    P5j(1,2) = -(Rbo^3/(24*(h-e2))) * ( 3 - 4*((h-e2)^2/Rbo^2));
    P5j(2:Nj,2) = +((sqrt(2)*Rbo)/(2*(h-e2)^2)).*f1(lambda_j,-h,-e2,h,h);  
    
    Hq5(:,2) = ((e1-b)/h).*M_ntau'*(Pi_1n.*P5n(:,2)) + ((h-e2)/h).*(L_jtau'*(Gamma_1j.*P5j(:,2))) - ((3*Rbo^2)/(8*h)).*M_ntau(1,:)' + (1/(2*h*(e1-b))).*I3_2'...
                 + ((3*Rbo^2)/(8*h)).*L_jtau(1,:)' - (1/(2*h*(h-e2))).*I4' + (1/h).*I6';

    Aq5(:,2) = Dm1\Hq5(:,2);
    Bq5(:,2) = M_ntau*Aq5(:,2) + P5n(:,2);
    Cq5(:,2) = L_jtau*Aq5(:,2) + P5j(:,2);

    tmp1 =  -( 1/(8*(e1-b)) ) * ( Rbo^6/6 );
    tmp2 = .25*Bq5(1,2)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^n.*besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq5(2:end,2);
    Zr(3,6,w) = - pi*rho*(Rbo*(I5*Aq5(:,2)) + tmp1 + tmp2);

    tmp1 = .25*Cq5(1,2)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Cq5(2:end,2);
    tmp2 = ( 1/(8*(h-e2)) ) * ( Rbo^6/6 - Rbo^4*(h-e2)^2);
    tmp3 = .25*Bq5(1,2)*Rbo^3 + sqrt(2)*Rbo^2*((besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq5(2:end,2);
    tmp4 = (- 1/(8*(e1-b)) ) * ( Rbo^6/6 - Rbo^4*(e1-b)^2);
    Zr(6,6,w) = - pi*rho*( Rbo*I6*Aq5(:,2) + tmp1 + tmp2 - (tmp3 + tmp4) );
    
    Zr(1,6,w) = -pi*rho*(Rbo*I1*Aq5(:,2));
    Zr(4,6,w) = -pi*rho*(Rbo*I2*Aq5(:,2));

Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1,1)) / (cosh(k0*h)*besselh(1,k0*Rbo));
Fe_Haskind(w,6) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1,2)) / (cosh(k0*h)*besselh(1,k0*Rbo));

%%% ***********************************************************************
%%% ***********************************************************************

%%%             Problème de DIFFRACTION

%%% ***********************************************************************
%%% ********************************************************************
a0 = 1i*omega;
%% m=0
B0=(-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));

Hq7(:,1) =  B0.*(besselj(0,k0*Rbo).*(((e1-b)/h).*(M_ntau'*(Pi_0n.*M_ntau(:,1))) + ((h-e2)/h).*(L_jtau'*(Gamma_0j.*L_jtau(:,1)))) + k0*besselj(1,k0*Rbo).*[1;zeros(Ni-1,1)]);

Aq7(:,1) = Dm0\Hq7(:,1);
Bq7(:,1) = M_ntau*Aq7(:,1)+B0.*besselj(0,k0*Rbo).*M_ntau(:,1);
Cq7(:,1) = L_jtau*Aq7(:,1)+B0.*besselj(0,k0*Rbo).*L_jtau(:,1);

tmp1 = ( ((-1).^n.*Rbo./lambda_n).*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo)))';
Fe(w,2) = 2*pi*a0*rho*( .5*Rbo^2*Bq7(1,1) + sqrt(2)*tmp1*Bq7(2:end,1));

tmp1 = (((-1).^j.*Rbo./lambda_j).*(besseli(1,lambda_j.*Rbo)./besseli(0,lambda_j.*Rbo)))' * Cq7(2:end,1);
tmp2 = ( (Rbo./lambda_n).* (besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo)) )' * Bq7(2:end,1);
Fe(w,5) = 2*pi*a0*rho*( .5*Rbo^2*(Cq7(1,1)-Bq7(1,1)) + sqrt(2)*(tmp1 - tmp2) );
		
%% m=1
B1 = 2*1i*B0;

Hq7(:,2) =  B1.*(besselj(1,k0*Rbo).*(((e1-b)/h).*(M_ntau'*(Pi_1n.*M_ntau(:,1))) + ((h-e2)/h).*(L_jtau'*(Gamma_1j.*L_jtau(:,1)))) - (k0*besselj(0,k0*Rbo)-besselj(1,k0*Rbo)/Rbo).*[1;zeros(Ni-1,1)]);

Aq7(:,2) = Dm1\Hq7(:,2);
Bq7(:,2) = M_ntau*Aq7(:,2)+B1.*besselj(1,k0*Rbo).*M_ntau(:,1);
Cq7(:,2) = L_jtau*Aq7(:,2)+B1.*besselj(1,k0*Rbo).*L_jtau(:,1);

Fe(w,1) = -pi*a0*rho*Rbo*(B1*besselj(1,k0*Rbo)*I1(1) + I1*Aq7(:,2));
Fe(w,4) = -pi*a0*rho*Rbo*(B1*besselj(1,k0*Rbo)*I2(1) + I2*Aq7(:,2));

tmp1 = .25*Bq7(1,2)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^n.*besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq7(2:end,2);
Fe(w,3) = - pi*a0*rho*( Rbo*(B1*besselj(1,k0*Rbo)*I5(1) + I5*Aq7(:,2)) + tmp1);

tmp1 = .25*Cq7(1,2)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Cq7(2:end,2);
tmp2 = .25*Bq7(1,2)*Rbo^3 + sqrt(2)*Rbo^2*((besseli(2,lambda_n.*Rbo))./(lambda_n.*besseli(1,lambda_n.*Rbo)))'*Bq7(2:end,2);
Fe(w,6) = -pi*a0*rho*( Rbo*(B1*besselj(1,k0*Rbo)*I6(1) + I6*Aq7(:,2)) + tmp1 - tmp2);

end% END OF FOR LOOP

%%	EXIT
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
