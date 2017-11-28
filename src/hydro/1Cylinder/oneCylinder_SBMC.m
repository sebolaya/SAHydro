function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK]=oneCylinder_SBMC(Omega, depth, WECstructure, options)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%				SUBMERGED BOTTOM MOUNTED CYLINDER
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the hydrodynamic coefficients for the submerged bottom
% mounted cylinder.
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
% - Fe : Matrix length(Omega)x3x1 of exciation forces (complex values)
% - A  : Matrix (3x1)x(3x1)xlength(Omega) of added mass coefficients
% - B  : Matrix (3x1)x(3x1)xlength(Omega) of radiation damping coefficients
% - A_inf : Matrix (3x1)x(3x1)xlength(Omega) of infinite added mass
% - Fe_Haskind : Fe computation with the Haskind's theorem 
%				 (for verification purposes)
% - Fe_FK : Froude-Krilov efforts -> Fe without diffraction phenomena
%
% *************************************************************************
%                          WEC Description
%
%
%           |------------------->| (Rbo)
%           .                                                              
% --------  +----------------------------------------- (z=0)
%           |                   ^       ^
%           .                   |       |
%           |                   |       |
%           .                   v (e1)  .
%           |--------------------+      .
%           .cccccccccccccccccccc|      .
%           |cccccccccccccccccccc|      |
%           .cccccccccccccccccccc|      |
%   ________|cccccccccccccccccccc+______v (h)_______
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
h = depth;

Rbo = WECstructure.Rbo;
e1 = WECstructure.e1;

if ~isfield(options,'Trunc')
	Ni=20; Nn=60;
else
	Ni = options.Trunc.Ni;
	Nn = options.Trunc.Nn;
end

if ~isfield(options,'Zc')
	Zc=0;
else
	Zc=options.Zc;
end

epsilon = 1e-7;

%%% Initialisation des sorties
Zr = zeros(3,3,length(Omega)+1); %%% matrice d'impédance --> A + (1i/omega)*B
A = zeros(3,3,length(Omega));
B = zeros(3,3,length(Omega));
Fe = zeros(length(Omega)+1,3); %%% force d'excitation
Fe_Haskind = zeros(length(Omega)+1,3); %%% effort d'excitation calculée à partir de la formulation de Haskind
Fe_FK = zeros(length(Omega)+1,3);%%% Froude-Krilov

Hq1=zeros(Ni,1); Hq3=zeros(Ni,1); Hq5=zeros(Ni,1); Hq7 = zeros(Ni,2);
Dm0=zeros(Ni,Ni); Dm1=zeros(Ni,Ni);


%% Calcul des solutions particulières
P3n=zeros(Nn,1);
P5n=zeros(Nn,1);

for w=1:length(Omega)+1
	
clc;
disp([num2str(round(w*100/(length(Omega)+1))),'%']);

if w<=length(Omega)
	omega = Omega(w);
	lambda_i=fcn_waveDispersion(Omega(w),h,Ni);
	lambda_n=fcn_waveDispersion(Omega(w),e1,Nn); lambda_n=conj(lambda_n)';
	k0 = -imag(lambda_i(1));
else
	omega = Inf;
	i=1:1:Ni; n=(1:1:Nn)';
	lambda_i=.5*(2*i-1)*pi/h;
	lambda_n=.5*(2*n-1)*pi/e1;
end

N_lambda_i=0.5.*(1+sin(2.*lambda_i.*h)./(2.*lambda_i.*h));
N_lambda_n=0.5.*(1+sin(2.*lambda_n.*e1)./(2.*lambda_n.*e1));

M_ntau = zeros(Nn,Ni);
for i=1:Ni
	for p=1:Nn
		if (lambda_n(p)/lambda_i(i) < 1+epsilon)&&(lambda_n(p)/lambda_i(i) > 1-epsilon)
			M_ntau(p,i) = (1/e1)*(N_lambda_i(i)^(-.5)*N_lambda_n(p).^(-.5))*...
				(.25*(sin(lambda_i(i)*(h+e1)) - sin(lambda_i(i)*(h-e1)))/lambda_i(i) + e1*.5*cos(lambda_i(i)*(h-e1)));
		else
			M_ntau(p,i) = (1/e1)*(N_lambda_i(i)^(-.5)*N_lambda_n(p).^(-.5))*...
				.5*( (sin(lambda_n(p)*e1-lambda_i(i)*h)-sin(lambda_n(p)*e1-lambda_i(i)*h - e1*(lambda_n(p)-lambda_i(i))))/(lambda_n(p)-lambda_i(i))...
						+ (sin(lambda_n(p)*e1+lambda_i(i)*h)-sin(lambda_n(p)*e1+lambda_i(i)*h - e1*(lambda_n(p)+lambda_i(i))))/(lambda_n(p)+lambda_i(i)));
		end
	end
end

Delta_0i=-lambda_i.*(besselk(1,lambda_i.*Rbo)./besselk(0,lambda_i.*Rbo));
Delta_1i=-lambda_i.*(besselk(2,lambda_i.*Rbo)./besselk(1,lambda_i.*Rbo)) + 1/Rbo;

Pi_0n=lambda_n.*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo));
Pi_1n=lambda_n.*(besseli(2,lambda_n.*Rbo)./besseli(1,lambda_n.*Rbo)) + 1/Rbo;

I1=N_lambda_i.^(-.5).*sin(lambda_i.*(h-e1))./lambda_i;	% Int_(-h)^(-e1){ Zi(z) dz}
I2=N_lambda_n.^(-.5).*f2(lambda_n,-e1,0,0,e1);			% Int_(-e1)^(0){ z*Zn(z) dz}
I3=N_lambda_i.^(-.5).*f2(lambda_i,-h,-e1,-Zc,h);		% Int_(-h)^(-e1){ (z-Zc)*Zi(z) dz}

%% Définition des matrices intrinsèques à la structure
Dm0=diag(Delta_0i)-(e1/h).*(M_ntau'*(diag(Pi_0n)*M_ntau));
Dm1=diag(Delta_1i)-(e1/h).*(M_ntau'*(diag(Pi_1n)*M_ntau));

%% Problème de RADIATION
% SURGE MODE
Hq1(:,1) = I1'/h;

Aq1=Dm1\Hq1(:,1);
Bq1=M_ntau*Aq1;
		
Zr(1,1,w)=-pi*rho*Rbo*(I1*Aq1);

tmp1=(Rbo^2*N_lambda_n.^(-.5).*(besseli(2,lambda_n.*Rbo)./besseli(1,lambda_n.*Rbo))./lambda_n)'*Bq1;
Zr(3,1,w) = -pi*rho*(Rbo*I3*Aq1 - tmp1);

%%% HASKIND's relation
Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1))/(cosh(k0*h)*besselh(1,k0*Rbo));

% HEAVE MODE
P3n(:,1)=-(I2 + (g/omega^2).*N_lambda_n.^(-.5).*sin(lambda_n.*e1)./lambda_n)./e1;    
Hq3(:,1)=(e1/h).*(M_ntau'*(Pi_0n.*P3n));

Aq3=Dm0\Hq3(:,1);
Bq3=M_ntau*Aq3+P3n(:,1);

tmp1=.5*(g/omega^2 - e1)*Rbo^2;
tmp2=(Rbo*N_lambda_n.^(-.5).*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo))./lambda_n)'*Bq3;
Zr(2,2,w) = -2*pi*rho*(tmp1+tmp2);

%%% HASKIND's relation
Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1)) / (cosh(k0*h)*besselh(0,k0*Rbo));

% PITCH MODE
P5n(:,1)=(Rbo/e1)*N_lambda_n.^(-.5).*f2(lambda_n,-e1,0,g/omega^2,e1);
tmp1 = N_lambda_i.^(-.5).*f2(lambda_i,-e1,0,g/omega^2,h);

Hq5(:,1)=(e1/h).*(M_ntau'*(Pi_1n.*P5n(:,1))) + (1/h).*I3' - (1/h).*tmp1';

Aq5=Dm1\Hq5(:,1);
Bq5=M_ntau*Aq5+P5n(:,1);

tmp1=-.25*(g/omega^2 - e1)*Rbo^4;
tmp2=(Rbo^2*N_lambda_n.^(-.5).*(besseli(2,lambda_n.*Rbo)./besseli(1,lambda_n.*Rbo))./lambda_n)'*Bq5;
Zr(3,3,w)=-pi*rho*(Rbo*I3*Aq5-(tmp1+tmp2));
Zr(1,3,w)=-pi*rho*Rbo*(I1*Aq5);

% *** HASKIND's relation
Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1)) / (cosh(k0*h)*besselh(1,k0*Rbo));


%% Problème de DIFFRACTION
a0 = 1i*omega;

% m=0
B0 = (-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));
Hq7(:,1)=B0.*(besselj(0,k0*Rbo)*(e1/h).*(M_ntau'*(Pi_0n.*M_ntau(:,1))) + k0*besselj(1,k0*Rbo).*[1;zeros(Ni-1,1)]);

Aq7(:,1) = Dm0\Hq7(:,1);
Bq7(:,1) = M_ntau*Aq7(:,1) + B0*besselj(0,k0*Rbo).*M_ntau(:,1);

tmp1=(Rbo*N_lambda_n.^(-.5).*(besseli(1,lambda_n.*Rbo)./besseli(0,lambda_n.*Rbo))./lambda_n)'*Bq7(:,1);
Fe(w,2)=-2*pi*rho*a0*tmp1;

%% m=1
B1 = 2*1i*B0;
Hq7(:,2) = B1*(besselj(1,k0*Rbo)*(e1/h).*(M_ntau'*(Pi_1n.*M_ntau(:,1))) - (k0*besselj(0,k0*Rbo)-besselj(1,k0*Rbo)/Rbo).*[1;zeros(Ni-1,1)]);

Aq7_m1 = Dm1\Hq7(:,2);
Bq7_m1 = M_ntau*Aq7_m1 + B1*besselj(1,k0*Rbo).*M_ntau(:,1);

Fe(w,1) = -pi*a0*rho*Rbo*(B1*besselj(1,k0*Rbo)*I1(1) + I1*Aq7_m1);

tmp1=(Rbo^2*N_lambda_n.^(-.5).*(besseli(2,lambda_n.*Rbo)./besseli(1,lambda_n.*Rbo))./lambda_n)'*Bq7_m1;
Fe(w,3)=-pi*a0*rho*(Rbo*(B1*besselj(1,k0*Rbo)*I3(1) + I3*Aq7_m1) - tmp1);

end% END OF FOR() LOOP

%% EXIT

for i=1:3
	for j=1:3
		A(i,j,:) = real(Zr(i,j,1:end-1));
		B(i,j,:) = squeeze(imag(Zr(i,j,1:end-1))).*Omega;
	end
end
A_inf = real(Zr(:,:,end));
Fe = conj(Fe(1:end-1,:));
Fe_Haskind = conj(Fe_Haskind(1:end-1,:));
Fe_FK = conj(Fe_FK(1:end-1,:));

save([DIR,filesep,'DATA/hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');

end% END OF FUNCTION

%% Additional functions
function out = f1(alpha, a, b, c, d )
%%% Int_(a)^(b){(z+c)^2*cos(alpha*(z+d))}
epsilon = 1e-8;
[n,m]=size(alpha);
out=zeros(n,m);
for i=1:length(alpha)
	if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
        out(i) = ((b+c)^3 - (a+c)^3) / 3;
	else
		out(i) = (((b+c)^2 - 2./alpha(i).^2).*sin(alpha(i).*(b+d)) - ((a+c)^2 - 2./alpha(i).^2).*sin((i).*(a+d)))./alpha(i)...
					+ 2.*((b+c).*cos(alpha(i).*(b+d)) - (a+c).*cos(alpha(i).*(a+d)))./alpha(i).^2;
	end
end
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


