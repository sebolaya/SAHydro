function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = oneCylinder(Omega, depth, WECstructure, options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the hydrodynamic coefficients for the one cylinder 
% structure depicted below.
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
%                         WEC Description
%
%           |-----------------------> (Rbo)
%           .                                                              
% --------  +-----------------------+ --------------------------- (z=0)
%           .ccccccccccccccccccccccc|  ^                    ^
%           |ccccccccccccccccccccccc|  |                    |
%           .ccccccccccccccccccccccc|  |                    |
%           |-----------------------+  v (b)                .
%           .                                               .
%           |                                               .
%                                                           |
%                                                           v (h)
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
b = WECstructure.b;

if ~isfield(options,'Trunc')
	Ni=80; Nj=100;
else
	Ni=options.Trunc.Ni;
	Nj=options.Trunc.Nj;
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
Fe_Haskind = zeros(length(Omega),3); %%% effort d'excitation calculée à partir de la formulation de Haskind
Fe_FK = zeros(length(Omega),3);%%% Froude-Krilov

j = (1:1:Nj-1)';
lambda_j = (pi/(h-b)).*j;

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
		lambda_i = .5*(2*i-1)*pi/h;
	end

	N_lambda_i = 0.5.*( 1 + sin(2.*lambda_i.*h)./(2.*lambda_i.*h));

	L_jtau = zeros(Nj,Ni);
	L_jtau(1,:) = N_lambda_i.^(-.5).*sin(lambda_i.*(h-b))./( (h-b).*lambda_i );
	for i=1:Ni
		for p = 1 : Nj-1
		    if (lambda_j(p)/lambda_i(i) < 1+epsilon)&&(lambda_j(p)/lambda_i(i) > 1-epsilon)
		        L_jtau(p+1,i) = .5*sqrt(2).*N_lambda_i(i)^(-.5);
		    else
		        L_jtau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*(sin((lambda_j(p)-lambda_i(i))*(h-b))/(lambda_j(p)-lambda_i(i))...
					+ sin((lambda_j(p)+lambda_i(i))*(h-b))/(lambda_j(p)+lambda_i(i)))/(h-b);
		    end
		end
	end
	
	Gamma_0j(1,1) = 0;
	Gamma_0j(2:Nj,1) = lambda_j.*Rbo.*besseli(1,lambda_j.*Rbo)./besseli(0,lambda_j.*Rbo);

	Gamma_1j(1,1) = 1;
	Gamma_1j(2:Nj,1) = lambda_j.*Rbo.*besseli(2,lambda_j.*Rbo)./besseli(1,lambda_j.*Rbo) + 1;

	Delta_0i = -lambda_i.*Rbo.*(besselk(1,lambda_i.*Rbo) ./ besselk(0,lambda_i.*Rbo));
	Delta_1i = -lambda_i.*Rbo.*(besselk(2,lambda_i.*Rbo)./besselk(1,lambda_i.*Rbo)) + 1;
	

	I1 = N_lambda_i.^(-.5).*( sin(lambda_i.*h)-sin(lambda_i.*(h-b)) ) ./ lambda_i;		%%% Int_(-b)^(0){ Z_i^I }
	I5 = N_lambda_i.^(-.5).*f2(lambda_i,-b,0,-Zc,h);									%%% Int_(-b)^(0){ (z-Zc)*Z_i^I }
	
	% *** Définition des matrices intrinsèques à la structure
	D0_tau_i = (diag(Delta_0i) - ((h-b)/h).*(L_jtau'*(diag(Gamma_0j)*L_jtau)));%*Nik;
	D1_tau_i = diag(Delta_1i) - ((h-b)/h).*(L_jtau'*(diag(Gamma_1j)*L_jtau));
	
% 	R1_tau_i=real(D1_tau_i); Im=1i*imag(D1_tau_i(1,1));

%%                     RADIATION problem
%%                         SURGE MODE
H_1tau_r1 = I1'.*(Rbo/h);

Aq1 = D1_tau_i\H_1tau_r1;
Bq1 = L_jtau*Aq1;


% Phi_1q=R1_tau_i\H_1tau_r1;
% alpha=Im/H_1tau_r1(1);
% Aq1_2=Phi_1q; Aq1_2(1)=Aq1_2(1)/(1+alpha*Phi_1q(1));

Zr(1,1,w)=-pi*rho*Rbo*(I1*Aq1);

tmp1 = .25*Bq1(1)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Bq1(2:end);
Zr(3,1,w) = -pi*rho*(Rbo*I5*Aq1 + tmp1);

%%% HASKIND's relation
Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1))/(cosh(k0*h)*besselh(1,k0*Rbo));

%% HEAVE MODE
%%% solution particulière pour r = R
p_03 = ( (h-b) / 12 ) * ( 2 - 3*(Rbo^2 / (h-b)^2 ) );
p_j3 = (sqrt(2)*(h-b).*(-1).^j) ./ (j.*pi).^2;
P_j3 = -[p_03;p_j3];

H_0tau_r3 = ((h-b)/h).*(L_jtau'*(Gamma_0j.*P_j3)) - (.5*Rbo^2/h).*L_jtau(1,:)';

% Ak = pinv(D0_tau_i)*H_0tau_r3;
% Aq3=Nik*Ak;

Aq3 = D0_tau_i\H_0tau_r3;
Bq3 = L_jtau*Aq3 + P_j3;

tmp1 = ( sqrt(2)*Rbo.*((-1).^j./lambda_j).*(besseli(1,lambda_j.*Rbo)./besseli(0,lambda_j.*Rbo)) )'*Bq3(2:end);
tmp2 = .25*(h-b)*( Rbo^2 - .25*Rbo^4/(h-b)^2 );                         % solution particulière

Zr(2,2,w) = 2*pi*rho*( tmp2 + .5*Bq3(1)*Rbo^2 + tmp1);

%%% HASKIND's relation
Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3(1)) / (cosh(k0*h)*besselh(0,k0*Rbo));

%%                         PITCH MODE
q_05 = ( Rbo^3 / (24*(h-b)) ) * ( 3 - 4*((h-b)^2 / Rbo^2) );
q_j5 = -(sqrt(2)*Rbo*(h-b).*(-1).^j) ./ (j.*pi).^2;
Q_j5 = -[q_05;q_j5];

I_i5_1 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-b,h,h);

H_1tau_r5 =  ((h-b)/h).*(L_jtau'*(Gamma_1j.*Q_j5)) + ((3*Rbo^3)/(8*h)).*L_jtau(1,:)' - (Rbo/(2*h*(h-b))).*I_i5_1' + (Rbo/h).*I5';

Aq5 = D1_tau_i\H_1tau_r5;
Bq5 = L_jtau*Aq5 + Q_j5;

tmp1 =  ( 1/(8*(h-b)) ) * ( Rbo^6/6 - Rbo^4*(h-b)^2);
tmp2 = .25*Bq5(1)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Bq5(2:end);

Zr(3,3,w) = - pi*rho*(Rbo*(I5*Aq5) + tmp1 + tmp2);
Zr(1,3,w) = - pi*rho*Rbo*(I1*Aq5);

% *** HASKIND's relation
Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1)) / (cosh(k0*h)*besselh(1,k0*Rbo));

%% DIFFRACTION problem
a0 = 1i*omega;
%% m=0
B0 = (-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));

H_0tau = B0.*( ((h-b)/h)*besselj(0,k0*Rbo).*(L_jtau'*(Gamma_0j.*L_jtau(:,1))) + k0*Rbo*besselj(1,k0*Rbo).*[1;zeros(Ni-1,1)] );

% Ak = pinv(D0_tau_i)*H_0tau;
% Aq7_m0=Nik*Ak;

Aq7_m0 = D0_tau_i\H_0tau;
Bq7_m0 = L_jtau*Aq7_m0 + B0*besselj(0,k0*Rbo).*L_jtau(:,1);

tmp1=(sqrt(2)*Rbo.*((-1).^j./(lambda_j)).*(besseli(1,lambda_j.*Rbo)./besseli(0,lambda_j.*Rbo)) )' * Bq7_m0(2:end);
Fe(w,2) = 2*pi*rho*a0*( .5*Bq7_m0(1)*Rbo^2 + tmp1);

%% m=1
B1 = 2*1i*B0;
H_1tau = B1*( ((h-b)/h)*besselj(1,k0*Rbo).*(L_jtau'*(Gamma_1j.*L_jtau(:,1))) - (k0*Rbo*besselj(0,k0*Rbo)-besselj(1,k0*Rbo)).*[1;zeros(Ni-1,1)] );

Aq7_m1 = D1_tau_i\H_1tau;
Bq7_m1 = L_jtau*Aq7_m1 + B1*besselj(1,k0*Rbo).*L_jtau(:,1);

Fe(w,1) = -pi*a0*rho*Rbo*( B1*besselj(1,k0*Rbo)*I1(1) + I1*Aq7_m1);

tmp1 = .25*Bq7_m1(1)*Rbo^3 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.*besseli(1,lambda_j.*Rbo)))'*Bq7_m1(2:end);
Fe(w,3) = -pi*a0*rho*( Rbo*(B1*besselj(1,k0*Rbo)*I5(1) + I5*Aq7_m1) + tmp1);



if isfield(options,'matchingConditions')
	if (options.matchingConditions==1)&&(omega==options.matchingFrequency)
		dz=.1; dr=.1;
		z = -10:dz:0;
		r = 0:dr:10;

		u_r = zeros(length(z),length(r));
		u_z = zeros(length(z),length(r));
		for i=1:length(z)
			for j=1:length(r)
				if ((z(i)>-b)&&(r(j)<Rbo))
					u_r(i,j) = NaN;% pour ne pas afficher de couleurs
					u_z(i,j) = NaN;
				elseif ((z(i) <= -b)&&(r(j) < Rbo))
					u_r(i,j) = abs( -.5*r(j)/(h-b) + sqrt(2)*( lambda_j'.*(besseli(1,lambda_j'*r(j))./besseli(0,lambda_j'*Rbo)).*cos(lambda_j'.*(z(i)+h)) )*Bq3(2:end));
					u_z(i,j) = abs( (z(i)+h)/(h-b) - sqrt(2)*( lambda_j'.*(besseli(0,lambda_j'*r(j))./besseli(0,lambda_j'*Rbo)).*sin(lambda_j'.*(z(i)+h)) )*Bq3(2:end));
				elseif (r(j) >= Rbo)
					u_r(i,j) = abs(-(lambda_i.*besselk(1,lambda_i*r(j))./besselk(0,lambda_i*Rbo).*sqrt(1./N_lambda_i).*cos(lambda_i*(z(i)+h)) )*Aq3);
					u_z(i,j) = abs((lambda_i.*besselk(0,lambda_i*r(j))./besselk(0,lambda_i*Rbo).*sqrt(1./N_lambda_i).*sin(lambda_i*(z(i)+h)) )*Aq3);
				end
			end
		end

		figure,  grid on, hold on;
% 			subplot(1,2,1), grid on, hold on;
		contourf(r,z,u_r,100);
		line([0 Rbo],[-b -b],'Color','k','LineWidth',2); line([Rbo Rbo],[0 -b],'Color','k','LineWidth',2);
		shading flat;
		colormap('Jet');
		colorbar;
		figure,  grid on, hold on;
% 			subplot(1,2,2), grid on, hold on;
		contourf(r,z,u_z,100);
		line([0 Rbo],[-b -b],'Color','k','LineWidth',2); line([Rbo Rbo],[0 -b],'Color','k','LineWidth',2);
		shading flat;
		colormap('Jet');
		colorbar;
	end
end


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

save([DIR,filesep,'DATA',filesep,'hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');

M = rho*pi*Rbo^2*b;
Kh = pi*rho*g*Rbo^2;

save([DIR,filesep,'DATA',filesep,'WECstructure.mat'],'Rbo','b','M','Kh');

if isfield(options,'IRF')
	if options.IRF==1	
		[irf,t]=IRFfcn(Omega,B,25,.1);
		save([DIR,filesep,'DATA',filesep,'IRF.mat'],'irf','t');
	end
end

end% END OF FUNCTION

% function out = f1( alpha, a, b, c, d )
% %%% Int_(a)^(b){(z+c)^2*cos(alpha*(z+d))}
%     if alpha == 0
%         out = ((b+c)^3 - (a+c)^3) / 3;
% 	else
% 		out = (((b+c)^2 - 2./alpha.^2).*sin(alpha.*(b+d)) - ((a+c)^2 - 2./alpha.^2).*sin(alpha.*(a+d)))./alpha...
% 					+ 2.*((b+c).*cos(alpha.*(b+d)) - (a+c).*cos(alpha.*(a+d)))./alpha.^2;
%     end
% end
%% Additional functions
function out = f1(alpha, a, b, c, d )
%%% Int_(a)^(b){(z+c)^2*cos(alpha*(z+d))}
[n,m]=size(alpha);
out=zeros(n,m);
out(1:end)=((b+c)^2-2./alpha.^2).*sin(alpha.*(b+d))./alpha + 2*(b+c).*cos(alpha.*(b+d))./alpha.^2-...
	      (((a+c)^2-2./alpha.^2).*sin(alpha.*(a+d))./alpha + 2*(a+c).*cos(alpha.*(a+d))./alpha.^2);
end

function [ out ] = f2( alpha, a, b, c, d )
%%% Int_(a)^(b){(z+c)*cos(alpha*(z+d))}
    if alpha == 0
		out = .5*((b+c)^2 - (a+c)^2);
	else
		out = ((b+c).*sin(alpha.*(b+d)) - (a+c).*sin(alpha.*(a+d)))./alpha  + (cos(alpha.*(b+d)) - cos(alpha.*(a+d)))./alpha.^2;
    end
end



