function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = oneCylinder_GALERKIN(Omega, depth, WECstructure, options)
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
% Date : 2016/06/05
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
	Ni=100; Nj=100; Nn=5;
else
	Ni=options.Trunc.Ni;
	Nj=options.Trunc.Nj;
	Nn=options.Trunc.Nn;
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
Fe = zeros(length(Omega),3); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),3); %%% effort d'excitation calculée à partir de la formulation de Haskind
Fe_FK = zeros(length(Omega),3);%%% Froude-Krilov

j=(1:1:Nj-1)';
lambda_j=(pi/(h-b)).*j;

for w=1:length(Omega)+1
	clc;
	disp([num2str(round(w*100/(length(Omega)+1))),'%']);
	
	if w<=length(Omega)
		omega=Omega(w);
		lambda_i=fcn_waveDispersion(Omega(w),h,Ni);
		k0=-imag(lambda_i(1));
	else
		omega = Inf;% omega >> 1
		i=1:1:Ni;
		lambda_i=.5*(2*i-1)*pi/h;
	end

	nu=1/6;
	
	N_lambda_i=.5.*(1+sin(2.*lambda_i.*h)./(2.*lambda_i.*h));

	Lni=zeros(Nn,Ni);
	Lnj=zeros(Nn,Nj); Lnj(1,1)=(1)/((2^nu)*nu*gamma(nu));
	for n=0:Nn-1
		if w<=length(Omega)
			Lni(n+1,1)=N_lambda_i(1).^(-.5).*besseli(nu+2*n,k0.*(h-b))./(k0.*(h-b)).^nu;%
			Lni(n+1,2:end)=(-1)^n.*N_lambda_i(2:end).^(-.5).*besselj(nu+2*n,lambda_i(2:end).*(h-b))./(lambda_i(2:end).*(h-b)).^nu;
		else
			Lni(n+1,:)=(-1)^n.*N_lambda_i.^(-.5).*besselj(nu+2*n,lambda_i.*(h-b))./(lambda_i.*(h-b)).^nu;
		end
		Lnj(n+1,2:end)=(-1)^n.*sqrt(2).*besselj(nu+2*n,lambda_j.*(h-b))./(lambda_j.*(h-b)).^nu;
	end
	
	gmj=zeros(2,Nj); Gmi=zeros(2,Ni);
	
	gmj(1,1)=0; gmj(2,1)=Rbo;
	gmj(1,2:Nj)=besseli(0,lambda_j.*Rbo)./(lambda_j.*(besseli(1,lambda_j.*Rbo)));
	gmj(2,2:Nj)=besseli(1,lambda_j.*Rbo)./(lambda_j.*(besseli(2,lambda_j.*Rbo) + besseli(1,lambda_j.*Rbo)./(lambda_j.*Rbo)));

	Gmi(1,:)=besselk(0,lambda_i.*Rbo)./(lambda_i.*(-besselk(1,lambda_i.*Rbo)));
	Gmi(2,:)=besselk(1,lambda_i.*Rbo)./(lambda_i.*(-besselk(2,lambda_i.*Rbo) + besselk(1,lambda_i.*Rbo)./(lambda_i.*Rbo)));
	
	Dm=zeros(Nn,Nn,2);
	
	Dm(:,:,1)=-(1/h)*Lni*conj((ones(Nn,1)*Gmi(1,:)).*Lni)'+(1/(h-b))*Lnj(:,2:end)*((ones(Nn,1)*gmj(1,2:end)).*Lnj(:,2:end))';
	Dm(:,:,2)=-(1/h)*Lni*conj((ones(Nn,1)*Gmi(2,:)).*Lni)'+(1/(h-b))*Lnj*((ones(Nn,1)*gmj(2,:)).*Lnj)';
	
% 	for r=1:Nn
% 		for n=1:Nn
% 			tmp(r,n)=sum(gmj(2,:).*Lnj(r,:).*Lnj(n,:))/(h-b)-sum(Gmi(2,:).*Lni(r,:).*Lni(n,:))/(h);
% 		end
% 	end
	
I1=N_lambda_i.^(-.5).*(sin(lambda_i.*h)-sin(lambda_i.*(h-b)))./lambda_i;   % Int_(-b)^(0){ Z_i^I }
I5=N_lambda_i.^(-.5).*f2(lambda_i,-b,0,-Zc,h);                             % Int_(-b)^(0){ (z-Zc)*Z_i^I }	

ami=zeros(Ni,3); bmj=zeros(Nj,3);
Ami=zeros(Ni,2); Bmj=zeros(Nj,2);
Hq=zeros(Nn,3);
%% RADIATION problem
% SURGE MODE
Hq(:,1)=Lni*conj(((I1./h).*Gmi(2,:))');

An=Dm(:,:,2)\Hq(:,1);

Umi=(1/h)*(Lni'*An) + (I1./h)';
% Umi=(1/h)*(conj(Lni')*An) + (I1./h)';
Umj=(1/(h-b))*(Lnj'*An);

ami(:,1)=Umi./conj(lambda_i.*(-besselk(2,lambda_i.*Rbo) + besselk(1,lambda_i.*Rbo)./(lambda_i.*Rbo)))';

Zr(1,1,w)=-pi*rho*Rbo*((I1.*Gmi(2,:))*Umi);

Fe_Haskind(w,1)=(4*rho*g*h*sqrt(N_lambda_i(1))*ami(1,1))/(cosh(k0*h));
% Fe_Haskind(w,1)=(4*rho*g*h*sqrt(N_lambda_i(1))*Umi(1))/(cosh(k0*h)*k0*(-besselh(2,k0*Rbo)+besselh(1,k0*Rbo)/(k0*Rbo)));

% HEAVE MODE
Hq(1,2)=-.25*(h-b)*(1/(nu+1)-(Rbo/(h-b))^2)/((2^nu)*nu*gamma(nu));
Hq(2,2)=-.25*(h-b)/(nu*(nu+1)*(2^nu)*(2+nu)*gamma(nu));

An=[0 conj(Lnj(:,1)');Lnj(:,1) Dm(:,:,1)]\[-.5*Rbo;Hq(:,2)];

Umi=(1/h)*(Lni'*An(2:end));
% Umi=(1/h)*(conj(Lni')*An(2:end)); 
Umj=(1/(h-b))*(Lnj(:,2:end)'*An(2:end));

ami(:,2)=Umi./conj(lambda_i.*(-besselk(1,lambda_i.*Rbo)))';
bmj(1,2)=An(1);
bmj(2:end,2)=Umj./(lambda_j.*besseli(1,lambda_j.*Rbo));

tmp1=.25*(h-b)*(Rbo^2-.25*Rbo^4/(h-b)^2);
% tmp2=sqrt(2)*Rbo.*(((-1).^j.*besseli(1,lambda_j.*Rbo))./(lambda_j))'*bmj(2:end,2);
tmp2=sqrt(2)*Rbo.*(((-1).^j.*besseli(1,lambda_j.*Rbo))./(lambda_j.^2.*besseli(1,lambda_j.*Rbo)))'*Umj;
Zr(2,2,w)=2*pi*rho*(tmp1 + .5*bmj(1,2)*Rbo^2 + tmp2);

% Fe_Haskind(w,2)=-(4*1i*rho*g*h*sqrt(N_lambda_i(1))*ami(1,2))/cosh(k0*h);
Fe_Haskind(w,2)=-(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Umi(1))/(cosh(k0*h)*k0*(-besselh(1,k0*Rbo)));

% PITCH MODE


%% DIFFRACTION problem
a0=1i*omega;
Bm(1)=(-1i*g*N_lambda_i(1)^.5)/(omega*cosh(k0*h)); Bm(2)=2*1i*Bm(1);

Fm(1)=-2*1i*Bm(1)/(pi*k0*Rbo*besselh(1,k0*Rbo));
Fm(2)=-2*1i*Bm(2)/(pi*k0*Rbo*(besselh(2,k0*Rbo)-besselh(1,k0*Rbo)/(k0*Rbo)));

% m=0
Hq7=zeros(Nn,2);
Hq7(:,1)=Fm(1)*Lni(:,1);

An=[0 conj(Lnj(:,1)');Lnj(:,1) Dm(:,:,1)]\[0;Hq7(:,1)];

Umi=(1/h)*(conj(Lni')*An(2:end)); Umj=(1/(h-b))*(Lnj(:,2:end)'*An(2:end));

Ami(:,1)=Umi./conj(lambda_i.*(-besselk(1,lambda_i.*Rbo)))';
Ami(1,1)=Ami(1,1) - Bm(1)*(besselj(1,k0*Rbo)/besselh(1,k0*Rbo));
Bmj(1,1)=An(1);% An=An(2:end);
Bmj(2:end,1)=Umj./(lambda_j.*(besseli(1,lambda_j.*Rbo)));

tmp1=sqrt(2)*Rbo.*(((-1).^j.*besseli(1,lambda_j.*Rbo))./(lambda_j.^2.*besseli(1,lambda_j.*Rbo)))'*Umj;
Fe(w,2)=2*pi*rho*a0*(.5*Bmj(1,1)*Rbo^2 + tmp1);

% m=1
Hq7(:,2)=Fm(2)*Lni(:,1);

An=Dm(:,:,2)\Hq7(:,2);
Umi=(1/h)*(conj(Lni')*An); Umj=(1/(h-b))*(Lnj'*An);

Fe(w,1) = -pi*a0*rho*Rbo*(Fm(2)*I1(1) + (I1.*Gmi(2,:))*Umi);

tmp1=.25*Umj(1)*Rbo^4 + sqrt(2)*Rbo^2*(((-1).^j.*besseli(2,lambda_j.*Rbo))./(lambda_j.^2.*(besseli(2,lambda_j.*Rbo) + besseli(1,lambda_j.*Rbo)./(lambda_j.*Rbo))))'*Umj(2:end);
Fe(w,3) = -pi*a0*rho*( Rbo*(Fm(2)*I5(1) + (I5.*Gmi(2,:))*Umi) + tmp1);

if isfield(options,'matchingConditions')
	if (options.matchingConditions==1)&&(omega==options.matchingFrequency)
		
		q=3;
		
		dz=.1; dr=.1;
		z = -10:dz:0;
		r = 0:dr:10;
		
% 		[R,Z]=meshgrid(r,z); [nz,nr]=size(R);
% 		
% 		phi=zeros(nz,nr); U=zeros(nz,nr); W=zeros(nz,nr);
% 		
% 		phi((Z>-b)&(R<Rbo))=NaN;
% 		U((Z>-b)&(R<Rbo))=NaN;
% 		W((Z>-b)&(R<Rbo))=NaN;
% 		
% 		phip=zeros(nz,nr); up=zeros(nz,nr); wp=zeros(nz,nr);
% 		% Solutions particulières fonction du problème de radiation
% 		if q==3% Heave
% 			m=0;
% 			% Domaine Omega_4
% % 			phip((R<=Rbo)&(Z<=-b))=((Z((R<=Rbo)&(R>=Rc)&(Z>=-e1)&(Z<=-b))+e1).^2-.5*R((R<=Rbo)&(R>=Rc)&(Z>=-e1)&(Z<=-b)).^2)./(2*(h-b));
% 			up((R<=Rbo)&(Z<=-b))=-.5*R((R<=Rbo)&(Z<=-b))/(h-b);
% 			wp((R<=Rbo)&(Z<=-b))=(Z((R<=Rbo)&(Z<=-b))+h)/(h-b);
% 		else
% 			m=1;
% 		end
		
		u_r = zeros(length(z),length(r));
		u_z = zeros(length(z),length(r));
		for i=1:length(z)
			for j=1:length(r)
				if ((z(i)>-b)&&(r(j)<Rbo))
					u_r(i,j) = NaN;% pour ne pas afficher de couleurs
					u_z(i,j) = NaN;
				elseif ((z(i) <= -b)&&(r(j) < Rbo))
					u_r(i,j) = abs(-.5*r(j)/(h-b) + sqrt(2)*(lambda_j.*(besseli(1,lambda_j*r(j))).*cos(lambda_j.*(z(i)+h)))'*bmj(2:end,2));
					u_z(i,j) = 1;%abs( (z(i)+h)/(h-b) - sqrt(2)*( lambda_j'.*(besseli(0,lambda_j'*r(j))./besseli(0,lambda_j'*Rbo)).*sin(lambda_j'.*(z(i)+h)) )*Bq3(2:end));
				elseif (r(j) >= Rbo)
					u_r(i,j) = abs(-(lambda_i.*besselk(1,lambda_i*r(j)).*N_lambda_i.^(-.5).*cos(lambda_i*(z(i)+h)))*ami(:,2));
					u_z(i,j) = 1;%abs((lambda_i.*besselk(0,lambda_i*r(j))./besselk(0,lambda_i*Rbo).*sqrt(1./N_lambda_i).*sin(lambda_i*(z(i)+h)) )*Aq3);
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

end% END OF FUNCTION


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



