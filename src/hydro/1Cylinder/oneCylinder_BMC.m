function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = oneCylinder_BMC(Omega, depth, WECstructure, options)%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                     BOTTOM MOUNTED CYLINDER
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the hydrodynamic coefficients for the bottom mounted
% cylinder.
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
%                          Description
%
%
%           |---->| (R)
%           .
%           |
%           .-----+
%           |ccccc|
%           .ccccc|
%           |ccccc|
% --------  +ccccc| ------------------------------------ (z=0)
%           |ccccc|      ^
%           .ccccc|      |
%           |ccccc|      |
%           .ccccc|      .
%           |ccccc|      .
%           .ccccc|      .
%           |ccccc|      |
%           .ccccc|      |
%   ________|ccccc+______v (h) _______
%														   
% 
% Written by S. OLAYA, IRDL/ENIB*, olaya@enib.fr
% * ENIB - École Nationale d'Ingénieurs de Brest, France
% Date : 2016/15/03
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

% ---- Parameters ---------------------------------------------------------
% Gravity
g = 9.81;
% Water density
rho = 1025;
% Water depth
h = depth;

% Cylinder radius
R = WECstructure.Rbo;


% ---- Some options -------------------------------------------------------
if nargin<4
	options = [];
end

% DATA sub-directory
if ~isfield(options,'DIR')
	DIR = pwd;
else
	DIR = options.DIR;
end

if ~isdir([DIR,filesep,'DATA'])
	mkdir([DIR,filesep,'DATA']);
	addpath(genpath(DIR));
	addpath(genpath([DIR,filesep,'DATA']));
end

if ~isfield(options,'Trunc')
	Ni = 50;
else
	Ni = options.Trunc.Ni;
end

if ~isfield(options,'Zc')
	Zc = 0;
else
	Zc = options.Zc;
end

% ---- Initialise output --------------------------------------------------
% Intrinsic matrix impedance -> A + (1i/omega)*B
Zr = zeros(3,3,length(Omega)+1);
% Added mass
A = zeros(3,3,length(Omega));
% Radiation damping
B = zeros(3,3,length(Omega));
% Wave excitation force computed from the diffraction problem
Fe = zeros(length(Omega)+1,3);
% Wave excitation force computed from the Haskind's relation
Fe_Haskind = zeros(length(Omega)+1,3); 
% Froude-Krilov (not implemented...)
Fe_FK = zeros(length(Omega)+1,3);

% ---- Some intermediate matrix definition --------------------------------
Hq1=zeros(Ni,1); Hq3=zeros(Ni,1); Hq5=zeros(Ni,1); Hq7 = zeros(Ni,2);
Dm0=zeros(Ni,Ni); Dm1=zeros(Ni,Ni);

% For each frequency, solve the radiation and the diffraction problem
% Add an extra value for the infinite problem
for w=1:length(Omega)+1

% Display where we are in computation
clc;
disp([num2str(round(w*100/(length(Omega)+1))),'%']);

% ---- Compute external eigenvalues ---------------------------------------
if w<=length(Omega)
	omega = Omega(w);
	lambda_i = fcn_waveDispersion(Omega(w),h,Ni);
	k0 = -imag(lambda_i(1));
else
	omega = Inf;
	i = 1:1:Ni;
	lambda_i = .5*(2*i-1)*pi/h;
end

% ---- Some intermediate definitions (need to be documented) --------------
N_lambda_i = 0.5.*(1+sin(2.*lambda_i.*h)./(2.*lambda_i.*h));

Delta_1i = -lambda_i.*(besselk(2,lambda_i.*R)./besselk(1,lambda_i.*R))+1/R;

% ---- Some pre-computed integrals ----------------------------------------
% Int_(-h)^(0){Zi(z)dz}
I1 = N_lambda_i.^(-.5).*sin(lambda_i.*h)./lambda_i;
% Int_(-h)^(0){(z-Zc)*Zi(z)dz}
I3 = N_lambda_i.^(-.5).*f2(lambda_i,-h,0,-Zc,h);

% ____ Intrinsic matrix definition ________________________________________
Dm1 = diag(Delta_1i);


%% ____ RADIATION PROBLEM _________________________________________________
% ---- SURGE MODE ---------------------------------------------------------
% Non-Homogeneous term for the surge problem
Hq1 = I1'/h;

% Unknown Fourier coefficients computation (most important)
% We are now able to evaluate the potential in the whole fluid domain
Aq1 = Dm1\Hq1;

Zr(1,1,w) = -pi*rho*R*(I1*Aq1);
Zr(3,1,w) = -pi*rho*R*(I3*Aq1);

Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1(1))...
                                  / (cosh(k0*h)*besselh(1,k0*R));


% ---- PITCH MODE ---------------------------------------------------------
% Non-Homogeneous term for the pitch problem
Hq5 = I3'/h;

% Unknown Fourier coefficients computation (most important)
% We are now able to evaluate the potential in the whole fluid domain
Aq5 = Dm1\Hq5;

Zr(3,3,w) = -pi*rho*R*(I3*Aq5);
Zr(1,3,w) = -pi*rho*R*(I1*Aq5);

Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5(1))...
                                    / (cosh(k0*h)*besselh(1,k0*R));


%% ____ DIFFRACTION PROBLEM _______________________________________________
% Some intermediate definitions
a0 = 1i*omega;
B0 = (-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));
B1 = 2*1i*B0;

% Non-Homogeneous term for the diffraction problem 
% (the incident potential is seen as non-homogeneous term)
Hq7 = -B1*(k0*besselj(0,k0*R)-besselj(1,k0*R)/R).*[1;zeros(Ni-1,1)];

% Unknown Fourier coefficients computation (most important)
% We are now able to evaluate the potential in the whole fluid domain
Aq7 = Dm1\Hq7;

Fe(w,1) = -pi*a0*rho*R*(B1*besselj(1,k0*R)*I1(1)+I1*Aq7);
Fe(w,3) = -pi*a0*rho*R*(B1*besselj(1,k0*R)*I3(1)+I3*Aq7);

end

%% ____ EXIT ______________________________________________________________
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

save([DIR,filesep,'DATA/hydroParameters.mat'],...
     'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');

end

%% ____ NESTED FUNCTIONS___________________________________________________
% ---- Integral -> Int_(a)^(b){(z+c)*cos(alpha*(z+d))} --------------------
function out = f2(alpha, a, b, c, d )
    epsilon = 1e-8;
    [n,m]=  size(alpha);
    out = zeros(n,m);
    for i = 1:length(alpha)
        if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
            out(i) = .5*((b+c)^2 - (a+c)^2);
        else
            out(i) = ((b+c).*sin(alpha(i).*(b+d)) - ...
                      (a+c).*sin(alpha(i).*(a+d)))./alpha(i)...
                   + (cos(alpha(i).*(b+d)) -...
                      cos(alpha(i).*(a+d)))./alpha(i).^2;
        end
    end
end





