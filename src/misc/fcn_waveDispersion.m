function Lambda = fcn_waveDispersion(w, depth, N, opts)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description : This code has been partially developed during my PhD at 
% the Institut de Recherche Dupuy de Lôme (http://www.irdl.fr/) where I
% worked on Self-Reacting Point Absorber control.
% PhD Supervisors:
%	- Prof. Mohamed El-Hachemi BENBOUZID - Mohamed.Benbouzid@univ-brest.fr
%	- Dr Jean-Matthieu BOURGEOT - bourgeot@enib.fr
% 
% Purpose : Compute the requiered eigenvalues from the dispersion relation 
% based on the article of  Newman, J. N. (1990). 
% "Numerical solutions of the water-wave dispersion relation"
% Applied Ocean Research, 12(1), 14–18. 
% http://doi.org/10.1016/S0141-1187(05)80013-6
% 
% Inputs :
% - w		: Wave frequencie (rad/s)
% - depth	: Water depth (m)
% - N		: Number of eigenvalues (infinite truncation)
% - options	: Some options for Newton-Raphson Algorithm
%		* kmax : Max iteration
%		* eps  : Tolerance
%
% Outputs :
% - Lambda : Eigenvalues
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

if nargin<4
	kmax = 25;
	eps = 1e-6;
else
	kmax = opts.kmax;
	eps = opts.eps;
end

g = 9.81;
cst = w^2*depth/g;

f1 = @(x,cst) .5*log((x+cst)/(x-cst))-x;
f1p = @(x,cst) -cst/(x^2-cst^2) -1;

f2 = @(x,n,cst) atan(cst/(n*pi-x))-x;
f2p = @(x,n,cst) cst/(cst^2+(n*pi-x)^2)-1;
f2y = @(x,n) (n*pi-x)/depth;

%% ____ Compute the solution of x=y*tanh(y) _______________________________
k0 = 0;
if cst<=2
	y = sqrt(cst)*(.9994+.1701*cst+.0305*cst^2);
else
	y = cst+2*cst*exp(-2*cst)-6*cst^2*exp(-4*cst);
end
if abs(y-cst)>eps
k = 1;
while k<kmax
	yprev = y;
	y = y-f1(y,cst)/f1p(y,cst);
	err = abs(yprev-y);
	if err < eps
		k0 = y/depth;
		break;
	end
	k = k+1;
end
else
	k0 = y/depth;
end

%% ____ Compute the solution of x=y*tan(y) ________________________________
ki = zeros(1,N-1);
u = 3*cst/(7+3*cst);
d = .0159+.1032*u+4.3152*u^2-2.8768*u^3;%
for i = 1:N-1
	y = f2y(d,i);
	k = 1;
	while k<kmax
		d = d-f2(d,i,cst)/f2p(d,i,cst);
		yprev = y;
		y = f2y(d,i);
		err = abs(yprev-y);
		if err < eps
			ki(1,i) = y;
			break;
		end
		k = k+1;
	end
	d = d-(pi*cst)/(cst^2+pi*(i+1)*(pi*i-d));
end

Lambda = [-1i*k0,ki];

end

