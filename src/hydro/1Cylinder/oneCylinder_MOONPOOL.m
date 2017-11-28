% --> function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = oneCylinder_MOONPOOL(Omega, depth, WECstructure, options)
%
% Purpose :
% 
% Inputs :
% - Omega				: Vector of wave frequencies (rad/s)
% - depth				: Water depth (m), 0 for deep water 
% - WECstructure		: 
% - options				:
%		*
%		*
%		*
%
% Outputs :
% - Fe : Matrix length(Omega)x3 of exciation forces (complex values)
% - A  : Matrix 3x3xlength(Omega) of added mass coefficients
% - B  : Matrix 3x3xlength(Omega) of radiation damping coefficients
%
% Written by S. OLAYA, LBMS, olaya@enib.fr
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%								  |----------------------> (Rbo)
%								  |---> (Rbi)
%                                 .                                                              
% -------- +-----------------+    |    +-----------------+ ---------- (z=0)
%          |ccccccccccccccccc|    .    |ccccccccccccccccc|  ^     ^
%          |ccccccccccccccccc|    |    |ccccccccccccccccc|  |     |
%          |ccccccccccccccccc|    .    |ccccccccccccccccc|  |     |
%          +-----------------+    |    +-----------------+  v (b) |
%								  .                               .
%                                 |                               .
%																  .
%																  |
%															      v (h)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = oneCylinder_MOONPOOL(Omega, depth, WECstructure, options)
%%% Création d'un répertoire de données
DIR = options.DIR;
if ~isdir([DIR,filesep,'DATA'])
	mkdir([DIR,filesep,'DATA']);
	addpath(genpath([DIR,filesep,'DATA']));
end

g = 9.81;
rho = 1025;    
h = depth;

Rbi = WECstructure.Rbi;
Rbo = WECstructure.Rbo;

b = WECstructure.b;

Ni = options.Truncate.Ni;
Nn = options.Truncate.Nn;

Zc = options.Zc;

epsilon = 1e-5;

n = (1:1:Nn-1)';
alpha_n = (pi/(h-b)).*n;

[T_0n_prim_R1,T_0n_tild_prim_R1] = Tprim_func(0,Rbi,Rbo,Rbi,alpha_n);
[T_0n_prim_R2,T_0n_tild_prim_R2] = Tprim_func(0,Rbo,Rbo,Rbi,alpha_n);

[T_1n_prim_R1,T_1n_tild_prim_R1] = Tprim_func(1,Rbi,Rbo,Rbi,alpha_n);
[T_1n_prim_R2,T_1n_tild_prim_R2] = Tprim_func(1,Rbo,Rbo,Rbi,alpha_n);
% 
[T_0n_int,T_0n_tild_int] = Ti_func(0,Rbi,Rbo,alpha_n);
[T_1n_int,T_1n_tild_int] = Ti_func(1,Rbi,Rbo,alpha_n);

%%% Initialisation des sorties
Zr = zeros(3,3,length(Omega)); %%% matrice d'impédance --> A + (1i/omega)*B
A = zeros(3,3,length(Omega));
A_inf = zeros(3,3);
B = zeros(3,3,length(Omega));
Fe = zeros(length(Omega),3); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),3); %%% effort d'excitation
Fe_FK = zeros(length(Omega),3);%%% Froude-Krilov


for w=1:length(Omega)+1
	
	clc;
	disp([num2str(round(w*100/(length(Omega)+1))),'%']);
	
	if w<=length(Omega)
		omega = Omega(w);
% 		lambda_i = Lambda_i(w,1:Ni);
		lambda_i=fcn_waveDispersion(Omega(w),h,Ni);
		k0 = -imag(lambda_i(1));
	else
		omega = Inf;
		i=1:1:Ni;
		lambda_i = .5*(2*i-1)*pi/h;
	end

	N_lambda_i = .5.*(1 + sin(2.*lambda_i.*h)./(2.*lambda_i.*h));

	M_ntau = zeros(Nn,Ni);
	M_ntau(1,:) = N_lambda_i.^(-.5).*sin(lambda_i*(h-b))./( lambda_i*(h-b) );
	for i = 1 : Ni
		for p = 1 : Nn-1
			if (alpha_n(p)/lambda_i(i) < 1+epsilon)&&(alpha_n(p)/lambda_i(i) > 1-epsilon)
				M_ntau(p+1,i) = .5*sqrt(2).*N_lambda_i(i)^(-.5);
			else
				M_ntau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*...
					(sin((alpha_n(p)-lambda_i(i))*(h-b))/(alpha_n(p)-lambda_i(i))+...
					 sin((alpha_n(p)+lambda_i(i))*(h-b))/(alpha_n(p)+lambda_i(i)))/(h-b);
			end
		end
	end

	Gamma_0i = real(lambda_i.*(besseli(1,lambda_i.*Rbi)./besseli(0,lambda_i.*Rbi)));
	Gamma_1i = lambda_i.*(besseli(2,lambda_i.*Rbi)./besseli(1,lambda_i.*Rbi)) + 1/Rbi;
	
	Delta_0i = -lambda_i.*(besselk(1,lambda_i.*Rbo)./besselk(0,lambda_i.*Rbo));
	Delta_1i = -lambda_i.*(besselk(2,lambda_i.*Rbo)./besselk(1,lambda_i.*Rbo)) + 1/Rbo;

    I1 = N_lambda_i.^(-.5).*(sin(lambda_i.*h) - sin(lambda_i.*(h-b))) ./ lambda_i;			%%% Int_(-b)^(0){ Z_i(z) dz }
	I2 = N_lambda_i.^(-.5).*f2(lambda_i,-b,0,-Zc,h)';										%%% Int_(-b)^(0){ (z-Zc)*Z_i(z) dz}
	I3 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-b,h,h)';										%%% Int_(-h)^(-b){ (z+h)^2*Z_i(z) dz}
	
	%%% ***********************************************************************
    %%% Définition des matrices intrinsèques à la structure
	%%% ***********************************************************************
	D_m0_tau = zeros(Ni+Ni,Ni+Ni);
	%%% interface r=R2
    D_m0_tau(1:Ni,1:Ni) = diag(Delta_0i) - ((h-b)/h).*(M_ntau'*(diag(T_0n_prim_R2)*M_ntau));
    D_m0_tau(1:Ni,1+Ni:end) = -((h-b)/h).*(M_ntau'*(diag(T_0n_tild_prim_R2)*M_ntau));
	%%% interface r=R1
    D_m0_tau(1+Ni:end,1:Ni) = -((h-b)/h).*(M_ntau'*(diag(T_0n_prim_R1)*M_ntau));
    D_m0_tau(1+Ni:end,1+Ni:end) = diag(Gamma_0i) - ((h-b)/h).*(M_ntau'*(diag(T_0n_tild_prim_R1)*M_ntau));

	D_m1_tau = zeros(Ni+Ni,Ni+Ni);
	%%% interface r=R2
    D_m1_tau(1:Ni,1:Ni) = diag(Delta_1i) - ((h-b)/h).*(M_ntau'*(diag(T_1n_prim_R2)*M_ntau));
    D_m1_tau(1:Ni,1+Ni:end) = -((h-b)/h).*(M_ntau'*(diag(T_1n_tild_prim_R2)*M_ntau));
	%%% interface r=R1
    D_m1_tau(1+Ni:end,1:Ni) = -((h-b)/h).*(M_ntau'*(diag(T_1n_prim_R1)*M_ntau));
    D_m1_tau(1+Ni:end,1+Ni:end) = diag(Gamma_1i) - ((h-b)/h).*(M_ntau'*(diag(T_1n_tild_prim_R1)*M_ntau));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                  Problème de RADIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         SURGE MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
%%% ************************************************************************

	H_q1_tau = zeros(Ni+Ni,1);
	%%% Interface r = R2
    H_q1_tau(1:Ni,1) = (1/h).*I1';
	%%% Interface r = R1
    H_q1_tau(Ni+1:end,1) = (1/h).*I1';
    
    coeff = D_m1_tau\H_q1_tau;

    a_r1 = coeff(1:Ni,1);
    d_r1 = coeff(1+Ni:end,1);

    b_r1 = M_ntau*a_r1;
    c_r1 = M_ntau*d_r1;
	
	Zr(1,1,w) = -pi*rho*(Rbo*I1*a_r1 - Rbi*I1*d_r1);

	tmp1 = T_1n_int(1)*b_r1(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*b_r1(2:end));
	tmp2 = T_1n_tild_int(1)*c_r1(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*c_r1(2:end));
	Zr(3,1,w) = - pi*rho*(Rbo*I2*a_r1 - Rbi*I2*d_r1 + tmp1 + tmp2);

	% SURGE Wave excitation force --> Haskind's theorem
	Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*a_r1(1))/(cosh(k0*h)*besselh(1,k0*Rbo));
end%%% END OF SURGE

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%                         HEAVE MODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(options.dof,'all') || strcmp(options.dof,'heave')
%%% ************************************************************************
    % % solution particulière pour r = R2
    P_3n(1,1) = ( (h-b) / 12 ) * ( 2 - 3*(Rbo^2 / (h-b)^2) );
    P_3n(2:Nn,1) = (sqrt(2)*(h-b).*(-1).^n) ./ (n.*pi).^2;
	
    % % solution particulière pour r = R
    Q_3n(1,1) = ( (h-b) / 12 ) * ( 2 - 3*(Rbi^2 / (h-b)^2) );
    Q_3n(2:Nn,1) = (sqrt(2)*(h-b).*(-1).^n) ./ (n.*pi).^2;

	H_q3_tau = zeros(Ni+Ni,1);
	%%% Interface r = R2
    H_q3_tau(1:Ni) = -((h-b)/h).*(M_ntau'*(T_0n_prim_R2.*P_3n+T_0n_tild_prim_R2.*Q_3n)) - (.5*Rbo/h).*M_ntau(1,:)';
	%%% Interface r = R1
    H_q3_tau(1+Ni:end) = -((h-b)/h).*(M_ntau'*(T_0n_prim_R1.*P_3n+T_0n_tild_prim_R1.*Q_3n)) - (.5*Rbi/h).*M_ntau(1,:)';

    %%% Calcul des coefficients A_mi, C_mn
    coeff = D_m0_tau\H_q3_tau;

    a_r3 = coeff(1:Ni);
    d_r3 = coeff(Ni+1:end);
   
    b_r3 = M_ntau*a_r3 - P_3n;
    c_r3 = M_ntau*d_r3 - Q_3n;

	tmp1 = T_0n_int(1)*b_r3(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*b_r3(2:end));
	tmp2 = T_0n_tild_int(1)*c_r3(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*c_r3(2:end));
	tmp3 = .25*(h-b)*( (Rbo^2-Rbi^2) - .25*( (Rbo^4-Rbi^4)/(h-b)^2 )); 
	Zr(2,2,w) = 2*pi*rho*(tmp1 + tmp2 + tmp3);

	% HEAVE Wave excitation force --> Haskind's theorem
	Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*a_r3(1)) / (cosh(k0*h)*besselh(0,k0*Rbo));

end% *** END OF HEAVE

%%% ***********************************************************************
%%%                         PITCH MODE
if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
%%% ***********************************************************************

    P_5n(1,1) = (1/(24*(h-b))) * (3*Rbo^3 - 4*Rbo*(h-b)^2);
    P_5n(2:Nn,1) = -(sqrt(2)*Rbo*(h-b).*(-1).^n) ./ (n.*pi).^2;
    
    Q_5n(1,1) = (1/(24*(h-b))) * ( 3*Rbi^3 - 4*Rbi*(h-b)^2);
    Q_5n(2:Nn,1) = -(sqrt(2)*Rbi*(h-b).*(-1).^n) ./ (n.*pi).^2;
	
	H_q5_tau = zeros(Ni+Ni,1);
	%%% Interface r = R2
    H_q5_tau(1:Ni,1) = -((h-b)/h).*(M_ntau'*(T_1n_prim_R2.*P_5n + T_1n_tild_prim_R2.*Q_5n)) + ((3*Rbo^2)/(8*h)).*M_ntau(1,:)' - .5/(h*(h-b))*I3' + (1/h)*I2';
    %%% Interface r = R1
    H_q5_tau(1+Ni:end,1) = -((h-b)/h).*(M_ntau'*(T_1n_prim_R1.*P_5n + T_1n_tild_prim_R1.*Q_5n)) + ((3*Rbi^2)/(8*h)).*M_ntau(1,:)' - .5/(h*(h-b))*I3' + (1/h)*I2';
                             
    coeff = D_m1_tau\H_q5_tau;

    a_r5 = coeff(1:Ni,1);
    d_r5 = coeff(1+Ni:end,1);
 
    b_r5 = M_ntau*a_r5 - P_5n;
    c_r5 = M_ntau*d_r5 - Q_5n;

	tmp1 = T_1n_int(1)*b_r5(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*b_r5(2:end));
	tmp2 = T_1n_tild_int(1)*c_r5(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*c_r5(2:end));
	tmp3 = ( (Rbo^6-Rbi^6)/6 - (Rbo^4-Rbi^4)*(h-b)^2 ) / (8*(h-b));
	Zr(3,3,w) = - pi*rho*(Rbo*I2*a_r5 - Rbi*I2*d_r5 + tmp1 + tmp2 + tmp3);

	Zr(1,3,w) = -pi*rho*(Rbo*I1*a_r5 - Rbi*I1*d_r5);

	% PITCH Wave excitation force --> Haskind's theorem
	Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*a_r5(1))/(cosh(k0*h)*besselh(1,k0*Rbo));

end% *** END OF PITCH

% %%% ***********************************************************************
% %%% ***********************************************************************
% %%%             Problème de DIFFRACTION
% %%% ***********************************************************************
% %%% ***********************************************************************

a0 = 1i*omega;
B0 = (-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));
B1 = 2*1i*B0;

%%% Wave Excitation force -- FROUDE-KRILOV

if strcmp(options.dof,'all') || strcmp(options.dof,'heave')

	H_d_0tau = zeros(Ni+Ni,1);
	H_d_0tau(1:Ni) = B0.*( besselj(0,k0*Rbo).*((h-b)/h).*(M_ntau'*(T_0n_prim_R2.*M_ntau(:,1))) + k0*besselj(1,k0*Rbo)*[1;zeros(Ni-1,1)]);
	H_d_0tau(1+Ni:end) = B0*besselj(0,k0*Rbo).*((h-b)/h).*(M_ntau'*(T_0n_prim_R1.*M_ntau(:,1)));

	coeff = D_m0_tau\H_d_0tau;

	a_0i = coeff(1:Ni);
	d_0i = coeff(1+Ni:end);

	b_0n = M_ntau*a_0i + B0.*besselj(0,k0*Rbo).*M_ntau(:,1);
	c_0n = M_ntau*d_0i;    

	tmp1 = T_0n_int(1)*b_0n(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*b_0n(2:end));
	tmp2 = T_0n_tild_int(1)*c_0n(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*c_0n(2:end));
	Fe(w,2) = 2*pi*a0*rho*(tmp1 + tmp2);

	Fe_FK(w,2) = 2*pi*a0*rho*B0*N_lambda_i(1)^(-.5)*cosh(k0*(h-b))*(Rbo*besselj(1,k0*Rbo) - Rbi*besselj(1,k0*Rbi))/k0;
end

if ~strcmp(options.dof,'heave')%%% tout sauf heave
	H_d_1tau = zeros(Ni+Ni,1);
	H_d_1tau(1:Ni) = B1.*( besselj(1,k0*Rbo).*((h-b)/h).*(M_ntau'*(T_1n_prim_R2.*M_ntau(:,1))) - (k0*besselj(0,k0*Rbo)-besselj(1,k0*Rbo)/Rbo)*[1;zeros(Ni-1,1)]);
	H_d_1tau(1+Ni:end) = B1*besselj(1,k0*Rbo).*((h-b)/h).*(M_ntau'*(T_1n_prim_R1.*M_ntau(:,1)));

	coeff = D_m1_tau\H_d_1tau;

	a_1i = coeff(1:Ni);
	d_1i = coeff(1+Ni:end);

	b_1n = M_ntau*a_1i + B1.*besselj(1,k0*Rbo).*M_ntau(:,1);
	c_1n = M_ntau*d_1i;      
end

if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
	Fe(w,1) = - pi*a0*rho*(Rbo*(B1*besselj(1,k0*Rbo)*I1(1) + I1*a_1i) - Rbi*I1*d_1i);

	Fe_FK(w,1) = -pi*a0*rho*B1*(Rbo*besselj(1,k0*Rbo)*I1(1) - Rbi*besselj(1,k0*Rbi)*I1(1));

end

if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')		
    tmp1 = T_1n_int(1)*b_1n(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*b_1n(2:end));
    tmp2 = T_1n_tild_int(1)*c_1n(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*c_1n(2:end));
	Fe(w,3) = -pi*a0*rho*( Rbo*(B1*besselj(1,k0*Rbo)*I2(1) + I2*a_1i) - Rbi*I2*d_1i + tmp1 + tmp2);

	tmp = N_lambda_i(1)^(-.5)*f2(lambda_i(1),-b,0,-Zc,h);%% int_(-b)^(0){ (z-Zc)Z_i(z) dz}
	tmp1 = Rbo*besselj(1,k0*Rbo)*tmp - Rbi*besselj(1,k0*Rbi)*tmp;
	tmp2 = N_lambda_i(1)^(-.5)*cosh(k0*(h-b))*(Rbo^2*besselj(2,k0*Rbo) - Rbi^2*besselj(2,k0*Rbi))/k0;
	Fe_FK(w,3) = -pi*a0*rho*B1*(tmp1 + tmp2);
end


end% END OF FOR() LOOP

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%								EXIT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% *** ATTENTION -> ici on reprend une dépendance en temps positive exp(-iwt) --> exp(iwt)
for i=1:size(A,1)
	for j=1:size(A,1)
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
epsilon = 1e-5;
for i=1:length(alpha)
	if (abs(alpha(i)) < epsilon)&&(abs(alpha(i)) > -epsilon)%alpha == 0
		out(i,1) = .5*((b+c)^2 - (a+c)^2);
	else
		out(i,1) = ((b+c).*sin(alpha(i).*(b+d)) - (a+c).*sin(alpha(i).*(a+d)))./alpha(i)  + (cos(alpha(i).*(b+d)) - cos(alpha(i).*(a+d)))./alpha(i).^2;
	end
end
end

function [T_prim,T_tild_prim] = Tprim_func(m,r,a,b,alpha)
	T_prim = zeros(length(alpha)+1,1);
	T_tild_prim = zeros(length(alpha)+1,1);
	
	cst = besseli(m,alpha.*a).*besselk(m,alpha.*b) - besseli(m,alpha.*b).*besselk(m,alpha.*a);
	switch m
		case 0
			T_prim(1,1) = (r*log(a/b))^(-1);
			T_tild_prim(1,1) = -(r*log(a/b))^(-1);

		case 1
			T_prim(1,1) = ((1/b + b/r^2)/(a/b-b/a));
			T_tild_prim(1,1) = ((-a/r^2 - 1/a)/(a/b-b/a));
	end
	
	T_prim(2:end,1) = (alpha./cst).*( besselk(m,alpha.*b).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r)) +...
								besseli(m,alpha.*b).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));

	T_tild_prim(2:end,1) = -(alpha./cst).*( besselk(m,alpha.*a).*(besseli(m+1,alpha.*r) + (m./(alpha.*r)).*besseli(m,alpha.*r)) +...
								besseli(m,alpha.*a).*(besselk(m+1,alpha.*r) - (m./(alpha.*r)).*besselk(m,alpha.*r)));
end