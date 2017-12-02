% --> function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = twoCylinders_COLUMN(Omega, depth, WECstructure, options)
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
% - Fe : Matrix length(Omega)x(3xnBodies) of exciation forces (complex values)
% - A  : Matrix (3x2)x(3x2)xlength(Omega) of added mass coefficients
% - B  : Matrix (3x2)x(3x2)xlength(Omega) of radiation damping coefficients
%
% Written by S. OLAYA, LBMS, olaya@enib.fr
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%           |-----------------------> (Rbo)
%           |-----> (Rbi=Rc)
%           |-----> (Rc)
%           .                                                              
% --------  +-----++----------------+ --------------------------- (z=0)
%           .ooooo||cccccccccccccccc|  ^      ^      ^
%           |ooooo||cccccccccccccccc|  |      |      |
%           .ooooo||cccccccccccccccc|  |      |      |
%           |ooooo|+----------------+  v (b)  |      |
%           .ooooo|           ^               |      |
%           |ooooo|           |               |      |
%           .ooooo|           |               |      |
%           |ooooo|           |               |      |
%           .ooooo|          (1)              |      |
%           |ooooo|                           |      |
%           .ooooo|                           |      |
%           |ooooo|                           |      |
%           .ooooo|                           |      |
%           |ooooo|                           |      |
%           .ooooo|<---(2)                    |      |
%           |ooooo|                           |      |
%           .ooooo|                           |      |
%           |ooooo|                           |      |
%           .ooooo|                           |      |
%           |ooooo|                           |      |
%           +-----+                           v (e2) .
%           .                                        .
%           |                                        .
%                                                    |
%                                                    |
%                                                    |
%                                                    v (h)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Fe, A, B, A_inf, Fe_Haskind, Fe_FK] = twoCylinders_COLUMN(Omega, depth, WECstructure, options)

%%% Création d'un répertoire de données
DIR = options.DIR;
if ~isdir([DIR,filesep,'DATA'])
	mkdir([DIR,filesep,'DATA']);
	addpath(genpath([DIR,filesep,'DATA']));
end

g = 9.81;
rho = 1025;    
h = depth;
    
Rbo = WECstructure.Rbo;
Rc = WECstructure.Rc;

b = WECstructure.b;
e2 = WECstructure.e2;

Ni = options.Truncate.Ni;
Nn = options.Truncate.Nn;
Nj = options.Truncate.Nj;

epsilon = 1e-7;

Zc = [0;0];%options.Zc;

%%% Initialisation des sorties
Zr = zeros(6,6,length(Omega)); %%% matrice d'impédance --> A + (1i/omega)*B
A = zeros(6,6,length(Omega));
B = zeros(6,6,length(Omega));
Fe = zeros(length(Omega),6); %%% force d'excitation
Fe_Haskind = zeros(length(Omega),6); %%% effort d'excitation
Fe_FK = zeros(length(Omega),6);%%% Froude-Krilov	

n = (1:1:Nn-1)'; lambda_n = (pi/(h-b)).*n;
j = (1:1:Nj-1)'; lambda_j = (pi/(h-e2)).*j;

Hq1=zeros(Ni+Nn,2);
Hq3=zeros(Ni+Nn,2);
Hq5=zeros(Ni+Nn,2);
Hq7 = zeros(Ni+Nn,2);

Dm0=zeros(Ni+Nn,Ni+Nn);
Dm1=zeros(Ni+Nn,Ni+Nn);

%% Calcul des solutions particulières
% *** solution particulière pour q=3 (Heave)
P3n=zeros(Nn,2);
P3n(1,1) = -( (h-b) / 12 ) * ( 2 - 3*(Rbo^2 / (h-b)^2) );
P3n(2:end,1) = -(sqrt(2)*(h-b).*(-1).^n) ./ (n.^2 .* pi^2);

P3j=zeros(Nj,2);
%k=1
P3j(1,1) = ( (h-e2)^2 / (12*(h-b)) ) * ( 2 - 3*(Rc^2 / (h-e2)^2 ));
P3j(2:end,1) = (sqrt(2)*(h-e2)^2.*(-1).^j) ./ ((j.^2 .* pi^2).*(h-b));
% k=2
P3j(1,2)=-( (h-e2) / 12 ) * ( 2 - 3*(Rc^2 / (h-e2)^2) );
P3j(2:end,2) = -(sqrt(2)*(h-e2).*(-1).^j) ./ (j.^2 .* pi^2);

% *** solution particulière pour q=5 (Pitch)
P5n=zeros(Nn,2);
P5n(1,1)=-(Rbo^3/(24*(h-b))) * ( 3 - 4*((h-b)/Rbo)^2 );
P5n(2:end,1)=((sqrt(2)*Rbo)/(2*(h-b)^2)).*f1(lambda_n,-h,-b,h,h);

% % solution particulière pour r = Rc
P5j=zeros(Nj,2);
% k=1
P5j(1,1)= (Rc^3/(24*(h-b))) * ( 3 - 4*((h-e2)/Rc)^2);
P5j(2:end,1)= -((sqrt(2)*Rc)/(2*(h-b)*(h-e2))).*f1(lambda_j,-h,-e2,h,h);
% k=2
P5j(1,2)=-(Rc^3/(24*(h-e2))) * ( 3 - 4*((h-e2)/Rc)^2);
P5j(2:end,2)=((sqrt(2)*Rc)/(2*(h-e2)^2)).*f1(lambda_j,-h,-e2,h,h);


%%
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

M_ntau = zeros(Nn,Ni);
M_ntau(1,:) = N_lambda_i.^(-.5).*sin(lambda_i.*(h-b))./( (h-b).*lambda_i );
for i=1:Ni
	for p=1:Nn-1
		if (lambda_n(p)/lambda_i(i) < 1+epsilon)&&(lambda_n(p)/lambda_i(i) > 1-epsilon)
			M_ntau(p+1,i) = .5*sqrt(2).*N_lambda_i(i)^(-.5);
		else
			M_ntau(p+1,i) = .5*sqrt(2)*N_lambda_i(i)^(-.5)*(sin((lambda_n(p)-lambda_i(i))*(h-b))/(lambda_n(p)-lambda_i(i))...
				+ sin((lambda_n(p)+lambda_i(i))*(h-b))/(lambda_n(p)+lambda_i(i)))/(h-b);
		end
	end
end

% L_jtau = (integraleInterface([0;lambda_j],[0;lambda_n],h,e2,h,h,epsilon))./(h-e2);
% L_jtau(2:end,1) = sqrt(2).*L_jtau(2:end,1);
% L_jtau(1,2:end) = sqrt(2).*L_jtau(1,2:end);
% L_jtau(2:end,2:end) = 2.*L_jtau(2:end,2:end);
	L_jtau=zeros(Nj,Nn);
	L_jtau(1,1)=1;
	L_jtau(1,2:end)=sqrt(2).*sin(lambda_n.*(h-e2))./((h-e2).*lambda_n);
	for i=1:Nn-1
		for p=1:Nj-1
		    if (lambda_j(p)/lambda_n(i) < 1+epsilon)&&(lambda_j(p)/lambda_n(i) > 1-epsilon)
		        L_jtau(p+1,i+1)=1+sin(2*lambda_j(i).*(h-e2))./(2*lambda_j(i).*(h-e2));
		    else
		        L_jtau(p+1,i+1)= ...
					( sin((lambda_j(p)-lambda_n(i))*(h-e2))/(lambda_j(p)-lambda_n(i))...
					+ sin((lambda_j(p)+lambda_n(i))*(h-e2))/(lambda_j(p)+lambda_n(i)))/(h-e2);
		    end
		end
	end

Gamma_0j(1,1)=0;
Gamma_0j(2:Nj,1) = lambda_j.*(besseli(1,lambda_j.*Rc) ./ besseli(0,lambda_j.*Rc));

Gamma_1j(1,1)=1/Rc;
Gamma_1j(2:Nj,1) = lambda_j.*(besseli(0,lambda_j.*Rc) ./ besseli(1,lambda_j.*Rc)) - 1/Rc;

Delta_0i(1,1:Ni) = -lambda_i.*(besselk(1,lambda_i.*Rbo) ./ besselk(0,lambda_i.*Rbo));
Delta_1i(1,1:Ni) = -lambda_i.*(besselk(2,lambda_i.*Rbo)./besselk(1,lambda_i.*Rbo)) + 1/Rbo;

[Tp1_0n_Rb,Tp2_0n_Rb] = Tp_func(0,Rbo,Rc,Rbo,[0;lambda_n]);
[Tp1_0n_Rc,Tp2_0n_Rc] = Tp_func(0,Rc,Rc,Rbo,[0;lambda_n]);

[Tp1_1n_Rb,Tp2_1n_Rb] = Tp_func(1,Rbo,Rc,Rbo,[0;lambda_n]);
[Tp1_1n_Rc,Tp2_1n_Rc] = Tp_func(1,Rc,Rc,Rbo,[0;lambda_n]);

[Ti1_0n,Ti2_0n] = Ti_func(0,Rc,Rbo,[0;lambda_n]);
[Ti1_1n,Ti2_1n] = Ti_func(1,Rc,Rbo,[0;lambda_n]);

I1 = N_lambda_i.^(-.5).*( sin(lambda_i.*h)-sin(lambda_i.*(h-b)) ) ./ lambda_i;              % % Int_(-b)^(0){ Z_i^I }
I2 = [(e2-b);sqrt(2).*(sin(lambda_n.*(h-b)) - sin(lambda_n.*(h-e2)))./lambda_n];               % % Int_(-e2)^(-e1){ Z_n^II } 
I3 = N_lambda_i.^(-.5).*f2(lambda_i,-b,0,Zc(1),h);                                      	% % Int_(-b)^(0){ (z-Zc)*Z_i^I }
I4 = [.5*((b+Zc(2))^2-(e2+Zc(2))^2);sqrt(2).*f2(lambda_n,-e2,-b,Zc(2),h)];                   % % Int_(-e1)^(-b){ (z-Zc)*Z_n^II }  
I6 = [(h-b)^3/3;sqrt(2).*f1(lambda_n,-h,-b,h,h)];                                            % % Int_{-h}^{-b}{(z+h)^2*Z_n^II}
I7 = N_lambda_i.^(-.5).*f1(lambda_i,-h,-b,h,h);                                             % % Int_(-h)^(-b){ (z+h)^2*Z_i^I }
I8 = [(h-e2)^3/3;sqrt(2).*f1(lambda_n,-h,-e2,h,h)];                                          % % Int_{-h}^{-e1}{(z+h)^2*Z_n^II}

%% Définition des matrices intrinsèques à la structure
Dm0(1:Ni,1:Ni)=diag(Delta_0i)-((h-b)/h).*(M_ntau'*(diag(Tp1_0n_Rb)*M_ntau));
Dm0(1:Ni,Ni+1:end)=-((h-b)/h).*(diag(Tp2_0n_Rb)*M_ntau)';
Dm0(Ni+1:end,1:Ni)=diag(Tp1_0n_Rc)*M_ntau ;
Dm0(Ni+1:end,Ni+1:end)= diag(Tp2_0n_Rc) - ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_0j)*L_jtau));

Dm1(1:Ni,1:Ni)=diag(Delta_1i) - ((h-b)/h).*(M_ntau'*(diag(Tp1_1n_Rb)*M_ntau));
Dm1(1:Ni,Ni+1:end)=-((h-b)/h).*(diag(Tp2_1n_Rb)*M_ntau)';
Dm1(Ni+1:end,1:Ni)=diag(Tp1_1n_Rc)*M_ntau ;
Dm1(Ni+1:end,Ni+1:end) = diag(Tp2_1n_Rc) - ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_1j)*L_jtau));


%% Problème de RADIATION
%% SURGE MODE
% *** Cas k=1, Buoy oscillating, column fixed
Hq1(1:Ni,1)=(1/h).*I1';
Hq1(Ni+1:end,1)=zeros(Nn,1);
              
coeff = Dm1\Hq1(:,1);

Aq1_k1 = coeff(1:Ni);
Cq1_k1 = coeff(Ni+1:end);

Bq1_k1 = M_ntau*Aq1_k1;
Fq1_k1 = L_jtau*Cq1_k1;

Zr(1,1,w) = - pi*rho*Rbo*(I1*Aq1_k1);
Zr(4,1,w) = - pi*rho*Rc*(I2'*Cq1_k1);

tmp1 = Ti1_1n(1)*Bq1_k1(1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Bq1_k1(2:end));
tmp2 = Ti2_1n(1)*Cq1_k1(1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Cq1_k1(2:end));
Zr(3,1,w) = - pi*rho*(Rbo*I3*Aq1_k1 + tmp1 + tmp2 );

tmp1 = .25*Fq1_k1(1)*Rc^3 + sqrt(2)*Rc^2*(((-1).^j.*besseli(2,lambda_j.*Rc))./(lambda_j.*besseli(1,lambda_j.*Rc)))'*Fq1_k1(2:end);   
Zr(6,1,w) = - pi*rho*( Rc*I4'*Cq1_k1 + tmp1);
	
% *** Cas k=2, Buoy fixed , column oscillating 
Hq1(1:Ni,2) = zeros(Ni,1);
Hq1(Ni+1:end,2)=(1/(h-b)).*I2; 

coeff = Dm1\Hq1(:,2);

Aq1_k2 = coeff(1:Ni);
Cq1_k2 = coeff(Ni+1:end);

Bq1_k2 = M_ntau*Aq1_k2;
Fq1_k2 = L_jtau*Cq1_k2;

Zr(1,4,w) = - pi*rho*Rbo*(I1*Aq1_k2);
Zr(4,4,w) = - pi*rho*Rc*(I2'*Cq1_k2);

tmp1 = Ti1_1n(1)*Bq1_k2(1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Bq1_k2(2:end));
tmp2 = Ti2_1n(1)*Cq1_k2(1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Cq1_k2(2:end));
Zr(3,4,w) = - pi*rho*(Rbo*I3*Aq1_k2 + tmp1 + tmp2 );

tmp1 = .25*Fq1_k2(1)*Rc^3 + sqrt(2)*Rc^2*(((-1).^j.*besseli(2,lambda_j.*Rc))./(lambda_j.*besseli(1,lambda_j.*Rc)))'*Fq1_k2(2:end);   
Zr(6,4,w) = - pi*rho*( Rc*I4'*Cq1_k2 + tmp1);

Fe_Haskind(w,1) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1_k1(1)) / (cosh(k0*h)*besselh(1,k0*Rbo));
Fe_Haskind(w,4) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq1_k2(1)) / (cosh(k0*h)*besselh(1,k0*Rbo));
	

%% HEAVE MODE
% ***     Cas k=1, Buoy oscillating, column fixed 
Hq3(1:Ni,1)=((h-b)/h).*M_ntau'*(Tp1_0n_Rb.*P3n(:,1)) -  (.5*Rbo/h).*M_ntau(1,:)';
Hq3(Ni+1:end,1)=((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*P3j(:,1))) - Tp1_0n_Rc.*P3n(:,1) + [(.5*Rc/(h-b));zeros(Nn-1,1)]; 

coeff = Dm0\Hq3(:,1);

Aq3_k1 = coeff(1:Ni);
Cq3_k1 = coeff(Ni+1:end);

Bq3_k1 = M_ntau*Aq3_k1 + P3n(:,1);
Fq3_k1 = L_jtau*Cq3_k1 + P3j(:,1);

tmp1 = Ti1_0n(1)*Bq3_k1(1) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Bq3_k1(2:end));
tmp2 = Ti2_0n(1)*Cq3_k1(1) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Cq3_k1(2:end));
tmp3 = .25*(h-b)*( (Rbo^2-Rc^2) - (1/4)*( (Rbo^4-Rc^4)/(h-b)^2 ) );
Zr(2,2,w) = 2*pi*rho*( tmp1 + tmp2 + tmp3 );

tmp1 = .5*Rc^2*Fq3_k1(1) + sqrt(2)*Rc*((-1).^j.*besseli(1,lambda_j.*Rc)./(lambda_j.*besseli(0,lambda_j.*Rc)))'*Fq3_k1(2:end);
Zr(5,2,w) = 2*pi*rho*tmp1;

% *** Cas k=2, Buoy fixed , column oscillating 
Hq3(1:Ni,2) = zeros(Ni,1);
Hq3(Ni+1:end,2)=((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*P3j(:,2))) - (.5*Rc/(h-b)).*L_jtau(1,:)'; 

coeff = Dm0\Hq3(:,2);

Aq3_k2 = coeff(1:Ni);
Cq3_k2 = coeff(Ni+1:end);

Bq3_k2 = M_ntau*Aq3_k2;
Fq3_k2 = L_jtau*Cq3_k2 + P3j(:,2);

tmp1 = Ti1_0n(1)*Bq3_k2(1) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Bq3_k2(2:end));
tmp2 = Ti2_0n(1)*Cq3_k2(1) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Cq3_k2(2:end));
Zr(2,5,w) = 2*pi*rho*(tmp1+tmp2);

tmp1 = .5*Rc^2*Fq3_k2(1) + sqrt(2)*Rc*((-1).^j.*besseli(1,lambda_j.*Rc)./(lambda_j.*besseli(0,lambda_j.*Rc)))'*Fq3_k2(2:end);
tmp2 = .5*(h-e2)*(.5*Rc^2 - (1/8)*( Rc^4/(h-e2)^2 ) );
Zr(5,5,w) = 2*pi*rho*( tmp1 + tmp2 );

Fe_Haskind(w,2) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3_k1(1)) / (cosh(k0*h)*besselh(0,k0*Rbo));
Fe_Haskind(w,5) = -(4*1i*rho*g*h*sqrt(N_lambda_i(1))*Aq3_k2(1)) / (cosh(k0*h)*besselh(0,k0*Rbo));
    
%%                         PITCH MODE

% *** Cas k=1, Buoy oscillating, plate fixed
Hq5(1:Ni,1)=((h-b)/h).*M_ntau'*(Tp1_1n_Rb.*P5n(:,1)) + ((3*Rbo^2)/(8*h)).*M_ntau(1,:)' - (1/(2*h*(h-b))).*I7' + (1/h).*I3';
Hq5(Ni+1:end,1)= ((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*P5j(:,1))) - Tp1_1n_Rc.*P5n(:,1) - [3*Rc^2/(8*(h-b));zeros(Nn-1,1)] + (1/(2*(h-b)^2)).*I6;

coeff=Dm1\Hq5(:,1);

Aq5_k1=coeff(1:Ni);
Cq5_k1=coeff(Ni+1:end);

Bq5_k1=M_ntau*Aq5_k1+P5n(:,1);
Fq5_k1=L_jtau*Cq5_k1+P5j(:,1);

tmp1 = Ti1_1n(1)*Bq5_k1(1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Bq5_k1(2:end));
tmp2 = Ti2_1n(1)*Cq5_k1(1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Cq5_k1(2:end));
tmp3 = ( (Rbo^6-Rc^6)/6 - (Rbo^4-Rc^4)*(h-b)^2 ) / (8*(h-b));
Zr(3,3,w) = - pi*rho*(Rbo*I3*Aq5_k1 + tmp1 + tmp2 + tmp3 );

tmp1 = (.5*Rc^3*(b^2-e2^2) - 4*Rc*( e2*(h-e2)^3 - b*(h-b)^3 - ((h-b)^4-(h-e2)^4)/4 )/3 ) / (8*(h-b));
tmp2 = .25*Fq5_k1(1)*Rc^3 + sqrt(2)*Rc^2*(((-1).^j.*besseli(2,lambda_j.*Rc))./(lambda_j.*besseli(1,lambda_j.*Rc)))'*Fq5_k1(2:end);
Zr(6,3,w) = - pi*rho*( Rc*(I4'*Cq5_k1 + tmp1) + tmp2 );

Zr(1,3,w) = - pi*rho*Rbo*(I1*Aq5_k1);

tmp1 = (Rc^3*(e2-b) - 4*Rc*((h-b)^3-(h-e2)^3)/3)/(8*(h-b));
Zr(4,3,w) = - pi*rho*Rc*(I2'*Cq5_k1 + tmp1);

% *** Cas k=2, Buoy fixed , column oscillating 

Hq5(1:Ni,2) = zeros(Ni,1);
Hq5(Ni+1:end,2)=((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*P5j(:,2))) + ((3*Rc^2)/(8*(h-b))).*L_jtau(1,:)' - (1/(2*(h-e2)*(h-b))).*I8 + (1/(h-b)).*I4 ;

coeff = Dm1\Hq5(:,2);

Aq5_k2 = coeff(1:Ni);
Cq5_k2 = coeff(Ni+1:end);

Bq5_k2 = M_ntau*Aq5_k2;
Fq5_k2 = L_jtau*Cq5_k2+P5j(:,2);

tmp1 = Ti1_1n(1)*Bq5_k2(1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Bq5_k2(2:end));
tmp2 = Ti2_1n(1)*Cq5_k2(1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Cq5_k2(2:end));
Zr(3,6,w) = - pi*rho*( Rbo*I3*Aq5_k2 + tmp1 + tmp2 );

tmp1 = .25*Fq5_k2(1)*Rc^3 + sqrt(2)*Rc^2*(((-1).^j.*besseli(2,lambda_j.*Rc))./(lambda_j.*besseli(1,lambda_j.*Rc)))'*Fq5_k2(2:end);
tmp2 = (Rc^6/6 - Rc^4*(h-e2)^2) / (8*(h-e2));
Zr(6,6,w) = - pi*rho*( Rc*I4'*Cq5_k2 + tmp1 + tmp2 );

Zr(1,6,w) = - pi*rho*Rbo*(I1*Aq5_k2);
Zr(4,6,w) = - pi*rho*Rc*(I2'*Cq5_k2);

Fe_Haskind(w,3) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5_k1(1)) / (cosh(k0*h)*besselh(1,k0*Rbo));
Fe_Haskind(w,6) = (4*rho*g*h*sqrt(N_lambda_i(1))*Aq5_k2(1)) / (cosh(k0*h)*besselh(1,k0*Rbo));

%% Problème de DIFFRACTION
a0 = 1i*omega;
%% m=0
B0=(-1i*g*N_lambda_i(1)^.5) / (omega*cosh(k0*h));

Hq7(1:Ni,1)=B0.*(besselj(0,k0*Rbo)*((h-b)/h).*(M_ntau'*(Tp1_0n_Rb.*M_ntau(:,1))) + k0*besselj(1,k0*Rbo).*[1;zeros(Ni-1,1)]);
Hq7(Ni+1:end,1)=-B0*besselj(0,k0*Rbo).*(Tp1_0n_Rc.*M_ntau(:,1));

coeff = Dm0\Hq7(:,1);

Aq7_m0 = coeff(1:Ni);
Cq7_m0 = coeff(Ni+1:end);

Bq7_m0 = M_ntau*Aq7_m0 + B0.*besselj(0,k0*Rbo).*M_ntau(:,1);
Fq7_m0 = L_jtau*Cq7_m0;

tmp1 = Ti1_0n(1)*Bq7_m0(1) + sqrt(2)*sum((-1).^n.*Ti1_0n(2:end).*Bq7_m0(2:end));
tmp2 = Ti2_0n(1)*Cq7_m0(1) + sqrt(2)*sum((-1).^n.*Ti2_0n(2:end).*Cq7_m0(2:end));
Fe(w,2) = 2*pi*a0*rho*(tmp1+tmp2);

tmp1 = .5*Rc^2*Fq7_m0(1) + sqrt(2)*Rc*((-1).^j.*besseli(1,lambda_j.*Rc)./(lambda_j.*besseli(0,lambda_j.*Rc)))'*Fq7_m0(2:end);
Fe(w,5) = 2*pi*a0*rho*tmp1;

Fe_FK(w,2)=(2*pi*rho*g*cosh(k0*(h-b))/(k0*cosh(k0*h)))*(Rbo*besselj(1,k0*Rbo) - Rc*besselj(1,k0*Rc));
Fe_FK(w,5)=(2*pi*rho*g*cosh(k0*(h-e2))/(k0*cosh(k0*h)))*(Rc*besselj(1,k0*Rc) );
	
%% m=1
B1 = 2*1i*B0;

Hq7(1:Ni,2)=B1.*(besselj(1,k0*Rbo)*((h-b)/h).*(M_ntau'*(Tp1_1n_Rb.*M_ntau(:,1))) - (k0*besselj(0,k0*Rbo)-besselj(1,k0*Rbo)/Rbo).*[1;zeros(Ni-1,1)]);
Hq7(Ni+1:end,2)=-B1*besselj(1,k0*Rbo).*(Tp1_1n_Rc.*M_ntau(:,1));

coeff = Dm1\Hq7(:,2);

Aq7_m1 = coeff(1:Ni);
Cq7_m1 = coeff(Ni+1:end);

Bq7_m1 = M_ntau*Aq7_m1 + B1*besselj(1,k0*Rbo).*M_ntau(:,1);
Fq7_m1 = L_jtau*Cq7_m1;

Fe(w,1)=-pi*a0*rho*Rbo*(B1.*besselj(1,k0*Rbo)*I1(1) + I1*Aq7_m1);
Fe(w,4)=- pi*a0*rho*Rc*(I2'*Cq7_m1);

tmp1 = Ti1_1n(1)*Bq7_m1(1) + sqrt(2)*sum((-1).^n.*Ti1_1n(2:end).*Bq7_m1(2:end));
tmp2 = Ti2_1n(1)*Cq7_m1(1) + sqrt(2)*sum((-1).^n.*Ti2_1n(2:end).*Cq7_m1(2:end));
Fe(w,3)=-pi*a0*rho*(Rbo*(B1.*besselj(1,k0*Rbo)*I3(1) + I3*Aq7_m1) + tmp1 + tmp2);

tmp1=.25*Fq7_m1(1)*Rc^3 + sqrt(2)*Rc^2*(((-1).^j.*besseli(2,lambda_j.*Rc))./(lambda_j.*besseli(1,lambda_j.*Rc)))'*Fq7_m1(2:end);   
Fe(w,6)=-pi*a0*rho*(Rc*I4'*Cq7_m1 + tmp1);
   
end% END OF FOR() LOOP


%% EXIT

% *** ATTENTION -> ici on reprend une dépendance en temps positive exp(-iwt) --> exp(iwt)
% Zr=A(w) i*B(w)/w;
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

save([DIR,filesep,'DATA/hydroParameters.mat'],'Omega','Fe','A','A_inf','B','Fe_Haskind','Fe_FK');
	
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
			Tp1(1,1) = -(m/a)*((r/a)^(-m-1)+(r/a)^(m-1))/(a/b-b/a);%((1/b + b/r^2)/(a/b-b/a));
			Tp2(1,1) = (m/b)*((r/b)^(-m-1)+(r/b)^(m-1))/(a/b-b/a);%((-a/r^2 - 1/a)/(a/b-b/a));
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
