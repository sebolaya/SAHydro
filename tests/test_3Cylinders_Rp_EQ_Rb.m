clear; close all; clc;

%%% ***********************************************************************
%%%                  Définition des paramètres de simulation
%%% ***********************************************************************

% Extension File
extension = '.dat';

% spécifie le répertoire de travail
DIR = '~/Documents/MATLAB/Project/THESE/HydroFunctions/3Cylinders_Rp_EQ_Rb/';
DIR_NEMOH = '~/Documents/MATLAB/Project/NEMOH/2Cylinders/';

REP0 = 'DATA/';


options.DIR = DIR;

h = 150;
rho = 1025;
g = 9.81;

WECstructure.Rbo = 15;          % diamètre exterieur - dDe
WECstructure.Rc = 5.5;%8.3         % diamètre de la colonne - dCeh  --> 8
WECstructure.b = 10;           % tirant d'eau flotteur
WECstructure.e1 = 25.5;          % tirant d'eau embase-colonne
WECstructure.e2 = WECstructure.e1+10;        % hauteur de l'embase heB

WECstructure.Rp=5.50001;%2.501;%

% Numerical truncation 
options.Truncate.Ni = 55;
options.Truncate.Nl = 150;
options.Truncate.Nn = 150;
options.Truncate.Nj = 150;

options.Zc = [0,0]';

options.Haskind = 0;
options.matchingConditions = 0; % *** ATTENTION penser à définir la pulsation omega dans le fichier correspondant au cas d'étude

Omega=(.1:.1:3)'; nw=length(Omega);

[Fe, A, B, A_inf, Fe_Haskind] = threeCylinders_Rp_INF_Rbo__Rbi_EQ_Rc(Omega, h, WECstructure, options);



if 1
% options.Haskind=0; options.dof='all';
% A2=zeros(6,6,nw); B2=zeros(6,6,nw);
% w = 1;
% while w <=nw
% 	omega = Omega(w);
% 	Lambda_i=fcn_waveDispersion(Omega(w),h,81);
% 	out = twoCylinder_Rp_inf_a_Rb_withColumn(WECstructure,omega,Lambda_i,options);
% 
% 	Fe2(w,1)=out.fx_1;Fe2(w,4)=out.fx_2;
% 	Fe2_Haskind(w,1)=out.fx_1_Haskind;Fe2_Haskind(w,4)=out.fx_2_Haskind;
% 
% 	Fe2(w,2)=out.fz_1;Fe2(w,5)=out.fz_2;
% 	Fe2_Haskind(w,2)=out.fz_1_Haskind;Fe2_Haskind(w,5)=out.fz_2_Haskind;
% 
% 	Fe2(w,3)=out.ty_1;Fe2(w,6)=out.ty_2;
% 	Fe2_Haskind(w,3)=out.ty_1_Haskind;Fe2_Haskind(w,6)=out.ty_2_Haskind;
% 
% 	% *** F_r(w) = - Z_r(w) * u(w) avec Z_r = B(w) + iw*M_a(w)
% 	A2(1,1,w) = imag(out.z_r11_11)./omega; B2(1,1,w) = real(out.z_r11_11);
% 	A2(1,4,w) = imag(out.z_r11_21)./omega; B2(1,4,w) = real(out.z_r11_21);
% 	A2(4,1,w) = imag(out.z_r11_12)./omega; B2(4,1,w) = real(out.z_r11_12);
% 	A2(4,4,w) = imag(out.z_r11_22)./omega; B2(4,4,w) = real(out.z_r11_22);
% 	
% 	A2(2,2,w) = imag(out.z_r33_11)./omega; B2(2,2,w) = real(out.z_r33_11);
% 	A2(2,5,w) = imag(out.z_r33_21)./omega; B2(2,5,w) = real(out.z_r33_21);
% 	A2(5,2,w) = imag(out.z_r33_12)./omega; B2(5,2,w) = real(out.z_r33_12);
% 	A2(5,5,w) = imag(out.z_r33_22)./omega; B2(5,5,w) = real(out.z_r33_22);
% 
% 	A2(3,3,w) = imag(out.z_r55_11)./omega; B2(3,3,w) = real(out.z_r55_11);
% 	A2(3,6,w) = imag(out.z_r55_21)./omega; B2(3,6,w) = real(out.z_r55_21);
% 	A2(6,3,w) = imag(out.z_r55_12)./omega; B2(6,3,w) = real(out.z_r55_12);
% 	A2(6,6,w) = imag(out.z_r55_22)./omega; B2(6,6,w) = real(out.z_r55_22);
% 
% 	A2(3,1,w) = imag(out.z_r15_11)./omega; B2(3,1,w) = real(out.z_r15_11);
% 	A2(3,4,w) = imag(out.z_r15_21)./omega; B2(3,4,w) = real(out.z_r15_21);
% 	A2(6,1,w) = imag(out.z_r15_12)./omega; B2(6,1,w) = real(out.z_r15_12);
% 	A2(6,4,w) = imag(out.z_r15_22)./omega; B2(6,4,w) = real(out.z_r15_22);
% 
% 	A2(1,3,w) = imag(out.z_r51_11)./omega; B2(1,3,w) = real(out.z_r51_11);
% 	A2(1,6,w) = imag(out.z_r51_21)./omega; B2(1,6,w) = real(out.z_r51_21);
% 	A2(4,3,w) = imag(out.z_r51_12)./omega; B2(4,3,w) = real(out.z_r51_12);
% 	A2(4,6,w) = imag(out.z_r51_22)./omega; B2(4,6,w) = real(out.z_r51_22);
% 	
% 	w = w +1;
% end



% [Fe2, A2, B2, A_inf2, Fe_Haskind2]=threeCylinders_Rp_EQ_Rbo__Rbi_EQ_Rc(Omega, h, WECstructure, options);
[Fe2, A2, B2, A_inf2, Fe_Haskind2]=twoCylinders_COLUMN(Omega, h, WECstructure, options);

% WECstructure.Rp = 15.0001;
% [Fe3, A3, B3, A_inf3, Fe_Haskind3]=threeCylinders_Rp_SUP_Rbo__Rbi_EQ_Rc(Omega, h, WECstructure, options);

% i=3; j=6;
% i=1; j=4;
i=2; j=5;

figure; hold on;
plot(Omega,squeeze(A(i,j,:)))
plot(Omega,squeeze(A(j,i,:)),'+')
plot(Omega,squeeze(A2(i,j,:)),'-.r')
plot(Omega,squeeze(A2(j,i,:)),'or')
% plot(Omega,squeeze(A3(i,j,:)),'-.k')
% plot(Omega,squeeze(A3(j,i,:)),'*k')

figure; hold on;
plot(Omega,squeeze(A(i,i,:)))
plot(Omega,squeeze(A2(i,i,:)),'+')
% plot(Omega,squeeze(A3(i,i,:)),'o')

figure; hold on;
plot(Omega,squeeze(A(j,j,:)))
plot(Omega,squeeze(A2(j,j,:)),'+')
% plot(Omega,squeeze(A3(j,j,:)),'o')
end

if 0

i=3; j=i; k=i;
	
figure, grid on, hold on;
[f1.hAx,f1.hL1,f1.hL2]=plotyy(Omega,[squeeze(A(i,j,:)) A_inf(i,j)*ones(length(Omega),1)],Omega,squeeze(B(i,j,:)));
hold(f1.hAx(1),'on'); hold(f1.hAx(2),'on');
set(f1.hL1(2),'Color','b','LineStyle','-.');
% ylim=get(f1.hAx(1),'Ylim'); set(f1.hAx(1),'Ylim',[0 ylim(2)]);
% ylim=get(f1.hAx(2),'Ylim'); set(f1.hAx(2),'Ylim',[0 ylim(2)]);

figure, grid on, hold on;
[f2.hAx,f2.hL1,f2.hL2]=plotyy(Omega,[abs(Fe(:,k)) abs(Fe_Haskind(:,k))],Omega,[angle(Fe(:,k)) angle(Fe_Haskind(:,k))]);
hold(f2.hAx(1),'on'); hold(f2.hAx(2),'on');
set(f2.hL1(2),'Marker','+','LineStyle','none'); set(f2.hL2(2),'Marker','+','LineStyle','none');


% [Fe2, A2, B2, A_inf2, Fe_Haskind2]=twoCylinders_COLUMN(Omega, h, WECstructure, options);
[Fe2, A2, B2, A_inf2, Fe_Haskind2]=threeCylinders_Rp_EQ_Rbo__Rbi_EQ_Rc(Omega, h, WECstructure, options);

plot(f1.hAx(1),Omega,[squeeze(A2(i,j,:))],'+');
plot(f1.hAx(2),Omega,squeeze(B2(i,j,:)),'+');

plot(f2.hAx(1),Omega,abs(Fe2(:,k)),'o');

WECstructure.Rp = 15.0001;
[Fe3, A3, B3, A_inf3, Fe_Haskind3]=threeCylinders_Rp_SUP_Rbo__Rbi_EQ_Rc(Omega, h, WECstructure, options);

plot(f1.hAx(1),Omega,[squeeze(A3(i,j,:))],'o');
plot(f1.hAx(2),Omega,squeeze(B3(i,j,:)),'o');

plot(f2.hAx(1),Omega,abs(Fe3(:,k)),'*');
end
