clear; close all; clc;

%%% ***********************************************************************
%%%                  Définition des paramètres de simulation
%%% ***********************************************************************

% Extension File
extension = '.dat';

% spécifie le répertoire de travail
DIR = '~/Documents/MATLAB/Project/THESE/HydroFunctions/3Cylinders_Rp_SUP_Rb/';
DIR_NEMOH = '~/Documents/MATLAB/Project/NEMOH/2Cylinders/';

REP0 = 'DATA/';


options.DIR = [DIR REP0];
if ~isdir(options.DIR); mkdir(options.DIR); end;


h = 150;
rho = 1025;
g = 9.81;

WECstructure.Rbo = .15;          % diamètre exterieur - dDe
WECstructure.Rc = .05;%8.3         % diamètre de la colonne - dCeh  --> 8
WECstructure.Rp = 5.0;          % diamètre de l'embase - deB
WECstructure.b = .05;           % tirant d'eau flotteur
WECstructure.e1 = 10;          % tirant d'eau embase-colonne
WECstructure.e2 = WECstructure.e1+5;        % hauteur de l'embase heB


options.Truncate.Ni = 80;%max 300;
options.Truncate.Nl = 80;%max 300
options.Truncate.Nn = 150;%100
options.Truncate.Nj = 150;%100

options.Zc = [0;0];%[0;-22.9];%

options.matchingConditions = 0;
options.matchingFrequency = 1;%rad/s


Omega = (.1:.1:3)';
options.dof = 'all';


[Fe, A, B, A_inf Fe_Haskind]=threeCylinders_Rp_SUP_Rbo__Rbi_EQ_Rc(Omega, h, WECstructure, options);

i=5; j=5; k=5;
figure, grid on, hold on;
[hAx,hL1,hL2]=plotyy(Omega,[squeeze(A(i,j,:)) A_inf(i,j)*ones(length(Omega),1)],Omega,squeeze(B(i,j,:)));
set(hL1(2),'Color','b','LineStyle','-.');

figure, grid on, hold on;
[hAx,hL1,hL2]=plotyy(Omega,[abs(Fe(:,k)) abs(Fe_Haskind(:,k))],Omega,[angle(Fe(:,k)) angle(Fe_Haskind(:,k))]);
set(hL1(2),'Marker','+','LineStyle','none'); set(hL2(2),'Marker','+','LineStyle','none');


if 0

% load([options.DIR 'hydroParameters.mat'],'Omega','Ma11','B11','Fx','Ma55','B55','Ty');
data_NEMOH = load([DIR_NEMOH 'hydroParameters_NEMOH_zG_mesh3.mat'],'Omega','A','B','Fe');
data_NEMOH.Fe = conj(data_NEMOH.Fe);

% 
figure;
subplot(2,2,1), grid on, hold on;
plot(Omega,squeeze(data1.Ma11(1,1,:)),'r');
plot(Omega,squeeze(data1.B11(1,1,:)));
plot(Omega,squeeze(A(1,1,:)),'*r');
plot(Omega,squeeze(B(1,1,:)),'*');
plot(Omega,squeeze(data_NEMOH.A(1,1,:)),'-g');
plot(Omega,squeeze(data_NEMOH.B(1,1,:)),'-g');
%%%
subplot(2,2,2), grid on, hold on;
plot(Omega,squeeze(data1.Ma11(1,2,:)),'r');
plot(Omega,squeeze(data1.B11(1,2,:)));
plot(Omega,squeeze(A(1,4,:)),'*r');
plot(Omega,squeeze(B(1,4,:)),'*');
%%%
subplot(2,2,3), grid on, hold on;
plot(Omega,squeeze(data1.Ma11(2,1,:)),'r');
plot(Omega,squeeze(data1.B11(2,1,:)));
plot(Omega,squeeze(A(4,1,:)),'*r');
plot(Omega,squeeze(B(4,1,:)),'*');
%%%
subplot(2,2,4), grid on, hold on;
plot(Omega,squeeze(data1.Ma11(2,2,:)),'r');
plot(Omega,squeeze(data1.B11(2,2,:)));
plot(Omega,squeeze(A(4,4,:)),'*r');
plot(Omega,squeeze(B(4,4,:)),'*');
plot(Omega,squeeze(data_NEMOH.A(4,4,:)),'-g');
plot(Omega,squeeze(data_NEMOH.B(4,4,:)),'-g');
% 
%

figure;
subplot(2,2,1), grid on, hold on;
plot(Omega,squeeze(data1.Ma33(1,1,:)),'r');
plot(Omega,squeeze(data1.B33(1,1,:)));
plot(Omega,squeeze(A(2,2,:)),'*r');
plot(Omega,squeeze(B(2,2,:)),'*');
plot(Omega,squeeze(data_NEMOH.A(2,2,:)),'-g');
plot(Omega,squeeze(data_NEMOH.B(2,2,:)),'-g');
%%%
subplot(2,2,2), grid on, hold on;
plot(Omega,squeeze(data1.Ma33(1,2,:)),'r');
plot(Omega,squeeze(data1.B33(1,2,:)));
plot(Omega,squeeze(A(5,2,:)),'*r');
plot(Omega,squeeze(B(5,2,:)),'*');
%%%
subplot(2,2,3), grid on, hold on;
plot(Omega,squeeze(data1.Ma33(2,1,:)),'r');
plot(Omega,squeeze(data1.B33(2,1,:)));
plot(Omega,squeeze(A(2,5,:)),'*r');
plot(Omega,squeeze(B(2,5,:)),'*');
%%%
subplot(2,2,4), grid on, hold on;
plot(Omega,squeeze(data1.Ma33(2,2,:)),'r');
plot(Omega,squeeze(data1.B33(2,2,:)));
plot(Omega,squeeze(A(5,5,:)),'*r');
plot(Omega,squeeze(B(5,5,:)),'*');
plot(Omega,squeeze(data_NEMOH.A(5,5,:)),'-g');
plot(Omega,squeeze(data_NEMOH.B(5,5,:)),'-g');

figure;
subplot(2,2,1), grid on, hold on;
plot(Omega,squeeze(data1.Ma55(1,1,:)),'r');
plot(Omega,squeeze(data1.B55(1,1,:)));
plot(Omega,squeeze(A(3,3,:)),'*r');
plot(Omega,squeeze(B(3,3,:)),'*');
plot(Omega,squeeze(data_NEMOH.A(3,3,:)),'-g');
plot(Omega,squeeze(data_NEMOH.B(3,3,:)),'-g');
%%%
subplot(2,2,2), grid on, hold on;
plot(Omega,squeeze(data1.Ma55(1,2,:)),'r');
plot(Omega,squeeze(data1.B55(1,2,:)));
plot(Omega,squeeze(A(6,3,:)),'*r');
plot(Omega,squeeze(B(6,3,:)),'*');
%%%
subplot(2,2,3), grid on, hold on;
plot(Omega,squeeze(data1.Ma55(2,1,:)),'r');
plot(Omega,squeeze(data1.B55(2,1,:)));
plot(Omega,squeeze(A(3,6,:)),'*r');
plot(Omega,squeeze(B(3,6,:)),'*');
%%%
subplot(2,2,4), grid on, hold on;
plot(Omega,squeeze(data1.Ma55(2,2,:)),'r');
plot(Omega,squeeze(data1.B55(2,2,:)));
plot(Omega,squeeze(A(6,6,:)),'*r');
plot(Omega,squeeze(B(6,6,:)),'*');
plot(Omega,squeeze(data_NEMOH.A(6,6,:)),'-g');
plot(Omega,squeeze(data_NEMOH.B(6,6,:)),'-g');
% 
% 



figure;
subplot(2,2,1), grid on, hold on;
plot(Omega,squeeze(abs(data1.Fx(1,1,:))));
plot(Omega,abs(Fe_Haskind(:,1)),'*r');
plot(Omega,abs(data_NEMOH.Fe(:,1)),'-g');
subplot(2,2,3), grid on, hold on;
plot(Omega,squeeze(angle(data1.Fx(1,1,:))));
plot(Omega,angle(Fe_Haskind(:,1)),'*r');
plot(Omega,angle(data_NEMOH.Fe(:,1)),'-g');
subplot(2,2,2), grid on, hold on;
plot(Omega,squeeze(abs(data1.Fx(2,1,:))));
plot(Omega,abs(Fe_Haskind(:,4)),'*r');
plot(Omega,abs(data_NEMOH.Fe(:,4)),'-g');
subplot(2,2,4), grid on, hold on;
plot(Omega,squeeze(angle(data1.Fx(2,1,:))));
plot(Omega,angle(Fe_Haskind(:,4)),'*r');
plot(Omega,angle(data_NEMOH.Fe(:,1)),'-g');


figure;
subplot(2,2,1), grid on, hold on;
plot(Omega,squeeze(abs(data1.Fz(1,1,:))));
plot(Omega,abs(Fe(:,2)),'r');
plot(Omega,abs(Fe_Haskind(:,2)),'*r');
plot(Omega,abs(data_NEMOH.Fe(:,2)),'-g');
subplot(2,2,3), grid on, hold on;
plot(Omega,squeeze(angle(data1.Fz(1,1,:))));
plot(Omega,angle(Fe(:,2)),'r');
plot(Omega,angle(Fe_Haskind(:,2)),'*r');
plot(Omega,angle(data_NEMOH.Fe(:,2)),'-g');
subplot(2,2,2), grid on, hold on;
plot(Omega,squeeze(abs(data1.Fz(2,1,:))));
plot(Omega,abs(Fe(:,5)),'r');
plot(Omega,abs(Fe_Haskind(:,5)),'*r');
plot(Omega,abs(data_NEMOH.Fe(:,5)),'-g');
subplot(2,2,4), grid on, hold on;
plot(Omega,squeeze(angle(data1.Fz(2,1,:))));
plot(Omega,angle(Fe(:,5)),'r');
plot(Omega,angle(Fe_Haskind(:,5)),'*r');
plot(Omega,angle(data_NEMOH.Fe(:,5)),'-g');

figure;
subplot(2,2,1), grid on, hold on;
plot(Omega,squeeze(abs(data1.Ty(1,1,:))));
plot(Omega,abs(Fe(:,3)),'r');
plot(Omega,abs(Fe_Haskind(:,3)),'*r');
plot(Omega,abs(data_NEMOH.Fe(:,3)),'-g');
subplot(2,2,3), grid on, hold on;
plot(Omega,squeeze(angle(data1.Ty(1,1,:))));
plot(Omega,angle(Fe(:,3)),'r');
plot(Omega,angle(Fe_Haskind(:,3)),'*r');
plot(Omega,angle(data_NEMOH.Fe(:,3)),'-g');
subplot(2,2,2), grid on, hold on;
plot(Omega,squeeze(abs(data1.Ty(2,1,:))));
plot(Omega,abs(Fe(:,6)),'r');
plot(Omega,abs(Fe_Haskind(:,6)),'*r');
plot(Omega,abs(data_NEMOH.Fe(:,6)),'-g');
subplot(2,2,4), grid on, hold on;
plot(Omega,squeeze(angle(data1.Ty(2,1,:))));
plot(Omega,angle(Fe(:,6)),'r');
plot(Omega,angle(Fe_Haskind(:,6)),'*r');
plot(Omega,angle(data_NEMOH.Fe(:,6)),'-g');
end