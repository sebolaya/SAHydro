function display_hydroData(hydroParameters, nBody, DoF)

if nargin<4
	DoF=[1 2 3];
end

W=hydroParameters.Omega; nw=length(W);
A=hydroParameters.A;
B=hydroParameters.B;
Fe=hydroParameters.Fe;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%								SURGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure;
j=1;
% subplot(2,2,1), grid on, hold on;
plotyy(W,abs(Fe(:,j)),W,angle(Fe(:,j))*(180/pi));
% plot(Omega,abs(Fe_FK(:,j)),'b');
% plot(Omega,abs(Fe_Haskind(:,j)),'or');
% subplot(2,2,3), grid on, hold on;
% plot(Omega,angle(Fe_Haskind(:,j)),'or');
% plot(Omega,angle(Fe(:,j)),'r');
% plot(Omega,angle(Fe_FK(:,j)),'b');

% j=5;
% subplot(2,2,2), grid on, hold on;
% plot(Omega,abs(Fe_Haskind(:,j)),'or');
% plot(Omega,abs(Fe(:,j)),'r');
% plot(Omega,abs(Fe_FK(:,j)),'b');
% 
% subplot(2,2,4), grid on, hold on;
% plot(Omega,angle(Fe_Haskind(:,j)),'or')
% plot(Omega,angle(Fe(:,j)),'r');
% plot(Omega,angle(Fe_FK(:,j)),'b');
% 
% 
% %%% SURGE
% figure;
% i=1;j=1;
% subplot(2,2,1), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,2), grid on, hold on;
% i=4;j=1;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,3), grid on, hold on;
% i=1;j=4;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% i=4;j=4;
% subplot(2,2,4), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% %%% SURGE --> PITCH
% figure;
% i=1;j=3;
% subplot(2,2,1), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,2), grid on, hold on;
% i=4;j=3;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,3), grid on, hold on;
% i=1;j=6;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% i=4;j=6;
% subplot(2,2,4), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');

% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%								HEAVE
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if 1
%

j=2;
figure, grid on, hold on;
plotyy(W,abs(Fe(:,j)),W,angle(Fe(:,j))*(180/pi));

j=5;
figure, grid on, hold on;
plotyy(W,abs(Fe(:,j)),W,angle(Fe(:,j))*(180/pi));

% figure;
% i=2;j=2;
% subplot(2,2,1), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,2), grid on, hold on;
% i=5;j=2;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,3), grid on, hold on;
% i=2;j=5;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% i=5;j=5;
% subplot(2,2,4), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% 
% figure;
% j=2;
% subplot(2,2,1), grid on, hold on;
% plot(Omega,abs(Fe(:,j)),'r');
% plot(Omega,abs(Fe_FK(:,j)),'b');
% plot(Omega,abs(Fe_Haskind(:,j)),'or');
% 
% subplot(2,2,3), grid on, hold on;
% plot(Omega,angle(Fe_Haskind(:,j)),'or');
% plot(Omega,angle(Fe(:,j)),'r');
% plot(Omega,angle(Fe_FK(:,j)),'b');
% 
% j=5;
% subplot(2,2,2), grid on, hold on;
% plot(Omega,abs(Fe_Haskind(:,j)),'or');
% plot(Omega,abs(Fe(:,j)),'r');
% plot(Omega,abs(Fe_FK(:,j)),'b');
% 
% subplot(2,2,4), grid on, hold on;
% plot(Omega,angle(Fe_Haskind(:,j)),'or')
% plot(Omega,angle(Fe(:,j)),'r');
% plot(Omega,angle(Fe_FK(:,j)),'b');
% 
% end
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%								PITCH
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% if 1
% 
% figure;
% j=3;
% subplot(2,2,1), grid on, hold on;
% plot(Omega,abs(Fe(:,j)),'r');
% plot(Omega,abs(Fe_FK(:,j)),'b');
% plot(Omega,abs(Fe_Haskind(:,j)),'or');
% 
% subplot(2,2,3), grid on, hold on;
% plot(Omega,angle(Fe_Haskind(:,j)),'or');
% plot(Omega,angle(Fe(:,j)),'r');
% plot(Omega,angle(Fe_FK(:,j)),'b');
% 
% j=6;
% subplot(2,2,2), grid on, hold on;
% plot(Omega,abs(Fe_Haskind(:,j)),'or');
% plot(Omega,abs(Fe(:,j)),'r');
% plot(Omega,abs(Fe_FK(:,j)),'b');
% 
% subplot(2,2,4), grid on, hold on;
% plot(Omega,angle(Fe_Haskind(:,j)),'or')
% plot(Omega,angle(Fe(:,j)),'r');
% plot(Omega,angle(Fe_FK(:,j)),'b');
% 
% figure;
% i=3;j=3;
% subplot(2,2,1), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,2), grid on, hold on;
% i=6;j=3;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,3), grid on, hold on;
% i=3;j=6;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% i=6;j=6;
% subplot(2,2,4), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% 
% 
% %%% PITCH --> SURGE
% figure;
% i=3;j=1;
% subplot(2,2,1), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,2), grid on, hold on;
% i=3;j=4;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% subplot(2,2,3), grid on, hold on;
% i=6;j=1;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% i=6;j=4;
% subplot(2,2,4), grid on, hold on;
% plot(Omega,squeeze(A(i,j,:)),'-r');
% plot(Omega,A_inf(i,j)*ones(length(Omega),1),'-.r');
% plot(Omega,squeeze(B(i,j,:)),'-b');
% 
% end


end

