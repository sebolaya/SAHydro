function [Wn,Lambda_i] = natFreq(Omega, M, A, K, DoF, B)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Il nous faut un système sous la forme [K+iw*B(w)-(M+A(w))*w^2]*X=Fe
w=Omega;
nw=length(w);
nDoF=length(DoF);
w_0=.5;
Wn=zeros(nDoF,1);

% Calcul des fréquences propres
for j=1:nw
	Lambda(j,:)=eig(K,M+A(:,:,j));
end
	
for i=1:nDoF%length(DoF)

	[w_n,F_n,flag] = fsolve(@(x) eigs(x,Lambda(:,DoF(i)),w),w_0,optimset('Display','off','TolFun',1e-10));
	
% 	idx=find(Omega>=w_n,1);
% 	[V,D]=eig(K,M+A(:,:,idx),'qz');
% 	Phi(:,i)=V(:,i);
% 	
	
	%k=Phi(:,i)'*K*Phi(:,i);
% 	for j=1:nw
% 		tmp=Phi(:,i)'*(M+A(:,:,j))*Phi(:,i);
% 		bc(j,1)=2*tmp*w_n;
% 	end
% % 	tmp=Phi(:,i)'*(M+A(:,:,idx))*Phi(:,i);
% 	2*tmp*w_n
	
	if flag~=1
		disp('something wrong...');
	else
		Wn(i,1)=w_n;
	end
	
end

Lambda_i=Lambda(:,DoF);


function F = eigs(x,lam,w)
F = x.^2 - interp1(w,lam,x);

