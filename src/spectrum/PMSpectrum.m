function S_w = PMSpectrum(Hs, Omega)
% *** Pierson-Moskovtiz Spectrum characteristic

W=Omega; dw=W(2)-W(1); [m,n]=size(W);

g=9.81;

alpha=8.1e-3; beta=.74;

U=sqrt(.5*Hs*g*sqrt(beta/alpha));

S_w=alpha*g^2.*W.^(-5).*exp(-beta.*((g/U)./W).^4);

m0=sum(S_w)*dw;
m1=sum(Omega.*S_w)*dw;
m2=sum(Omega.^2.*S_w)*dw;

end

