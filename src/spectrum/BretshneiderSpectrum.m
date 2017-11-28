function S_w = BretshneiderSpectrum(Hs, Tp, Omega)

W=Omega; dw=W(2)-W(1); [m,n]=size(W);
	
w0=2*pi/Tp;%(3/5)^.25*(2*pi/Tp);

S_w=1.25*(Hs^2/4)*w0^4.*W.^(-5).*exp(-1.25*(w0./W).^4);

m0=sum(S_w)*dw
m1=sum(W.*S_w)*dw;
m2=sum(W.^2.*S_w)*dw;

end

