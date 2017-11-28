function S_w = JONSWAPSpectrum(Hs, Tp, Omega, gamma)

	if nargin<4
		gamma=3.3;
	end
	
	W=Omega; dw=W(2)-W(1); [m,n]=size(W);
	
    w0=2*pi/Tp;%(3/5)^.25*(2*pi/Tp); 
	sigma=zeros(m,n); sigma(W<w0)=.07; sigma(W>=w0)=.09;
    
	a=exp(-(W-w0).^2./(2*sigma.^2*w0^2));

    S_w = 1.25*(Hs^2/4)*w0^4.*(W.^-5).*exp(-1.25.*(w0./W).^4).*gamma.^a;

    alpha = Hs^2/(16*sum(S_w)*dw);
    S_w = alpha*S_w;

    m0=sum(S_w)*dw;
    m1=sum(W.*S_w)*dw;
    m2=sum(W.^2.*S_w)*dw;

end

