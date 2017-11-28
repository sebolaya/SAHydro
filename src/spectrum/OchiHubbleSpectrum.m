function S_w = OchiHubbleSpectrum( Hs, Tp, Omega, lambda )

np=length(Hs);
nw=length(Omega);

S_w=zeros(nw,1);
for j=1:np

    omega_p = 2*pi/Tp(j);

    tmp1 = (((lambda(j)+.25)*omega_p^4)^lambda(j))/gamma(lambda(j));
    tmp2 = (Hs(j)^2./(Omega.^(4*lambda(j)+1))).*exp(-(lambda(j)+.25)*(omega_p./Omega).^4);
    
    S_w = S_w + .25.*tmp1.*tmp2;
end
	
end

