function [ out ] = f1( alpha, a, b, c, d )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
a=-a; b=-b;
    % % Int_(-a)^(-b){(z+c)^2*cos(alpha*(z+d))}

    if alpha == 0
        out = ((a+c)^3 - (b+c)^3) / 3;
    else
        out = ( (alpha).^2.*( (c-b)^2.*sin(alpha.*(d-b)) - (c-a)^2.*sin(alpha.*(d-a)) )...
            + 2.*(alpha).*( (c-b).*cos(alpha.*(d-b)) - (c-a).*cos(alpha.*(d-a)) )...
            - 2.*( sin(alpha.*(d-b)) - sin(alpha.*(d-a)) ) ) ./ ((alpha).^3);
    end

end
