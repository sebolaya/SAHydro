function [ out ] = f2( alpha, a, b, c, d )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    % % Int_(a)^(b){(z-c)*cos(alpha*(z+d))}
    
    if alpha == 0
        out = .5*((b-c)^2 - (a-c)^2);
    else
%         i1 = ( (alpha).*( b.*sin(alpha.*(d+b)) - a.*sin(alpha.*(d+a)) )...
%                 + ( cos(alpha.*(d+b)) - cos(alpha.*(d+a)) ) ) ./ ((alpha).^2);
%         i2 = c.*( sin(alpha.*(d+b)) - sin(alpha.*(d+a)) )./(alpha);
%         out = i1 - i2;
    out = ( (alpha).*((b-c).*sin(alpha.*(d+b)) - (a-c).*sin(alpha.*(d+a)))...
        +  ( cos(alpha.*(d+b)) - cos(alpha.*(d+a)) )) ./ alpha.^2;
    end

end
