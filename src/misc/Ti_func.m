function [T_int,T_tild_int] = Tint_func(m,a,b,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch m
    case 0
        
        cst = besseli(0,alpha.*b).*besselk(0,alpha.*a)...
            - besseli(0,alpha.*a).*besselk(0,alpha.*b);
        
        t_int_0 = (.5*b^2*log(b/a) - .25*(b^2-a^2))/log(b/a);
        t_tild_int_0 = (-.5*a^2*log(b/a) + .25*(b^2-a^2))/log(b/a);

        t_int_n = (1./ (alpha.*cst)) .*...
                ( besselk(0,alpha.*a).*(b.*besseli(1,alpha.*b) - a.*besseli(1,alpha.*a)) +...
                  besseli(0,alpha.*a).*(b.*besselk(1,alpha.*b) - a.*besselk(1,alpha.*a)));
        t_tild_int_n = (-1./ (alpha.*cst)) .*...
                ( besselk(0,alpha.*b).*(b.*besseli(1,alpha.*b) - a.*besseli(1,alpha.*a)) +...
                  besseli(0,alpha.*b).*(b.*besselk(1,alpha.*b) - a.*besselk(1,alpha.*a)));        

    case 1
        
        cst = besseli(1,alpha.*b).*besselk(1,alpha.*a) - besseli(1,alpha.*a).*besselk(1,alpha.*b);

        t_int_0 = (.25*(b^4-a^4)/a - .5*(b^2-a^2)*a) / (b/a - a/b);
        t_tild_int_0 = (.5*(b^2-a^2)*b - .25*(b^4-a^4)/b) / (b/a - a/b);
        
        t_int_n = (1./ (alpha.*cst)) .*...
                ( (b^2.*besseli(2,alpha.*b) - a^2.*besseli(2,alpha.*a)).*besselk(1,alpha.*a)...
                    - besseli(1,alpha.*a).*(-b^2.*besselk(2,alpha.*b) + a^2.*besselk(2,alpha.*a)) );
        t_tild_int_n = (1./ (alpha.*cst)) .*...
                ( besseli(1,alpha.*b).*(-b^2.*besselk(2,alpha.*b) + a^2.*besselk(2,alpha.*a))...
                    - (b^2.*besseli(2,alpha.*b) - a^2.*besseli(2,alpha.*a)).*besselk(1,alpha.*b) );
end

T_int = [t_int_0;t_int_n];
T_tild_int = [t_tild_int_0;t_tild_int_n];

end