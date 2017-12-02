function [S_int,S_tild_int] = Sint_func(m,a,b,kappa)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

kappa_0 = imag(kappa(1));
kappa_l = kappa(2:end);

switch m
    case 0
        
        cst1 = besselh(0,kappa_0.*b).*besselh(0,2,kappa_0.*a) - besselh(0,kappa_0.*a).*besselh(0,2,kappa_0.*b);
        cst2 = besseli(0,kappa_l.*b).*besselk(0,kappa_l.*a) - besseli(0,kappa_l.*a).*besselk(0,kappa_l.*b);
        
        s_int_0 = (1./ (kappa_0.*cst1)) .*...
            ( (b.*besselh(1,kappa_0.*b) - a.*besselh(1,kappa_0.*a)).*besselh(0,2,kappa_0.*a)...
              - besselh(0,kappa_0.*a).*(b.*besselh(1,2,kappa_0.*b) - a.*besselh(1,2,kappa_0.*a)) );

        s_int_l = (1./ (kappa_l.*cst2)) .*...
            ( (b.*besseli(1,kappa_l.*b) - a.*besseli(1,kappa_l.*a)).*besselk(0,kappa_l.*a)...
              + besseli(0,kappa_l.*a).*(b.*besselk(1,kappa_l.*b) - a.*besselk(1,kappa_l.*a)) );     

        s_tild_int_0 = (1./ (kappa_0.*cst1)) .*...
            ( besselh(0,kappa_0.*b).*(b.*besselh(1,2,kappa_0.*b) - a.*besselh(1,2,kappa_0.*a))...
                - (b.*besselh(1,kappa_0.*b) - a.*besselh(1,kappa_0.*a)).*besselh(0,2,kappa_0.*b) );

        s_tild_int_l = (1./ (kappa_l.*cst2)) .*...
            ( - besseli(0,kappa_l.*b).*(b.*besselk(1,kappa_l.*b) - a.*besselk(1,kappa_l.*a))...
                -(b.*besseli(1,kappa_l.*b) - a.*besseli(1,kappa_l.*a)).* besselk(0,kappa_l.*b) );    

    case 1
        
        cst1 = besselh(1,kappa_0.*b).*besselh(1,2,kappa_0.*a) - besselh(1,kappa_0.*a).*besselh(1,2,kappa_0.*b);
        cst2 = besseli(1,kappa_l.*b).*besselk(1,kappa_l.*a) - besseli(1,kappa_l.*a).*besselk(1,kappa_l.*b);
 
        s_int_0 = (1./ (kappa_0.*cst1)) .*...
            ( (b^2.*besselh(2,kappa_0.*b) - a^2.*besselh(2,kappa_0.*a)).*besselh(1,2,kappa_0.*a)...
              - besselh(1,kappa_0.*a).*(b^2.*besselh(2,2,kappa_0.*b) - a^2.*besselh(2,2,kappa_0.*a)) );
        s_int_l = (1./ (kappa_l.*cst2)) .*...
            ( (b^2.*besseli(2,kappa_l.*b) - a^2.*besseli(2,kappa_l.*a)).*besselk(1,kappa_l.*a)...
              - besseli(1,kappa_l.*a).*(-b^2.*besselk(2,kappa_l.*b) + a^2.*besselk(2,kappa_l.*a)));
          
        s_tild_int_0 = (1./ (kappa_0.*cst1)) .*...
            ( besselh(1,kappa_0.*b).*(b^2.*besselh(2,2,kappa_0.*b) - a^2.*besselh(2,2,kappa_0.*a))...
             - (b^2.*besselh(2,kappa_0.*b) - a^2.*besselh(2,kappa_0.*a)).*besselh(1,2,kappa_0.*b) );

        s_tild_int_l = (1./ (kappa_l.*cst2)) .*...
            ( besseli(1,kappa_l.*b).*(-b^2.*besselk(2,kappa_l.*b) + a^2.*besselk(2,kappa_l.*a))...
                - (b^2.*besseli(2,kappa_l.*b) - a^2.*besseli(2,kappa_l.*a)).*besselk(1,kappa_l.*b) );
    
end    
    S_int = [s_int_0;s_int_l];
    S_tild_int = [s_tild_int_0;s_tild_int_l];
    
end