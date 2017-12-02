function [S_prim,S_tild_prim] = Sprim_func(m,r,R1,R2,kappa)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

kappa_0 = imag(kappa(1));
kappa_l = kappa(2:end);

switch m
    case 0

        cst = besselh(0,kappa_0.*R1).*besselh(0,2,kappa_0.*R2) - besselh(0,kappa_0.*R2).*besselh(0,2,kappa_0.*R1);

        s_prim_0 = (r.*kappa_0./cst).*...
            ( - besselh(1,kappa_0.*r).*besselh(0,2,kappa_0.*R2) + besselh(0,kappa_0.*R2).*besselh(1,2,kappa_0.*r) );

        s_tild_prim_0 = (r.*kappa_0./cst).*...
            ( -besselh(0,kappa_0.*R1).*besselh(1,2,kappa_0.*r) + besselh(1,kappa_0.*r).*besselh(0,2,kappa_0.*R1) );

        cst2 = besseli(0,kappa_l.*R1).*besselk(0,kappa_l.*R2) - besseli(0,kappa_l.*R2).*besselk(0,kappa_l.*R1);

        s_prim_l = (r.*kappa_l./cst2).*...
            (besseli(1,kappa_l.*r).*besselk(0,kappa_l.*R2) + besseli(0,kappa_l.*R2).*besselk(1,kappa_l.*r));

        s_tild_prim_l = (-r.*kappa_l./cst2).*...
            (besseli(1,kappa_l.*r).*besselk(0,kappa_l.*R1) + besseli(0,kappa_l.*R1).*besselk(1,kappa_l.*r) );
    
    case 1
        
        cst = besselh(1,kappa_0.*R1).*besselh(1,2,kappa_0.*R2) - besselh(1,kappa_0.*R2).*besselh(1,2,kappa_0.*R1);

        s_prim_0 = (r.*kappa_0./cst).*...
            ( (besselh(0,kappa_0.*r) - besselh(1,kappa_0.*r)./(kappa_0.*r)).*besselh(1,2,kappa_0.*R2)...
                    - besselh(1,kappa_0.*R2).*(besselh(0,2,kappa_0.*r) - besselh(1,2,kappa_0.*r)./(kappa_0.*r)) );

        s_tild_prim_0 = (r.*kappa_0./cst).*...
            ( besselh(1,kappa_0.*R1).*(besselh(0,2,kappa_0.*r) - besselh(1,2,kappa_0.*r)./(kappa_0.*r))...
                    - (besselh(0,kappa_0.*r) - besselh(1,kappa_0.*r)./(kappa_0.*r)).*besselh(1,2,kappa_0.*R1) );

        cst2 = besseli(1,kappa_l.*R1).*besselk(1,kappa_l.*R2) - besseli(1,kappa_l.*R2).*besselk(1,kappa_l.*R1);

        s_prim_l = (r.*kappa_l./cst2).*...
            ( (besseli(0,kappa_l.*r) - besseli(1,kappa_l.*r)./(kappa_l.*r)).*besselk(1,kappa_l.*R2)...
                    - besseli(1,kappa_l.*R2).*(- besselk(0,kappa_l.*r) - besselk(1,kappa_l.*r)./(kappa_l.*r)) );

        s_tild_prim_l = (r.*kappa_l./cst2).*...
            ( besseli(1,kappa_l.*R1).*(- besselk(0,kappa_l.*r) - besselk(1,kappa_l.*r)./(kappa_l.*r))...
                    - (besseli(0,kappa_l.*r) - besseli(1,kappa_l.*r)./(kappa_l.*r)).*besselk(1,kappa_l.*R1)  );     
end    
    S_prim = [s_prim_0;s_prim_l];
    S_tild_prim = [s_tild_prim_0;s_tild_prim_l];
    
end

