function [T_prim,T_tild_prim] = Tprim_func(m,r,R1,R2,alpha)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

switch m
    case 0
        
        cst = besseli(0,alpha.*R1).*besselk(0,alpha.*R2)...
            - besseli(0,alpha.*R2).*besselk(0,alpha.*R1);

        t_prim_0 = 1/log(R1/R2);
        t_prim_n = ((r.*alpha)./cst).*...
            (besseli(1,alpha.*r).*besselk(0,alpha.*R2) + besseli(0,alpha.*R2).*besselk(1,alpha.*r));

        t_tild_prim_0 = -1/log(R1/R2);
        t_tild_prim_n = (-(r.*alpha)./cst).*...
            (besseli(1,alpha.*r).*besselk(0,alpha.*R1) + besseli(0,alpha.*R1).*besselk(1,alpha.*r) );
        
    case 1
        
        cst = besseli(1,alpha.*R1).*besselk(1,alpha.*R2)...
            - besseli(1,alpha.*R2).*besselk(1,alpha.*R1);
        
        t_prim_0 = r*((1/R2 + R2/r^2)/(R1/R2-R2/R1));
        t_prim_n = (r.*alpha./cst).*...
            ( (besseli(0,alpha.*r) - besseli(1,alpha.*r)./(alpha.*r)).*besselk(1,alpha.*R2)...
                    - besseli(1,alpha.*R2).*(- besselk(0,alpha.*r) - besselk(1,alpha.*r)./(alpha.*r)) );

        t_tild_prim_0 = r*((-R1/r^2 - 1/R1)/(R1/R2-R2/R1));
        t_tild_prim_n = (r.*alpha./cst).*...
            ( besseli(1,alpha.*R1).*(- besselk(0,alpha.*r) - besselk(1,alpha.*r)./(alpha.*r))...
                    - (besseli(0,alpha.*r) - besseli(1,alpha.*r)./(alpha.*r)).*besselk(1,alpha.*R1)  );
                
end
T_prim = [t_prim_0;t_prim_n];
T_tild_prim = [t_tild_prim_0;t_tild_prim_n];
end
