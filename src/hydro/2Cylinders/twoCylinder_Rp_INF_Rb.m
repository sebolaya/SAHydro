function out = twoCylinder_Rp_INF_Rb(WECstructure,omega,Lambda_i,options)

    g = 9.81;
    rho = 1025;    
    h = 150;
    
    Rbo = WECstructure.Rbo;
    Rp = WECstructure.Rp;
    
    b = WECstructure.b;
    e1 = WECstructure.e1;
    e2 = WECstructure.e2;
    
    Ni = options.Truncate.Ni;
    Nn = options.Truncate.Nn;
    Nj = options.Truncate.Nj;
    Nl = options.Truncate.Nl;
    
    epsilon = 1e-7;

    n = (1:1:Nn)';
    j = (1:1:Nj)';
    l = (1:1:Nl)';
    
    beta_j = (pi/(h-e2)).*j;
    alpha_n = (pi/(e1-b)).*n;   
    gamma_l = (pi/(h-b)).*l;

    lambda_0 = -imag(Lambda_i(1));
    lambda_i = Lambda_i(2:Ni+1);
    
    N_lambda_i = 0.5.*( 1 + sin(2.*Lambda_i.*h)./(2.*Lambda_i.*h));

    K_ltau = (integraleInterface([0;gamma_l],Lambda_i,h,b,h,h,epsilon)*diag(N_lambda_i.^(-.5)))./(h-b);
    K_ltau(2:end,:) = sqrt(2).*K_ltau(2:end,:); 

    L_jtau = (integraleInterface([0;beta_j],[0;gamma_l],h,e2,h,h,epsilon))./(h-e2);
    L_jtau(2:end,1) = sqrt(2).*L_jtau(2:end,1);
    L_jtau(1,2:end) = sqrt(2).*L_jtau(1,2:end);
    L_jtau(2:end,2:end) = 2.*L_jtau(2:end,2:end);
    
    M_ntau = (integraleInterface([0;alpha_n],[0;gamma_l],e1,b,e1,h,epsilon))./(e1-b);
    M_ntau(2:end,1) = sqrt(2).*M_ntau(2:end,1);
    M_ntau(1,2:end) = sqrt(2).*M_ntau(1,2:end);
    M_ntau(2:end,2:end) = 2.*M_ntau(2:end,2:end);

    gamma_00 = 0;
    gamma_0j = beta_j.*Rp.*(besseli(1,beta_j.*Rp) ./ besseli(0,beta_j.*Rp));
    
    gamma_10 = 1;
    gamma_1j = beta_j.*Rp.*(besseli(0,beta_j.*Rp) ./ besseli(1,beta_j.*Rp)) - 1;
    
    pi_00 = 0;
    pi_0n = alpha_n.*Rp.*(besseli(1,alpha_n.*Rp) ./ besseli(0,alpha_n.*Rp));
    
    pi_10 = 1;
    pi_1n = alpha_n.*Rp.*(besseli(0,alpha_n.*Rp) ./ besseli(1,alpha_n.*Rp)) - 1;
    
    delta_00 = lambda_0*Rbo*(besselh(1,lambda_0*Rbo) / besselh(0,lambda_0*Rbo));
    delta_0i = lambda_i.*Rbo.*(besselk(1,lambda_i.*Rbo) ./ besselk(0,lambda_i.*Rbo));
    
    delta_10 = 1 - lambda_0*Rbo*(besselh(0,lambda_0*Rbo) / besselh(1,lambda_0*Rbo));
    delta_1i = 1 + lambda_i.*Rbo.*( besselk(0,lambda_i.*Rbo) ./ besselk(1,lambda_i.*Rbo) );

    Gamma_0j = [gamma_00;gamma_0j];
    Gamma_1j = [gamma_10;gamma_1j];
    
    Delta_0i = [delta_00 delta_0i];
    Delta_1i = [delta_10 delta_1i];

    Pi_0n = [pi_00;pi_0n];
    Pi_1n = [pi_10;pi_1n];
    
    [T_0l_prim_Rb,T_0l_tild_prim_Rb] = Tp_func(0,Rbo,Rbo,Rp,gamma_l);
    [T_0l_prim_Rp,T_0l_tild_prim_Rp] = Tp_func(0,Rp,Rbo,Rp,gamma_l);
    
    [T_1l_prim_Rb,T_1l_tild_prim_Rb] = Tp_func(1,Rbo,Rbo,Rp,gamma_l);
    [T_1l_prim_Rp,T_1l_tild_prim_Rp] = Tp_func(1,Rp,Rbo,Rp,gamma_l);
    
    [T_0l_int,T_0l_tild_int] = Ti_func(0,Rp,Rbo,gamma_l);
    [T_1l_int,T_1l_tild_int] = Ti_func(1,Rp,Rbo,gamma_l);
    
    
    I1 = N_lambda_i.^(-.5).*( sin(Lambda_i.*h)-sin(Lambda_i.*(h-b)) ) ./ Lambda_i;              % % Int_(-b)^(0){ Z_i^I }
    I2 = [(e2-e1);sqrt(2).*(sin(gamma_l.*(h-e1)) - sin(gamma_l.*(h-e2)))./gamma_l];             % % Int_(-e2)^(-e1){ Z_l^II }
    I3 = N_lambda_i.^(-.5).*f2(Lambda_i,-b,0,0,h);                                                 % % Int_(-b)^(0){ z*Z_i^I }
    I4 = [.5*(e1^2-e2^2);sqrt(2).*f2(gamma_l,-e2,-e1,0,h)];                                         % % Int_(-e2)^(-e1){ z*Z_l^II }
    I5_1 = [(e1-b)^3/3;sqrt(2).*f1(gamma_l,e1,b,e1,h)];                                         % % Int_{-e1}^{-b}{(z+e1)^2*Z_l^II}
    I5_2 = [(e1-b)^3/3;sqrt(2).*f1(gamma_l,e1,b,b,h)];                                          % % Int_{-e1}^{-b}{(z+b)^2*Z_l^II}
    I6 = [(h-b)^3/3;sqrt(2).*f1(gamma_l,h,b,h,h)];                                              % % Int_{-h}^{-b}{(z+h)^2*Z_l^II}
    I7 = N_lambda_i.^(-.5).*f1(Lambda_i,h,b,h,h);                                               % % Int_(-h)^(-b){ (z+h)^2*Z_i^I }
    I8 = [(h-e2)^3/3;sqrt(2).*f1(gamma_l,h,e2,h,h)];                                            % % Int_{-h}^{-e2}{(z+h)^2*Z_l^II}

%%% ***********************************************************************
%%% ***         Définition des matrices intrinsèques à la structure
%%% ***********************************************************************
% *** m = 0 -> z direction
    d_r3_tau_11 = diag(Delta_0i) + ((h-b)/h).*(K_ltau'*(diag(T_0l_prim_Rb)*K_ltau));
    d_r3_tau_12 = ((h-b)/h).*(diag(T_0l_tild_prim_Rb)*K_ltau)';
    
    d_r3_tau_21 = - diag(T_0l_prim_Rp)*K_ltau ;
    d_r3_tau_22 = ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_0j)*L_jtau))...
                 + ((e1-b)/(h-b)).*(M_ntau'*(diag(Pi_0n)*M_ntau)) - diag(T_0l_tild_prim_Rp);
   
    D_r3_tau = [d_r3_tau_11 d_r3_tau_12;...
                d_r3_tau_21 d_r3_tau_22];
% *** m = 1 -> x and beta direction
    d_r1_tau_11 = diag(Delta_1i) + ((h-b)/h).*(K_ltau'*(diag(T_1l_prim_Rb)*K_ltau));
    d_r1_tau_12 = ((h-b)/h).*(diag(T_1l_tild_prim_Rb)*K_ltau)';
    
    d_r1_tau_21 = - diag(T_1l_prim_Rp)*K_ltau ;
    d_r1_tau_22 = ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_1j)*L_jtau))...
                 + ((e1-b)/(h-b)).*(M_ntau'*(diag(Pi_1n)*M_ntau)) - diag(T_1l_tild_prim_Rp);
   
    D_r1_tau = [d_r1_tau_11 d_r1_tau_12;...
                d_r1_tau_21 d_r1_tau_22];

%%% ***********************************************************************
%%% ***********************************************************************

%%%                         Problème de RADIATION

%%% ***********************************************************************
%%% ***********************************************************************

%%% ***********************************************************************

%%%                         HEAVE MODE

%%% ***********************************************************************
if strcmp(options.dof,'all') || strcmp(options.dof,'heave')
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, column fixed
%%% ***********************************************************************    
    % % solution particulière pour r = R
    p_03 = ( (e1-b) / 12 ) * ( 2 - 3*(Rp^2 / (e1-b)^2) );
    p_n3 = sqrt(2).*f1(alpha_n,e1,b,e1,e1)./(2*(e1-b)^2);
    P_n3 = [p_03;p_n3];
    
    s_03 = ( 2*((h-b)^3-(h-e1)^3) - 3*Rp^2*(e1-b) ) / (12*(h-b)*(e1-b)) ;
    s_n3 = sqrt(2).*f1(alpha_n,e1,b,h,e1)./(2*(h-b)*(e1-b));
    S_n3 = [s_03;s_n3];
    
    r_03 = ( (h-b) / 12 ) * ( 2 - 3*(Rbo^2 / (h-b)^2) );
    r_l3 = sqrt(2).*f1(gamma_l,h,b,h,h)./(2*(h-b)^2);
    R_l3 = [r_03;r_l3];

    q_03 = ( (h-e2)^2 / (12*(h-b)) ) * ( 2 - 3*(Rp^2 / (h-e2)^2 ));
    q_j3 = sqrt(2).*f1(beta_j,h,e2,h,h)./(2*(h-b)*(h-e2));
    Q_j3 = [q_03;q_j3];
     
    % % Interface r = Rb
    h_r3_tau_1_k1 = ((h-b)/h) .* K_ltau'*(T_0l_prim_Rb.*R_l3) +  (.5*Rbo^2/h).*K_ltau(1,:)';
%     d_r3_tau_11 = diag(Delta_0i) + ((h-b)/h).*(K_ltau'*(diag(T_0l_prim_Rb)*K_ltau));
%     d_r3_tau_12 = ((h-b)/h).*(diag(T_0l_tild_prim_Rb)*K_ltau)';
    
    % % Interface r = Rp
    h_r3_tau_2_k1 = - T_0l_prim_Rp.*R_l3 - [(.5*Rp^2/(h-b));zeros(Nl,1)] - ((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*Q_j3))...
            + ((e1-b)/(h-b)).*(M_ntau'*(Pi_0n.*(P_n3-S_n3))) + (.5*Rp^2/(h-b)).*M_ntau(1,:)'; 
    
%     d_r3_tau_21 = - diag(T_0l_prim_Rp)*K_ltau ;
%     d_r3_tau_22 = ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_0j)*L_jtau))...
%                  + ((e1-b)/(h-b)).*(M_ntau'*(diag(Pi_0n)*M_ntau)) - diag(T_0l_tild_prim_Rp);
   
    H_r3_tau_k1 = [h_r3_tau_1_k1;h_r3_tau_2_k1];
%     D_r3_tau = [d_r3_tau_11 d_r3_tau_12;...
%                 d_r3_tau_21 d_r3_tau_22];
%     
    coeff = D_r3_tau\H_r3_tau_k1;

    A_r3_k1 = coeff(1:Ni+1);
    C_r3_k1 = coeff(Ni+2:end);
    
    B_r3_k1 = K_ltau*A_r3_k1 - R_l3;
    D_r3_k1 = M_ntau*C_r3_k1 - (P_n3-S_n3);
    F_r3_k1 = L_jtau*C_r3_k1 + Q_j3;
          
    tmp1 = T_0l_int(1)*B_r3_k1(1) + sqrt(2)*sum((-1).^l.*T_0l_int(2:end).*B_r3_k1(2:end));
    tmp2 = T_0l_tild_int(1)*C_r3_k1(1) + sqrt(2)*sum((-1).^l.*T_0l_tild_int(2:end).*C_r3_k1(2:end));
	tmp3 = .5*(h-b)*(.5*(Rbo^2 - Rp^2) - (1/8)*( (Rbo^4-Rp^4)/(h-b)^2 ));
    tmp4 = .5*Rp^2*D_r3_k1(1) + sqrt(2)*Rp*((-1).^n.*besseli(1,alpha_n.*Rp)./(alpha_n.*besseli(0,alpha_n.*Rp)))'*D_r3_k1(2:end);
    tmp5 = .5*(e1-b)*(.5*Rp^2 - (1/8)*( Rp^4/(e1-b)^2 ));
    z_r33_11 = 2*pi*rho*( tmp1 + tmp2 + tmp3 + tmp4 + tmp5);
   
    tmp1 = .5*Rp^2*F_r3_k1(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k1(2:end);
    tmp2 = .5*Rp^2*D_r3_k1(1) + sqrt(2)*Rp*( besseli(1,alpha_n.*Rp)./(alpha_n.*besseli(0,alpha_n.*Rp)) )'*D_r3_k1(2:end);
    z_r33_12 = 2*pi*rho*( tmp1 - (tmp2 - Rp^4/(16*(e1-b))) );
            
    
%%% ***********************************************************************    
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************    
    % % solution particulière pour r = Rc
    p_03 = - ( (e1-b) / 12 ) * ( 2 - 3*( Rp^2/(e1-b)^2 ) );
    p_n3 = - sqrt(2).*f1(alpha_n,e1,b,b,e1)./(2*(e1-b)^2);
    P_n3 = [p_03;p_n3];
    
    q_03 = ( (h-e2) / 12 ) * ( 2 - 3*(Rp/(h-e2))^2);
    q_j3 = sqrt(2).*f1(beta_j,h,e2,h,h)./(2*(h-e2)^2);
    Q_j3 = [q_03;q_j3];
    
    h_r3_tau_1_k2 = zeros(Ni+1,1);
    h_r3_tau_2_k2 = ((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*Q_j3)) + ((e1-b)/(h-b)).*(M_ntau'*(Pi_0n.*P_n3))...
            + (.5*Rp^2/(h-b)).*(L_jtau(1,:)'-M_ntau(1,:)');
    H_r3_tau_k2 = [h_r3_tau_1_k2;h_r3_tau_2_k2];
    
    % % Calcul des coefficients A_mi & C_ml
    coeff = D_r3_tau\H_r3_tau_k2;

    A_r3_k2 = coeff(1:Ni+1);
    C_r3_k2 = coeff(Ni+2:end);
   
    B_r3_k2 = K_ltau*A_r3_k2;
    D_r3_k2 = M_ntau*C_r3_k2 - P_n3;
    F_r3_k2 = L_jtau*C_r3_k2 - Q_j3;
    
    tmp1 = T_0l_int(1)*B_r3_k2(1) + sqrt(2)*sum((-1).^l.*T_0l_int(2:end).*B_r3_k2(2:end));
    tmp2 = T_0l_tild_int(1)*C_r3_k2(1) + sqrt(2)*sum((-1).^l.*T_0l_tild_int(2:end).*C_r3_k2(2:end));     
    tmp3 = .5*Rp^2*D_r3_k2(1) + sqrt(2)*Rp*((-1).^n.*besseli(1,alpha_n.*Rp)./(alpha_n.*besseli(0,alpha_n.*Rp)))'*D_r3_k2(2:end);
    z_r33_21 = 2*pi*rho*( tmp1 + tmp2 + tmp3 + Rp^4/(16*(e1-b)) );


    tmp1 = .5*(h-e2)*(.5*Rp^2 - (1/8)*(Rp^4/(h-e2)^2));
    tmp2 = .5*Rp^2*F_r3_k2(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k2(2:end);
    tmp3 = .5*Rp^2*D_r3_k2(1) + sqrt(2)*Rp*( besseli(1,alpha_n.*Rp)./(alpha_n.*besseli(0,alpha_n.*Rp)))'*D_r3_k2(2:end);
    tmp4 = - .5*(e1-b)*(.5*Rp^2 - (1/8)*( Rp^4/(e1-b)^2 ));
    z_r33_22 = 2*pi*rho*( tmp1 + tmp2 - (tmp3 + tmp4));
else
    
    z_r33_11 = 0;
    z_r33_12 = 0;
    z_r33_21 = 0;
    z_r33_22 = 0;
    
end% *** END OF HEAVE

%%% ***********************************************************************

%%%                         SURGE MODE

%%% ***********************************************************************
if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, column fixed
%%% ***********************************************************************

    % % Interface r = Rb
    h_r1_tau_1_k1 = -(Rbo/h).*I1';
%     d_r1_tau_11 = diag(Delta_1i) + ((h-b)/h).*(K_ltau'*(diag(T_1l_prim_Rb)*K_ltau));
%     d_r1_tau_12 = ((h-b)/h).*(diag(T_1l_tild_prim_Rb)*K_ltau)';
    
    % % Interface r = Rp
    h_r1_tau_2_k1 = zeros(Nl+1,1); 
    
%     d_r1_tau_21 = - diag(T_1l_prim_Rp)*K_ltau ;
%     d_r1_tau_22 = ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_1j)*L_jtau))...
%                  + ((e1-b)/(h-b)).*(M_ntau'*(diag(Pi_1n)*M_ntau)) - diag(T_1l_tild_prim_Rp);
   
    H_r1_tau_k1 = [h_r1_tau_1_k1;h_r1_tau_2_k1];
%     D_r1_tau = [d_r1_tau_11 d_r1_tau_12;...
%                 d_r1_tau_21 d_r1_tau_22];
    
    coeff = D_r1_tau\H_r1_tau_k1;

    A_r1_k1 = coeff(1:Ni+1);
    C_r1_k1 = coeff(Ni+2:end);
    
    B_r1_k1 = K_ltau*A_r1_k1;
    D_r1_k1 = M_ntau*C_r1_k1;
    F_r1_k1 = L_jtau*C_r1_k1;

    z_r11_11 = - pi*rho*Rbo*(I1*A_r1_k1);
    z_r11_12 = - pi*rho*Rp*(I2'*C_r1_k1);
    
    tmp1 = T_1l_int(1)*B_r1_k1(1) + sqrt(2)*sum((-1).^l.*T_1l_int(2:end).*B_r1_k1(2:end));
    tmp2 = T_1l_tild_int(1)*C_r1_k1(1) + sqrt(2)*sum((-1).^l.*T_1l_tild_int(2:end).*C_r1_k1(2:end));
    tmp3 = .25*D_r1_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^n.*besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r1_k1(2:end);
    z_r15_11 = - pi*rho*(Rbo*I3*A_r1_k1 + tmp1 + tmp2 + tmp3);
   
    tmp1 = .25*F_r1_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k1(2:end);
    tmp2 = .25*D_r1_k1(1)*Rp^3 + sqrt(2)*Rp^2*((besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r1_k1(2:end);   
    z_r15_12 = - pi*rho*( Rp*I4'*C_r1_k1 + tmp1 - tmp2);
    

%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************

    % % Interface r = Rb
    h_r1_tau_1_k2 = zeros(Ni+1,1);
    
    % % Interface r = Rp
    h_r1_tau_2_k2 = -(Rp/(h-b)).*I2; 
    
    H_r1_tau_k2 = [h_r1_tau_1_k2;h_r1_tau_2_k2];
    
    coeff = D_r1_tau\H_r1_tau_k2;

    A_r1_k2 = coeff(1:Ni+1);
    C_r1_k2 = coeff(Ni+2:end);
    
    B_r1_k2 = K_ltau*A_r1_k2;
    D_r1_k2 = M_ntau*C_r1_k2;
    F_r1_k2 = L_jtau*C_r1_k2;

    z_r11_21 = - pi*rho*Rbo*(I1*A_r1_k2);
    z_r11_22 = - pi*rho*Rp*(I2'*C_r1_k2);
    
    tmp1 = T_1l_int(1)*B_r1_k2(1) + sqrt(2)*sum((-1).^l.*T_1l_int(2:end).*B_r1_k2(2:end));
    tmp2 = T_1l_tild_int(1)*C_r1_k2(1) + sqrt(2)*sum((-1).^l.*T_1l_tild_int(2:end).*C_r1_k2(2:end));
    tmp3 = .25*D_r1_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^n.*besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r1_k2(2:end);
    z_r15_21 = - pi*rho*(Rbo*I3*A_r1_k2 + tmp1 + tmp2 + tmp3);
    
    tmp1 = .25*F_r1_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k2(2:end);
    tmp2 = .25*D_r1_k2(1)*Rp^3 + sqrt(2)*Rp^2*((besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r1_k2(2:end);   
    z_r15_22 = - pi*rho*( Rp*I4'*C_r1_k2 + tmp1 - tmp2);
else
    z_r11_11 = 0;
    z_r11_12 = 0;
    z_r15_11 = 0;
    z_r15_12 = 0;
    
    z_r11_21 = 0;
    z_r11_22 = 0;
    z_r15_21 = 0;
    z_r15_22 = 0;
end% *** END OF SURGE

%%% ***********************************************************************

%%%                         PITCH MODE

%%% ***********************************************************************
if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, plate fixed
%%% *********************************************************************** 


    % % solution particulière pour r = R
    p_05 = (Rp^3/(24*(e1-b))) * ( 3 - 4*((e1-b)/Rp)^2);
    p_n5 = -((sqrt(2)*Rp)/(2*(e1-b)^2)).*f1(alpha_n,e1,b,e1,e1);
    P_n5 = [p_05;p_n5];
    
    s_05 =  (Rp^3/(24*(e1-b)*(h-b)))*( 3*(e1-b) - 4*((h-b)^3-(h-e1)^3)/Rp^2 ) ;
    s_n5 = -((sqrt(2)*Rp)/(2*(e1-b)*(h-b))).*f1(alpha_n,e1,b,h,e1);
    S_n5 = [s_05;s_n5];
    
    r_05 = (Rbo^3/(24*(h-b))) * ( 3 - 4*((h-b)/Rbo)^2 );
    r_l5 = -((sqrt(2)*Rbo)/(2*(h-b)^2)).*f1(gamma_l,h,b,h,h);
    R_l5 = [r_05;r_l5];

    q_05 = (Rp^3/(24*(h-b))) * ( 3 - 4*((h-e2)/Rp)^2);
    q_j5 = -((sqrt(2)*Rp)/(2*(h-b)*(h-e2))).*f1(beta_j,h,e2,h,h);
    Q_j5 = [q_05;q_j5];
    
    h_r5_tau_1_k1 = ((h-b)/h).*K_ltau'*(T_1l_prim_Rb.*R_l5) - ((3*Rbo^3)/(8*h)).*K_ltau(1,:)' + (Rbo/(2*h*(h-b))).*I7' - (Rbo/h).*I3';
    
    h_r5_tau_2_k1 = - T_1l_prim_Rp.*R_l5 + [3*Rp^3/(8*(h-b));zeros(Nl,1)] - (.5*Rp/(h-b)^2).*I6...
            + ((e1-b)/(h-b)).*(M_ntau'*(Pi_1n.*(P_n5-S_n5))) - ((3*Rp^3)/(8*(h-b))).*M_ntau(1,:)' + (Rp/(2*(e1-b)*(h-b))).*I5_1...
            - ((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*Q_j5)); 
    

    H_r5_tau_k1 = [h_r5_tau_1_k1;h_r5_tau_2_k1];
    D_r5_tau = D_r1_tau;
    
    coeff = D_r5_tau\H_r5_tau_k1;

    A_r5_k1 = coeff(1:Ni+1);
    C_r5_k1 = coeff(Ni+2:end);
    
    B_r5_k1 = K_ltau*A_r5_k1 - R_l5;
    D_r5_k1 = M_ntau*C_r5_k1 - (P_n5 - S_n5);
    F_r5_k1 = L_jtau*C_r5_k1 + Q_j5;

    
    tmp1 = T_1l_int(1)*B_r5_k1(1) + sqrt(2)*sum((-1).^l.*T_1l_int(2:end).*B_r5_k1(2:end));
    tmp2 = T_1l_tild_int(1)*C_r5_k1(1) + sqrt(2)*sum((-1).^l.*T_1l_tild_int(2:end).*C_r5_k1(2:end));
    tmp3 = ( (Rbo^6-Rp^6)/6 - (Rbo^4-Rp^4)*(h-b)^2 ) / (8*(h-b));
    tmp4 = .25*D_r5_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^n.*besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r5_k1(2:end);
    tmp5 = (Rp^6/6 - Rp^4*(e1-b)^2) / (8*(e1-b));
    z_r55_11 = - pi*rho*(Rbo*I3*A_r5_k1 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5);
    
    tmp1 = (.5*Rp^3*(e1^2-e2^2) - 4*Rp*( e2*(h-e2)^3 - e1*(h-e1)^3 - ((h-e1)^4-(h-e2)^4)/4 )/3 ) / (8*(h-b));
    tmp2 = .25*F_r5_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k1(2:end);
    tmp3 = .25*D_r5_k1(1)*Rp^3 + sqrt(2)*Rp^2*((besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r5_k1(2:end);
    tmp4 = Rp^6/(48*(e1-b));
    z_r55_12 = - pi*rho*( Rp*(I4'*C_r5_k1 + tmp1) + tmp2 - (tmp3 + tmp4) );
    
    z_r51_11 = - pi*rho*Rbo*(I1*A_r5_k1);
    
    tmp1 = (Rp^3*(e2-e1) - 4*Rp*((h-e1)^3-(h-e2)^3)/3)/(8*(h-b));
    z_r51_12 = - pi*rho*Rp*(I2'*C_r5_k1 + tmp1);


%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************   
    p_05 = -(Rp^3/(24*(e1-b))) * ( 3 - 4*((e1-b)/Rp)^2);
    p_n5 = ((sqrt(2)*Rp)/(2*(e1-b)^2)).*f1(alpha_n,e1,b,b,e1);
    P_n5 = [p_05;p_n5];

    q_05 = (Rp^3/(24*(h-e2))) * ( 3 - 4*((h-e2)/Rp)^2);
    q_j5 = -((sqrt(2)*Rp)/(2*(h-e2)^2)).*f1(beta_j,h,e2,h,h);
    Q_j5 = [q_05;q_j5];
    
    h_r5_tau_1_k2 = zeros(Ni+1,1);
    h_r5_tau_2_k2 = ((e1-b)/(h-b)).*(M_ntau'*(Pi_1n.*P_n5)) + ((3*Rp^3)/(8*(h-b))).*M_ntau(1,:)' - (Rp/(2*(e1-b)*(h-b))).*I5_2...
                        + ((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*Q_j5)) - ((3*Rp^3)/(8*(h-b))).*L_jtau(1,:)' + (Rp/(2*(h-e2)*(h-b))).*I8 - (Rp/(h-b)).*I4;
    

    H_r5_tau_k2 = [h_r5_tau_1_k2;h_r5_tau_2_k2];
    
    coeff = D_r5_tau\H_r5_tau_k2;

    A_r5_k2 = coeff(1:Ni+1);
    C_r5_k2 = coeff(Ni+2:end);
    
    B_r5_k2 = K_ltau*A_r5_k2;
    D_r5_k2 = M_ntau*C_r5_k2 - P_n5 ;
    F_r5_k2 = L_jtau*C_r5_k2 - Q_j5;

    tmp1 = T_1l_int(1)*B_r5_k2(1) + sqrt(2)*sum((-1).^l.*T_1l_int(2:end).*B_r5_k2(2:end));
    tmp2 = T_1l_tild_int(1)*C_r5_k2(1) + sqrt(2)*sum((-1).^l.*T_1l_tild_int(2:end).*C_r5_k2(2:end));
    tmp3 = .25*D_r5_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^n.*besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r5_k2(2:end);
    tmp4 = -(Rp^6/6) / (8*(e1-b));
    z_r55_21 = - pi*rho*(Rbo*I3*A_r5_k2 + tmp1 + tmp2 + tmp3 + tmp4);
    
    tmp1 = .25*F_r5_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k2(2:end);
    tmp2 = (Rp^6/6 - Rp^4*(h-e2)^2) / (8*(h-e2));
    tmp3 = .25*D_r5_k2(1)*Rp^3 + sqrt(2)*Rp^2*((besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*D_r5_k2(2:end);
    tmp4 = -(Rp^6/6 - Rp^4*(e1-b)^2) / (8*(e1-b));
    z_r55_22 = - pi*rho*( Rp*I4'*C_r5_k2 + tmp1 + tmp2 - (tmp3 + tmp4) );

    z_r51_21 = - pi*rho*Rbo*(I1*A_r5_k2);
    z_r51_22 = - pi*rho*Rp*(I2'*C_r5_k2);
else
    z_r55_11 = 0;
    z_r55_12 = 0;
    z_r51_11 = 0;
    z_r51_12 = 0;
    
    z_r55_21 = 0;
    z_r55_22 = 0;
    z_r51_21 = 0;
    z_r51_22 = 0;
end% *** END OF PITCH   
%%% ***********************************************************************
%%% ***********************************************************************
%%%             Problème de DIFFRACTION
%%% ***********************************************************************
%%% ***********************************************************************
    
    a0 = 1i*omega;
    
    Z0_lambda0 = (N_lambda_i(1)^(-.5))*cosh(lambda_0*h);
    B = (-1i*g) / (omega*Z0_lambda0);

    if strcmp(options.dof,'all') || strcmp(options.dof,'heave')
        Chi_0 = -lambda_0*Rbo*besselj(1,lambda_0*Rbo);

        h_d_0tau_1 = B.*( [Chi_0;zeros(Ni,1)] -  besselj(0,lambda_0*Rbo)*((h-b)/h).*(K_ltau'*(T_0l_prim_Rb.*K_ltau(:,1))) );
        h_d_0tau_2 = B*besselj(0,lambda_0*Rbo).*(T_0l_prim_Rp.*K_ltau(:,1));

%         d_d_0tau_11 = B.*d_r3_tau_11;
%         d_d_0tau_21 = B.*d_r3_tau_21;

        H_d_0tau = [h_d_0tau_1;h_d_0tau_2];
%         D_d_0tau = [d_d_0tau_11 d_r3_tau_12;...
%                     d_d_0tau_21 d_r3_tau_22];

        coeff = D_r3_tau\H_d_0tau;

        a_0i = coeff(1:Ni+1);
        c_0l = coeff(Ni+2:end);

        b_0l = B.*besselj(0,lambda_0*Rbo).*K_ltau(:,1) + K_ltau*a_0i;
        d_0n = M_ntau*c_0l;
        f_0j = L_jtau*c_0l;    

        tmp1 = T_0l_int(1)*b_0l(1) + sqrt(2)*sum((-1).^l.*T_0l_int(2:end).*b_0l(2:end));
        tmp2 = T_0l_tild_int(1)*c_0l(1) + sqrt(2)*sum((-1).^l.*T_0l_tild_int(2:end).*c_0l(2:end)); 
        tmp3 = .5*Rp^2*d_0n(1) + sqrt(2)*Rp*((-1).^n.*besseli(1,alpha_n.*Rp)./(alpha_n.*besseli(0,alpha_n.*Rp)) )'*d_0n(2:end);
        fz_1 = 2*pi*a0*rho*(tmp1 + tmp2 + tmp3);

        tmp1 = .5*Rp^2*f_0j(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*f_0j(2:end);
        tmp2 = .5*Rp^2*d_0n(1) + sqrt(2)*Rp*( besseli(1,alpha_n.*Rp)./(alpha_n.*besseli(0,alpha_n.*Rp)) )'*d_0n(2:end);
        fz_2 = 2*pi*a0*rho*(tmp1 - tmp2);

        % *** HASKIND (GREEN's Second theorem on Sb)
        if options.Haskind
            step = .0001;
            V_r = Rp:step:Rbo;

            C = 2*pi*rho*g*A*lambda_0/cosh(lambda_0*h);

            tmp1 = ( C*cosh(lambda_0*(h-b))/lambda_0^2 )*Rbo*besselj(1,lambda_0*Rbo);

            cst1 = C*sinh(lambda_0*(h-b));
            tmp21 = (e1-b)*Rp*besselj(1,lambda_0*Rp)/(2*lambda_0);
            tmp22 = (1/(4*(e1-b)*lambda_0))*( Rp^3*besselj(1,lambda_0*Rp) - 2*Rp^2*besselj(2,lambda_0*Rp)/lambda_0);
            tmp23 = D_r3_k1(1)*Rp*besselj(1,lambda_0*Rp)/lambda_0;
            tmp24 = (lambda_0*Rp.*besseli(0,alpha_n.*Rp).*besselj(1,lambda_0*Rp)...
                        + alpha_n.*Rp.*besseli(1,alpha_n.*Rp).*besselj(0,lambda_0*Rp))./ (alpha_n.^2 + lambda_0^2);
            tmp25 = sqrt(2)*sum((-1).^n.*D_r3_k1(2:end).*tmp24./besseli(0,alpha_n.*Rp));
            tmp2 = cst1*(tmp21 - tmp22 + tmp23 + tmp25);

            cst2 = C*sinh(lambda_0*(h-e1));
            tmp35 = sqrt(2)*sum(D_r3_k1(2:end).*tmp24./besseli(0,alpha_n.*Rp));
            tmp3 = cst2*(-tmp22 + tmp23 + tmp35);

            % ***
            tmp41 = ((h-b)/(2*lambda_0))*( Rbo*besselj(1,lambda_0*Rbo) - Rp*besselj(1,lambda_0*Rp) );
            tmp42 = (1/(4*(h-b)))*(sum(V_r.^3.*besselj(0,lambda_0.*V_r))*step);
            tmp43 = B_r3_k1(1)*(sum(log(V_r./Rp).*besselj(0,lambda_0.*V_r).*V_r)*step)/log(Rbo/Rp);
            tmp44 = C_r3_k1(1)*(sum(log(Rbo./V_r).*besselj(0,lambda_0.*V_r).*V_r)*step)/log(Rbo/Rp);

            cst = besseli(0,gamma_l.*Rbo).*besselk(0,gamma_l.*Rp) - besseli(0,gamma_l.*Rp).*besselk(0,gamma_l.*Rbo);
            tmp451 = besselk(0,gamma_l.*Rp).*sum(besseli(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step...
                - besseli(0,gamma_l.*Rp).*sum(besselk(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step;
            tmp452 = besseli(0,gamma_l.*Rbo).*sum(besselk(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step...
                - besselk(0,gamma_l.*Rbo).*sum(besseli(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step;
            tmp45 = sqrt(2)*sum((-1).^l.*( B_r3_k1(2:end).*tmp451 + C_r3_k1(2:end).*tmp452 )./cst);
            tmp4 = cst1*(tmp41 - tmp42 + tmp43 + tmp44 + tmp45);
            % *** 

            cst4 = C*sinh(lambda_0*(h-e2));
            tmp51 = F_r3_k1(1)*Rp*besselj(1,lambda_0*Rp)/lambda_0;
            tmp52 = (lambda_0*Rp.*besseli(0,beta_j.*Rp).*besselj(1,lambda_0*Rp)...
                        + beta_j.*Rp.*besseli(1,beta_j.*Rp).*besselj(0,lambda_0*Rp))./ (beta_j.^2 + lambda_0^2);
            tmp53 = sqrt(2)*sum((-1).^j.*F_r3_k1(2:end).*tmp52./besseli(0,beta_j.*Rp));
            tmp5 = cst4*(tmp51 + tmp53);

            tmp6 = C*Rbo*besselj(1,lambda_0*Rbo)*(N_lambda_i.^(-.5)*(A_r3_k1.*integraleInterface( Lambda_i, Lambda_i(1),b, 0, h, h, epsilon )));
            tmp71 = (C_r3_k1(1) - Rp^2/(4*(h-b)))*integraleInterface(0,Lambda_i(1),e2,e1,0,h,epsilon);
            tmp72 = sqrt(2)*sum( C_r3_k1(2:end).*integraleInterface(gamma_l,Lambda_i(1),e2,e1,h,h,epsilon) );
            tmp73 = f1(Lambda_i(1),e2,e1,h,h)/(2*(h-b));
            tmp7 = C*Rp*besselj(1,lambda_0*Rp)*(tmp71 + tmp72 + tmp73);
            fz_1_Haskind = ( tmp1 - (tmp2 - tmp3 + tmp4 + tmp5 + tmp6 + tmp7 ) ) ;

        else
            fz_1_Haskind = 0;
        end% *** END OF HASKIND
    else
        fz_1 = 0;
        fz_2 = 0;

        fz_1_Haskind = 0;
    end
    
    if ~strcmp(options.dof,'heave')
        Chi_1 = 2*1i*(lambda_0*Rbo*besselj(0,lambda_0*Rbo)-besselj(1,lambda_0*Rbo));
        
        h_d_1tau_1 = B.*( [Chi_1;zeros(Ni,1)] -  (2*1i*besselj(1,lambda_0*Rbo))*((h-b)/h).*(K_ltau'*(T_1l_prim_Rb.*K_ltau(:,1))) );
        h_d_1tau_2 = B*(2*1i*besselj(1,lambda_0*Rbo)).*(T_1l_prim_Rp.*K_ltau(:,1));

        H_d_1tau = [h_d_1tau_1;h_d_1tau_2];

    %     d_d_1tau_11 = B.*d_r1_tau_11;
    %     d_d_1tau_21 = B.*d_r1_tau_21;

    %     D_d_1tau = [d_d_1tau_11 d_r1_tau_12;...
    %                 d_d_1tau_21 d_r1_tau_22];

        coeff = D_r1_tau\H_d_1tau;   

        a_1i = coeff(1:Ni+1);
        c_1l = coeff(Ni+2:end);

        b_1l = B*2*1i*besselj(1,lambda_0*Rbo).*K_ltau(:,1) + K_ltau*a_1i;
        d_1n = M_ntau*c_1l;
        f_1j = L_jtau*c_1l;
        
    end

    if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
        fx_1 = -B*pi*a0*rho*Rbo*( 2*1i*besselj(1,lambda_0*Rbo)*I1(1) + I1*a_1i );
        fx_2 = -pi*a0*rho*Rp*( I2'*c_1l );
    else
        fx_1 = 0;
        fx_2 = 0;
    end

    if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
        tmp1 = T_1l_int(1)*b_1l(1) + sqrt(2)*sum((-1).^l.*T_1l_int(2:end).*b_1l(2:end));
        tmp2 = T_1l_tild_int(1)*c_1l(1) + sqrt(2)*sum((-1).^l.*T_1l_tild_int(2:end).*c_1l(2:end));
        tmp3 = .25*d_1n(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^n.*besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*d_1n(2:end);
        ty_1 = - pi*a0*rho*( B*Rbo*(2*1i*besselj(1,lambda_0*Rbo)*I3(1) + I3*a_1i) + tmp1 + tmp2 + tmp3 );

        tmp1 = .25*f_1j(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*f_1j(2:end);
        tmp2 = .25*d_1n(1)*Rp^3 + sqrt(2)*Rp^2*((besseli(2,alpha_n.*Rp))./(alpha_n.*besseli(1,alpha_n.*Rp)))'*d_1n(2:end);   
        ty_2 = - pi*a0*rho*( Rp*I4'*c_1l + tmp1 - tmp2);
    else
        ty_1 = 0;
        ty_2 = 0;
    end
   

%%% ***********************************************************************
%%% ***********************************************************************


% % z1 = -h:0.001:0;
% % z2 = -h:0.001:-b;
% % z3 = -e1:0.001:-b;
% % z4 = -h:0.001:-e2;
% % 
% % % % Condition sur la dérivée du potentiel de vitesse dans Omega I.
% % 
% % Z_i_I = zeros(Ni+1,length(z1));
% % for i = 1 : Ni+1
% %     Z_i_I(i,:) = N_lambda_i(i)^(-.5).*cos(Lambda_i(i).*(z1+h));
% % end
% % 
% % z_0_II = ones(1,length(z2));
% % z_l_II = zeros(Nl,length(z2));
% % for i = 1 : Nl
% %     z_l_II(i,:) = sqrt(2).*cos(gamma_l(i).*(z2+h));
% % end
% % Z_l_II = [z_0_II;z_l_II];
% % 
% % z_0_III = ones(1,length(z3));
% % z_n_III = zeros(Nn,length(z3));
% % for i = 1 : Nn
% %     z_n_III(i,:) = sqrt(2).*cos(alpha_n(i).*(z3+e1));
% % end
% % Z_n_III = [z_0_III;z_n_III];
% % 
% % z_0_IV = ones(1,length(z4));
% % z_j_IV = zeros(Nj,length(z4));
% % for i = 1 : Nj
% %     z_j_IV(i,:) = sqrt(2).*cos(beta_j(i).*(z4+h));
% % end
% % Z_j_IV = [z_0_IV;z_j_IV];
% % 
% % 
% % phi_1_Rb_r5_k1 = conj(A_r5_k1)'*Z_i_I;
% % d_phi_1_Rb_r5_k1 =  - (conj(A_r5_k1)'.*Delta_1i)*Z_i_I;
% % 
% % phi_2_Rb_r5_k1 = conj(B_r5_k1)'*Z_l_II + (Rb^3 - 4*Rb.*(z2+h).^2)./(8*(h-b)) ;
% % d_phi_2_Rb_r5_k1 = conj(B_r5_k1.*T_1l_prim_Rb + C_r5_k1.*T_1l_tild_prim_Rb)'*Z_l_II + (3*Rb^3 - 4*Rb.*(z2+h).^2)./(8*(h-b)) ;
% % 
% % phi_2_Rp_r5_k1 = conj(C_r5_k1)'*Z_l_II + (Rp^3 - 4*Rp.*(z2+h).^2)./(8*(h-b)) ;
% % d_phi_2_Rp_r5_k1 = conj(B_r5_k1.*T_1l_prim_Rp + C_r5_k1.*T_1l_tild_prim_Rp)'*Z_l_II + (3*Rp^3 - 4*Rp.*(z2+h).^2)./(8*(h-b)) ;
% % 
% % phi_3_Rp_r5_k1 = conj(D_r5_k1)'*Z_n_III + (Rp^3 - 4*Rp.*(z3+e1).^2)./(8*(e1-b));
% % d_phi_3_Rp_r5_k1 = conj(D_r5_k1.*Pi_1n)'*Z_n_III + (3*Rp^3 - 4*Rp.*(z3+e1).^2)./(8*(e1-b)) ;
% % 
% % phi_4_Rp_r5_k1 = conj(F_r5_k1)'*Z_j_IV;
% % d_phi_4_Rp_r5_k1 = (conj(F_r5_k1).*Gamma_1j)'*Z_j_IV;


% % phi_1_Rb_r5_k2 = conj(A_r5_k2)'*Z_i_I;
% % d_phi_1_Rb_r5_k2 =  - (conj(A_r5_k2)'.*Delta_1i)*Z_i_I;
% % 
% % phi_2_Rb_r5_k2 = conj(B_r5_k2)'*Z_l_II;
% % d_phi_2_Rb_r5_k2 = conj(B_r5_k2.*T_1l_prim_Rb + C_r5_k2.*T_1l_tild_prim_Rb)'*Z_l_II;
% % 
% % phi_2_Rp_r5_k2 = conj(C_r5_k2)'*Z_l_II;
% % d_phi_2_Rp_r5_k2 = conj(B_r5_k2.*T_1l_prim_Rp + C_r5_k2.*T_1l_tild_prim_Rp)'*Z_l_II;
% % 
% % phi_3_Rp_r5_k2 = conj(D_r5_k2)'*Z_n_III - (Rp^3 - 4*Rp.*(z3+b).^2)./(8*(e1-b));
% % d_phi_3_Rp_r5_k2 = conj(D_r5_k2.*Pi_1n)'*Z_n_III - (3*Rp^3 - 4*Rp.*(z3+b).^2)./(8*(e1-b)) ;
% % 
% % phi_4_Rp_r5_k2 = conj(F_r5_k2)'*Z_j_IV + (Rp^3 - 4*Rp.*(z4+h).^2)./(8*(h-e2));
% % d_phi_4_Rp_r5_k2 = (conj(F_r5_k2).*Gamma_1j)'*Z_j_IV + (3*Rp^3 - 4*Rp.*(z4+h).^2)./(8*(h-e2));



% % figure(1); hold on;
% % subplot(2,1,1), hold on;
% % plot(z2,abs(phi_2_Rb_r5_k1),'r','LineWidth',1);
% % plot(z1,abs(phi_1_Rb_r5_k1),'--');
% % subplot(2,1,2), hold on;
% % plot(z2,abs(d_phi_2_Rb_r5_k1),'r','LineWidth',1);
% % plot(z1,abs(d_phi_1_Rb_r5_k1),'--');
% % 
% % 
% % figure(2); hold on;
% % subplot(2,1,1), hold on;
% % plot(z4,abs(phi_4_Rp_r5_k1),'r','LineWidth',1);
% % plot(z3,abs(phi_3_Rp_r5_k1),'r','LineWidth',1);
% % plot(z2,abs(phi_2_Rp_r5_k1),'--');
% % subplot(2,1,2), hold on;
% % plot(z4,abs(d_phi_4_Rp_r5_k1),'r','LineWidth',1);
% % plot(z3,abs(d_phi_3_Rp_r5_k1),'r','LineWidth',1);
% % plot(z2,abs(d_phi_2_Rp_r5_k1),'--');
% % 
% % figure(1); hold on;
% % subplot(2,1,1), hold on;
% % plot(z2,abs(phi_2_Rb_r5_k2),'r','LineWidth',1);
% % plot(z1,abs(phi_1_Rb_r5_k2),'--');
% % subplot(2,1,2), hold on;
% % plot(z2,abs(d_phi_2_Rb_r5_k2),'r','LineWidth',1);
% % plot(z1,abs(d_phi_1_Rb_r5_k2),'--');
% % 
% % 
% % figure(2); hold on;
% % subplot(2,1,1), hold on;
% % plot(z4,abs(phi_4_Rp_r5_k2),'r','LineWidth',1);
% % plot(z3,abs(phi_3_Rp_r5_k2),'r','LineWidth',1);
% % plot(z2,abs(phi_2_Rp_r5_k2),'--');
% % subplot(2,1,2), hold on;
% % plot(z4,abs(d_phi_4_Rp_r5_k2),'r','LineWidth',1);
% % plot(z3,abs(d_phi_3_Rp_r5_k2),'r','LineWidth',1);
% % plot(z2,abs(d_phi_2_Rp_r5_k2),'--');



    out.fz_1 = fz_1;
    out.fz_2 = fz_2;
    
    out.fz_1_Haskind = fz_1_Haskind;
    out.fz_2_Haskind = 0;
    
    out.fx_1 = fx_1;
    out.fx_2 = fx_2;
    
    out.fx_1_Haskind = 0;
    out.fx_2_Haskind = 0;
    
    out.ty_1 = ty_1;
    out.ty_2 = ty_2;
    
    out.ty_1_Haskind = 0;
    out.ty_2_Haskind = 0;
    
    out.z_r33_11 = z_r33_11;
    out.z_r33_12 = z_r33_12;
    out.z_r33_22 = z_r33_22;
    out.z_r33_21 = z_r33_21;
    
    out.z_r11_11 = z_r11_11;
    out.z_r11_12 = z_r11_12;
    out.z_r11_22 = z_r11_22;
    out.z_r11_21 = z_r11_21;
    
    out.z_r15_11 = z_r15_11;
    out.z_r15_12 = z_r15_12;
    out.z_r15_22 = z_r15_22;
    out.z_r15_21 = z_r15_21;
    
    out.z_r55_11 = z_r55_11;
    out.z_r55_12 = z_r55_12;
    out.z_r55_22 = z_r55_22;
    out.z_r55_21 = z_r55_21;
    
    out.z_r51_11 = z_r51_11;
    out.z_r51_12 = z_r51_12;
    out.z_r51_22 = z_r51_22;
    out.z_r51_21 = z_r51_21; 


end