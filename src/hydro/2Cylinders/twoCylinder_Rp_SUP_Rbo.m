function out = twoCylinder_Rp_sup_a_Rb(EnvParam,WECstructure,omega,Lambda_i,Gamma_l,options)

    A = EnvParam.A;
    g = EnvParam.g;
    rho = EnvParam.rho;    
    h = EnvParam.h;
    
    Rb = WECstructure.Rb;
    Rp = WECstructure.Rp;
    
    b = WECstructure.b;
    e1 = WECstructure.e1;
    e2 = WECstructure.e2;
    
    Ni = options.Truncate.Ni;
    Nn = options.Truncate.Nn;
    Nj = options.Truncate.Nj;
    Nl = options.Truncate.Nl;
    
    Zc = options.Zc;
    
    epsilon = options.epsilon;

    n = (1:1:Nn)';
    j = (1:1:Nj)';
    
    beta_j = (pi/(h-e2)).*j;
    alpha_n = (pi/(e1-b)).*n;   
   
% %     options=optimset('Display','off');   % NO option to display output
% %     lambda_0 = fsolve(@(x)(omega^2 - x*g*tanh(x*h)),omega^2/g,options);
% %     gamma_0 = fsolve(@(x)(omega^2 - x*g*tanh(x*e1)),omega^2/g,options);
% % 
% %     cst = lambda_0.*tanh(lambda_0*h);
% % 
% %     lambda_i = zeros(1,Ni);
% %     options=optimset('Display','off');   % NO option to display output
% %     for i = 1 : Ni
% %         lambda_i(1,i) = fsolve(@(x)(cst + x*tan(x*h)),(pi/2+(i-1)*pi)*(1/h),options);
% %     end
% %     
% %     Lambda_i = [-1i*lambda_0 lambda_i];

    Lambda_i = Lambda_i(1:Ni+1);
    Gamma_l = Gamma_l(1:Nl+1);
    
    lambda_0 = -imag(Lambda_i(1));
    lambda_i = Lambda_i(2:Ni+1);
    
    gamma_0 = -imag(Gamma_l(1));
    gamma_l = Gamma_l(2:Nl+1);
    
    N_lambda_i = 0.5.*( 1 + sin(2.*Lambda_i.*h)./(2.*Lambda_i.*h));
    n_lambda_0 = N_lambda_i(1);
    n_lambda_i = N_lambda_i(2:end);

% %     cst = kappa_0.*tanh(kappa_0*e1);
    
% %     gamma_l = zeros(Nl,1);
% %     options=optimset('Display','off');   % NO option to display output
% %     for l = 1 : Nl
% %         kappa_l(l,1) = fsolve(@(x)(cst + x*tan(x*e1)),(pi/2+(l-1)*pi)*(1/e1),options);
% %     end
% %     
%     Gamma_l = [-1i*kappa_0;kappa_l];
    
    N_gamma_l =.5.*( 1 + sin(2.*Gamma_l.*e1)./(2.*Gamma_l.*e1));
    n_gamma_0 = N_gamma_l(1);
    n_gamma_l = N_gamma_l(2:end);
    
    L_jlambda0 = ( (-1).^[0;j].*n_lambda_0^(-.5)*lambda_0*sinh(lambda_0*(h-e2)) )...
                    ./ ( (h-e2).*(lambda_0^2 + [0;beta_j].^2) );
    L_jlambda_i = zeros(Nj,Ni);
    for i = 1 : Ni
        for p = 1 : Nj
            if (beta_j(p)/lambda_i(i) < 1+epsilon)&&(beta_j(p)/lambda_i(i) > 1-epsilon)
                L_jlambda_i(p,i) = n_lambda_i(i)^(-.5)/2;
            else
                L_jlambda_i(p,i) = ( (-1)^p*n_lambda_i(i)^(-.5)*lambda_i(i)*sin(lambda_i(i)*(h-e2)) )...
                    / ( (h-e2)*(lambda_i(i)^2 - beta_j(p)^2) );
            end
        end
    end
    L_0lambda_i = (n_lambda_i.^(-.5).*sin(lambda_i.*(h-e2))) ./ ((h-e2).*lambda_i);
    L_jtau = [L_jlambda0(1) L_0lambda_i;L_jlambda0(2:end) L_jlambda_i];
    L_jtau(2:end,:) = sqrt(2).*L_jtau(2:end,:);


    M_nkappa_0 = ( (-1).^[0;n].*n_gamma_0^(-.5)*gamma_0*sinh(gamma_0*(e1-b)) )...
                    ./ ( (e1-b).*(gamma_0^2 + [0;alpha_n].^2) );
    M_nkappa_l = zeros(Nn,Nl);
    for i = 1 : Nl
        for p = 1 : Nn
            if (alpha_n(p)/gamma_l(i) < 1+epsilon)&&(alpha_n(p)/gamma_l(i) > 1-epsilon)
                M_nkappa_l(p,i) = n_gamma_l(i)^(-.5)/2;
            else
                M_nkappa_l(p,i) = ( (-1)^p*n_gamma_l(i)^(-.5)*gamma_l(i)*sin(gamma_l(i)*(e1-b)) )...
                    / ( (e1-b)*(gamma_l(i)^2 - alpha_n(p)^2) );
            end
        end
    end
    M_0kappa_l = (n_gamma_l'.^(-.5).*sin(gamma_l'.*(e1-b))) ./ ((e1-b).*gamma_l');
    M_ntau = [M_nkappa_0(1) M_0kappa_l;M_nkappa_0(2:end) M_nkappa_l];
    M_ntau(2:end,:) = sqrt(2).*M_ntau(2:end,:);
    
    K_ltau = zeros(Nl+1,Ni+1);
    for i = 1 : Ni+1
        for p = 1 : Nl+1
            if (Gamma_l(p)/Lambda_i(i) < 1+epsilon)&&(Gamma_l(p)/Lambda_i(i) > 1-epsilon)
                K_ltau(p,i) = (N_lambda_i(i)^(-.5)*N_gamma_l(p).^(-.5))*...
                    .5*((sin(Lambda_i(i)*(h+e1)) - sin(Lambda_i(i)*(h-e1)))/(2*Lambda_i(i)) + e1*cos(Lambda_i(i)*(h-e1)))/e1;
%             end % SO -- modification du 23/06/2014
                  % rajout du "else"
            else
                K_ltau(p,i) = (N_lambda_i(i)^(-.5)*N_gamma_l(p).^(-.5))*...
                    ( Gamma_l(p)*(cos(Lambda_i(i)*h)*sin(Gamma_l(p)*e1))...
                     - Lambda_i(i)*(sin(Lambda_i(i)*h)*cos(Gamma_l(p)*e1) - sin(Lambda_i(i)*(h-e1))) )/...
                     (e1*(Gamma_l(p)^2 - Lambda_i(i)^2));
            end
        end
    end
    
    
    gamma_00 = 0;
    gamma_0j = beta_j.*Rp.*(besseli(1,beta_j.*Rp) ./ besseli(0,beta_j.*Rp));
    
    gamma_10 = 1;
    gamma_1j = beta_j.*Rp.*(besseli(0,beta_j.*Rp) ./ besseli(1,beta_j.*Rp)) - 1;

    pi_00 = 0;
    pi_0n = alpha_n.*Rb.*(besseli(1,alpha_n.*Rb) ./ besseli(0,alpha_n.*Rb));
    
    pi_10 = 1;
    pi_1n = alpha_n.*Rb.*(besseli(0,alpha_n.*Rb) ./ besseli(1,alpha_n.*Rb)) - 1;
    
    delta_00 = lambda_0*Rp*(besselh(1,lambda_0*Rp) / besselh(0,lambda_0*Rp));
    delta_0i = lambda_i.*Rp.*(besselk(1,lambda_i.*Rp) ./ besselk(0,lambda_i.*Rp));
    
    delta_10 = 1 - lambda_0*Rp*(besselh(0,lambda_0*Rp) / besselh(1,lambda_0*Rp));
    delta_1i = 1 + lambda_i.*Rp.*( besselk(0,lambda_i.*Rp) ./ besselk(1,lambda_i.*Rp) );

    Gamma_0j = [gamma_00;gamma_0j];
    Gamma_1j = [gamma_10;gamma_1j];
    
    Delta_0i = [delta_00 delta_0i];
    Delta_1i = [delta_10 delta_1i];

    Pi_0n = [pi_00;pi_0n];
    Pi_1n = [pi_10;pi_1n];

    [S_prim_0l_Rp,S_tild_prim_0l_Rp] = Sprim_func(0,Rp,Rp,Rb,Gamma_l);
    [S_prim_0l_Rb,S_tild_prim_0l_Rb] = Sprim_func(0,Rb,Rp,Rb,Gamma_l);
    
    [S_prim_1l_Rp,S_tild_prim_1l_Rp] = Sprim_func(1,Rp,Rp,Rb,Gamma_l);
    [S_prim_1l_Rb,S_tild_prim_1l_Rb] = Sprim_func(1,Rb,Rp,Rb,Gamma_l);

    [S_0l_int,S_0l_tild_int] = Sint_func(0,Rb,Rp,Gamma_l);
    [S_1l_int,S_1l_tild_int] = Sint_func(1,Rb,Rp,Gamma_l);
    

    I1 = N_gamma_l.^(-.5).*( sin(Gamma_l.*e1)-sin(Gamma_l.*(e1-b)) ) ./ Gamma_l;                % % Int_(-b)^(0){ Z_l^II }
    I2 = N_lambda_i.^(-.5).*( sin(Lambda_i.*(h-e1)) - sin(Lambda_i.*(h-e2)) ) ./ Lambda_i;      % % Int_(-e2)^(-e1){ Z_i^I }
    I3_1 = N_gamma_l.^(-.5).*f1(Gamma_l,e1,b,e1,e1);                                            % % Int_(-e1)^(-b){ (z+e1)^2*Z_l^II }
    I3_2 = N_gamma_l.^(-.5).*f1(Gamma_l,e1,b,b,e1);                                             % % Int_(-e1)^(-b){ (z+b)^2*Z_l^II }
    I4 = N_lambda_i.^(-.5).*f1(Lambda_i,h,e2,h,h);                                              % % Int_(-h)^(-e2){ (z+h)^2*Z_i^I }
    I5 = N_gamma_l.^(-.5).*f2(Gamma_l,-b,0,Zc(1),e1);                                           % % Int_(-b)^(0){ (z-Zc)*Z_l^II }
    I6 = N_lambda_i.^(-.5).*f2(Lambda_i,-e2,-e1,Zc(2),h);                                       % % Int_(-e2)^(-e1){ (z-Zc)*Z_i^I }
    I7 = N_gamma_l.^(-.5).*f2(Gamma_l,-e1,0,0,e1);                                              % % Int_(-e1)^(0){ z*Z_l^II }
    
%%% ***********************************************************************
%%% ***         Définition des matrices intrinsèques à la structure
%%% ***********************************************************************
% *** m = 0 -> z direction
    d_r3_tau_11 = diag(Delta_0i) + ((h-e2)/h).*(L_jtau'*(diag(Gamma_0j)*L_jtau))...
                        + (e1/h).*(K_ltau'*(diag(S_prim_0l_Rp)*K_ltau));
    d_r3_tau_12 = (e1/h).*( diag(S_tild_prim_0l_Rp)*K_ltau )';

    d_r3_tau_21 = - diag(S_prim_0l_Rb)*K_ltau;
    d_r3_tau_22 = ((e1-b)/e1).*(M_ntau'*(diag(Pi_0n)*M_ntau)) - diag(S_tild_prim_0l_Rb);

    D_r3_tau = [d_r3_tau_11 d_r3_tau_12;...
                d_r3_tau_21 d_r3_tau_22];
% *** m = 1 -> x and beta direction
    d_r1_tau_11 = diag(Delta_1i) + ((h-e2)/h).*(L_jtau'*(diag(Gamma_1j)*L_jtau))...
                        + (e1/h).*(K_ltau'*(diag(S_prim_1l_Rp)*K_ltau));
    d_r1_tau_12 = (e1/h).*( diag(S_tild_prim_1l_Rp)*K_ltau )';
    d_r1_tau_21 = - diag(S_prim_1l_Rb)*K_ltau;
    d_r1_tau_22 = ((e1-b)/e1).*(M_ntau'*(diag(Pi_1n)*M_ntau)) - diag(S_tild_prim_1l_Rb);

    D_r1_tau = [d_r1_tau_11 d_r1_tau_12;...
                d_r1_tau_21 d_r1_tau_22];
%%% ***********************************************************************
%%% ***********************************************************************

%%%                      Problème de RADIATION

%%% ***********************************************************************
%%% ***********************************************************************

%%%                         HEAVE MODE

%%% ***********************************************************************
if strcmp(options.dof,'all') || strcmp(options.dof,'heave')
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, column fixed
%%% ***********************************************************************    
    % % solution particulière pour r = R
    p_03 = ( (e1-b) / 12 ) * ( 2 - 3*(Rb^2 / (e1-b)^2) );
    p_n3 = (sqrt(2)*(e1-b).*(-1).^n) ./ (n.^2 .* pi^2);
    P_n3 = [p_03;p_n3];   
    
    h_r3_tau_1_k1 = zeros(Ni+1,1);
%     d_r3_tau_11 = diag(Delta_0i) + ((h-e2)/h).*(L_jtau'*(diag(Gamma_0j)*L_jtau))...
%                         + (e1/h).*(K_ltau'*(diag(S_prim_0l_Rp)*K_ltau));
%     d_r3_tau_12 = (e1/h).*( diag(S_tild_prim_0l_Rp)*K_ltau )';
% 
%     d_r3_tau_21 = - diag(S_prim_0l_Rb)*K_ltau;
%     d_r3_tau_22 = ((e1-b)/e1).*(M_ntau'*(diag(Pi_0n)*M_ntau)) - diag(S_tild_prim_0l_Rb);
  
    h_r3_tau_2_k1 = ((e1-b)/e1).*M_ntau'*(Pi_0n.*P_n3) +  (.5*Rb^2/e1).*M_ntau(1,:)';
   
    H_r3_tau_k1 = [h_r3_tau_1_k1;...
                h_r3_tau_2_k1];

%     D_r3_tau = [d_r3_tau_11 d_r3_tau_12;...
%                 d_r3_tau_21 d_r3_tau_22];
    
    coeff = D_r3_tau\H_r3_tau_k1;

    A_r3_k1 = coeff(1:Ni+1);
    C_r3_k1 = coeff(Ni+2:end);
    
    B_r3_k1 = K_ltau*A_r3_k1;
    D_r3_k1 = M_ntau*C_r3_k1 - P_n3;
    F_r3_k1 = L_jtau*A_r3_k1;

	tmp1 = .5*(e1-b)*(.5*Rb^2 - (1/8)*( Rb^4/(e1-b)^2 ));
    tmp2 = .5*Rb^2*D_r3_k1(1) + sqrt(2)*Rb*((-1).^n.*besseli(1,alpha_n.*Rb)./(alpha_n.*besseli(0,alpha_n.*Rb)))'*D_r3_k1(2:end);
    z_r33_11 = 2*pi*rho*(tmp1 + tmp2);

    tmp1 = .5*Rp^2*F_r3_k1(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k1(2:end);
    tmp2 = sum(N_gamma_l.^(-.5).*S_0l_int.*B_r3_k1);
    tmp3 = sum(N_gamma_l.^(-.5).*S_0l_tild_int.*C_r3_k1);
    tmp4 = .5*Rb^2*D_r3_k1(1) + sqrt(2)*Rb*( besseli(1,alpha_n.*Rb)./(alpha_n.*besseli(0,alpha_n.*Rb)) )'*D_r3_k1(2:end);
    z_r33_12 = 2*pi*rho*( tmp1 - (tmp2 + tmp3 + tmp4 - Rb^4/(16*(e1-b)) ) );
           
    
%%% ***********************************************************************    
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************    
    % % solution particulière pour r = Rc
    p_03 = - ( (e1-b) / 12 ) * ( 2 - 3*( Rb^2/(e1-b)^2) );
    p_n3 = - (sqrt(2)*(e1-b)) ./ (n.^2 .* pi^2);
    P_n3 = [p_03;p_n3];
    
    q_03 = ( (h-e2) / 12 ) * ( 2 - 3*(Rp/(h-e2))^2);
    q_j3 = (sqrt(2)*(h-e2).*(-1).^j) ./ (j.^2 .* pi^2);
    Q_j3 = [q_03;q_j3];
    
    R_l3 = ( I7 + (g/omega^2).*N_gamma_l.^(-.5).*sin(Gamma_l.*e1)./Gamma_l ) ./ e1;    

    s_03 = .5*(b^2-e1^2)/(e1-b) + (g/omega^2) ;
    s_n3 = sqrt(2).*f2(alpha_n,-e1,-b,0,e1)./(e1-b);
    S_n3 = [s_03;s_n3];


    h_r3_tau_1_k2 = ((h-e2)/h).*(L_jtau'*(Gamma_0j.*Q_j3)) + (e1/h).*(K_ltau'*(S_prim_0l_Rp.*R_l3)) + (.5*Rp^2/h).*L_jtau(1,:)';
    h_r3_tau_2_k2 = -S_prim_0l_Rb.*R_l3 + ((e1-b)/e1).*M_ntau'*(Pi_0n.*(P_n3-S_n3)) - (.5*Rb^2/e1).*M_ntau(1,:)';    
    H_r3_tau_k2 = [h_r3_tau_1_k2;h_r3_tau_2_k2];
    
    % % Calcul des coefficients A_mi, C_mn & E_mn
    coeff = D_r3_tau\H_r3_tau_k2;

    A_r3_k2 = coeff(1:Ni+1);
    C_r3_k2 = coeff(Ni+2:end);
   
    B_r3_k2 = K_ltau*A_r3_k2 - R_l3;
    D_r3_k2 = M_ntau*C_r3_k2 - P_n3 + S_n3;
    F_r3_k2 = L_jtau*A_r3_k2 - Q_j3;
    
    tmp1 = .5*Rb^2*D_r3_k2(1) + sqrt(2)*Rb*((-1).^n.*besseli(1,alpha_n.*Rb)./(alpha_n.*besseli(0,alpha_n.*Rb)))'*D_r3_k2(2:end);
    z_r33_21 = 2*pi*rho*( tmp1 + Rb^4/(16*(e1-b)) );

    
    tmp1 = .5*Rp^2*F_r3_k2(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k2(2:end);
    tmp2 = .5*(h-e2)*(.5*Rp^2 - .125*(Rp^4/(h-e2)^2));
    tmp3 = sum(N_gamma_l.^(-.5).*S_0l_int.*B_r3_k2);
    tmp4 = sum(N_gamma_l.^(-.5).*S_0l_tild_int.*C_r3_k2);
    tmp5 = .5*(g/omega^2 - e1)*(Rp^2-Rb^2);
    tmp6 = .5*Rb^2*D_r3_k2(1) + sqrt(2)*Rb*( besseli(1,alpha_n.*Rb)./(alpha_n.*besseli(0,alpha_n.*Rb)) )'*D_r3_k2(2:end);
    tmp7 = - .5*(e1-b)*(.5*Rb^2 - .125*(Rb^4/(e1-b)^2));
    z_r33_22 = 2*pi*rho*( tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7) );
    
    fz_1_Haskind = -(4*1i*rho*g*h*sqrt(n_lambda_0)*A_r3_k1(1))...
                            / (cosh(lambda_0*h)*besselh(0,lambda_0*Rp));
    fz_2_Haskind = -(4*1i*rho*g*h*sqrt(n_lambda_0)*A_r3_k2(1))...
                            / (cosh(lambda_0*h)*besselh(0,lambda_0*Rp));

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

    h_r1_tau_1_k1 = zeros(Ni+1,1);
%     d_r1_tau_11 = diag(Delta_1i) + ((h-e2)/h).*(L_jtau'*(diag(Gamma_1j)*L_jtau))...
%                         + (e1/h).*(K_ltau'*(diag(S_prim_1l_Rp)*K_ltau));
%     d_r1_tau_12 = (e1/h).*( diag(S_tild_prim_1l_Rp)*K_ltau )';

    h_r1_tau_2_k1 = -(Rb/e1).*I1;
%     d_r1_tau_21 = - diag(S_prim_1l_Rb)*K_ltau;
%     d_r1_tau_22 = ((e1-b)/e1).*(M_ntau'*(diag(Pi_1n)*M_ntau)) - diag(S_tild_prim_1l_Rb);
  
   
    H_r1_tau_k1 = [h_r1_tau_1_k1;h_r1_tau_2_k1];
%     D_r1_tau = [d_r1_tau_11 d_r1_tau_12;...
%                 d_r1_tau_21 d_r1_tau_22];
    
    coeff = D_r1_tau\H_r1_tau_k1;

    A_r1_k1 = coeff(1:Ni+1);
    C_r1_k1 = coeff(Ni+2:end);
    
    B_r1_k1 = K_ltau*A_r1_k1;
    D_r1_k1 = M_ntau*C_r1_k1;
    F_r1_k1 = L_jtau*A_r1_k1;

    z_r11_11 = -pi*rho*Rb*(I1'*C_r1_k1);
    z_r11_12 = -pi*rho*Rp*(I2*A_r1_k1);
    
    tmp1 = .25*D_r1_k1(1)*Rb^3 + sqrt(2)*Rb^2*(((-1).^n.*besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r1_k1(2:end);
    z_r15_11 = - pi*rho*(Rb*(I5'*C_r1_k1) + tmp1);

    tmp1 = .25*F_r1_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k1(2:end);
    tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r1_k1);
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r1_k1);
    tmp4 = .25*D_r1_k1(1)*Rb^3 + sqrt(2)*Rb^2*((besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r1_k1(2:end);
    z_r15_12 = - pi*rho*( Rp*I6*A_r1_k1 + tmp1 - (tmp2 + tmp3 + tmp4));
    
%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************

    h_r1_tau_1_k2 =  -(Rp/h).*I2';
    h_r1_tau_2_k2 = zeros(Nl+1,1);
    H_r1_tau_k2 = [h_r1_tau_1_k2;h_r1_tau_2_k2];

    coeff = D_r1_tau\H_r1_tau_k2;
    
    A_r1_k2 = coeff(1:Ni+1);
    C_r1_k2 = coeff(Ni+2:end);
    
    B_r1_k2 = K_ltau*A_r1_k2;
    D_r1_k2 = M_ntau*C_r1_k2;
    F_r1_k2 = L_jtau*A_r1_k2;
    
    z_r11_21 = -pi*rho*Rb*(I1'*C_r1_k2);
    z_r11_22 = -pi*rho*Rp*(I2*A_r1_k2);
    
    tmp1 = .25*D_r1_k2(1)*Rb^3 + sqrt(2)*Rb^2*(((-1).^n.*besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r1_k2(2:end);
    z_r15_21 = - pi*rho*(Rb*(I5'*C_r1_k2) + tmp1);
 
    tmp1 = .25*F_r1_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k2(2:end);
    tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r1_k2);
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r1_k2);
    tmp4 = .25*D_r1_k2(1)*Rb^3 + sqrt(2)*Rb^2*((besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r1_k2(2:end);
    z_r15_22 = - pi*rho*( Rp*I6*A_r1_k2 + tmp1 - (tmp2 + tmp3 + tmp4));
    
    fx_1_Haskind = (4*rho*g*h*sqrt(n_lambda_0)*A_r1_k1(1))...
                            / (cosh(lambda_0*h)*besselh(1,lambda_0*Rp));
    fx_2_Haskind = (4*rho*g*h*sqrt(n_lambda_0)*A_r1_k2(1))...
                            / (cosh(lambda_0*h)*besselh(1,lambda_0*Rp));

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

%%%                             PITCH MODE

%%% ***********************************************************************
if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
%%% ***********************************************************************
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, plate fixed
%%% *********************************************************************** 
    p_05 = (Rb^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rb^2));
    p_n5 = -((sqrt(2)*Rb)/(2*(e1-b)^2)).*f1(alpha_n,e1,b,e1,e1);
    P_n5 = [p_05;p_n5];
    
    h_r5_tau_1_k1 = zeros(Ni+1,1);
    h_r5_tau_2_k1 = ((e1-b)/e1).*M_ntau'*(Pi_1n.*P_n5) - ((3*Rb^3)/(8*e1)).*M_ntau(1,:)' + (Rb/(2*e1*(e1-b))).*I3_1 - (Rb/e1).*I5;

    H_r5_tau_k1 = [h_r5_tau_1_k1;h_r5_tau_2_k1];

    D_r5_tau = D_r1_tau;
    
    coeff = D_r5_tau\H_r5_tau_k1;

    A_r5_k1 = coeff(1:Ni+1);
    C_r5_k1 = coeff(Ni+2:end);
    
    B_r5_k1 = K_ltau*A_r5_k1;
    D_r5_k1 = M_ntau*C_r5_k1 - P_n5;
    F_r5_k1 = L_jtau*A_r5_k1;

    tmp1 =  ( 1/(8*(e1-b)) ) * ( Rb^6/6 - Rb^4*(e1-b)^2);
    tmp2 = .25*D_r5_k1(1)*Rb^3 + sqrt(2)*Rb^2*(((-1).^n.*besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r5_k1(2:end);
    z_r55_11 = - pi*rho*(Rb*(I5'*C_r5_k1) + tmp1 + tmp2);

    tmp1 = .25*F_r5_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k1(2:end);
    tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r5_k1);
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r5_k1);
    tmp4 = .25*D_r5_k1(1)*Rb^3 + sqrt(2)*Rb^2*((besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r5_k1(2:end);
    tmp5 = Rb^6/(48*(e1-b));
    z_r55_12 = - pi*rho*( Rp*I6*A_r5_k1 + tmp1 - (tmp2 + tmp3 + tmp4 + tmp5));
    
    z_r51_11 = -pi*rho*Rb*(I1'*C_r5_k1);
    z_r51_12 = -pi*rho*Rp*(I2*A_r5_k1);

%%% ***********************************************************************
% %     Cas k = 2       Buoy fixed, plate oscillating 
%%% ***********************************************************************

    p_05 = -(Rb^3/(24*(e1-b))) * ( 3 - 4*((e1-b)^2/Rb^2));
    p_n5 = ((sqrt(2)*Rb)/(2*(e1-b)^2)).*f1(alpha_n,e1,b,b,e1);
    P_n5 = [p_05;p_n5];
    
    s_05 = - (Rb/(e1-b))*(.5*(b^2-e1^2) + (g/omega^2)*(e1-b));
    s_n5 = - (sqrt(2)*Rb/(e1-b))*f2(alpha_n,-e1,-b,0,e1);
    S_n5 = [s_05;s_n5];
    
    q_05 = (Rp^3/(24*(h-e2))) * ( 3 - 4*((h-e2)^2/Rp^2));
    q_j5 = -((sqrt(2)*Rp)/(2*(h-e2)^2)).*f1(beta_j,h,e2,h,h);
    Q_j5 = [q_05;q_j5];
    
    R_l5 = - (Rp/e1).*(I7 + (g/omega^2).*N_gamma_l.^(-.5).*sin(Gamma_l.*e1) ./ Gamma_l);
    
    tmp1 = N_lambda_i.^(-.5).*f2(Lambda_i,-e1,0,0,h);
    tmp2 = (g/omega^2).*N_lambda_i.^(-.5).*( sin(Lambda_i.*h) - sin(Lambda_i.*(h-e1)) ) ./ Lambda_i;
    h_r5_tau_1_k2 =  ((h-e2)/h).*(L_jtau'*(Gamma_1j.*Q_j5)) - ((3*Rp^3)/(8*h)).*L_jtau(1,:)' + (.5*Rp/(h*(h-e2))).*I4'...
        + (e1/h).*(K_ltau'*(S_prim_1l_Rp.*R_l5)) + (Rp/h).*(tmp1+tmp2)' - (Rp/h).*I6';
    
    tmp2 = (g/omega^2).*N_gamma_l.^(-.5).*sin(Gamma_l.*e1) ./ Gamma_l;
    h_r5_tau_2_k2 = ((e1-b)/e1).*M_ntau'*(Pi_1n.*(P_n5-S_n5)) + ((3*Rb^3)/(8*e1)).*M_ntau(1,:)' - (.5*Rb/(e1*(e1-b))).*I3_2...
                                - S_prim_1l_Rb.*R_l5 - (Rb/e1).*(I7 + tmp2) ;
    H_r5_tau_k2 = [h_r5_tau_1_k2;h_r5_tau_2_k2];
    
    coeff = D_r5_tau\H_r5_tau_k2;
    
    A_r5_k2 = coeff(1:Ni+1);
    C_r5_k2 = coeff(Ni+2:end);
    
    B_r5_k2 = K_ltau*A_r5_k2 - R_l5;
    D_r5_k2 = M_ntau*C_r5_k2 - P_n5 + S_n5;
    F_r5_k2 = L_jtau*A_r5_k2 - Q_j5;

    tmp1 = Rb*(-b^3/3 + .5*(g/omega^2)*b^2);
    tmp2 = .25*D_r5_k2(1)*Rb^3 + sqrt(2)*Rb^2*(((-1).^n.*besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r5_k2(2:end);
    tmp3 = - Rb^6 / (48*(e1-b));
    z_r55_21 = - pi*rho*( Rb*(I5'*C_r5_k2 + tmp1) + tmp2 + tmp3);
    
    tmp1 = .25*F_r5_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k2(2:end);
    tmp2 = ( 1/(8*(h-e2)) ) * ( Rp^6/6 - Rp^4*(h-e2)^2);
    tmp3 = sum(N_gamma_l.^(-.5).*S_1l_int.*B_r5_k2);
    tmp4 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*C_r5_k2);
    tmp5 = -.25*(g/omega^2 - e1)*(Rp^4-Rb^4);
    tmp6 = .25*D_r5_k2(1)*Rb^3 + sqrt(2)*Rb^2*((besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*D_r5_k2(2:end);
    tmp7 = - ( 1/(8*(e1-b)) ) * ( Rb^6/6 - Rb^4*(e1-b)^2);
    z_r55_22 = - pi*rho*( Rp*I6*A_r5_k2 + tmp1 + tmp2 - (tmp3 + tmp4 + tmp5 + tmp6 + tmp7));
    
    tmp1 = Rb*(.5*b^2 - (g/omega^2)*b);
    z_r51_21 = - pi*rho*Rb*(I1'*C_r5_k2 + tmp1);
    
    z_r51_22 = - pi*rho*Rp*(I2*A_r5_k2);
    
    ty_1_Haskind = (4*rho*g*h*sqrt(n_lambda_0)*A_r5_k1(1))...
                            / (cosh(lambda_0*h)*besselh(1,lambda_0*Rp));
    ty_2_Haskind = (4*rho*g*h*sqrt(n_lambda_0)*A_r5_k2(1))...
                            / (cosh(lambda_0*h)*besselh(1,lambda_0*Rp));

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
    
    Z0_lambda0 = (n_lambda_0^(-.5))*cosh(lambda_0*h);
    B = (-1i*g*A) / (omega*Z0_lambda0);

    if strcmp(options.dof,'all') || strcmp(options.dof,'heave')
        Chi_0 = -lambda_0*Rp*besselj(1,lambda_0*Rp);
        
        h_d_0tau_1 = B.*( [Chi_0;zeros(Ni,1)] -  besselj(0,lambda_0*Rp).*...
                    ( (e1/h).*(K_ltau'*(S_prim_0l_Rp.*K_ltau(:,1)))...
                    + ((h-e2)/h).*(L_jtau'*(Gamma_0j.*L_jtau(:,1))) ) );
        h_d_0tau_2 = B*besselj(0,lambda_0*Rp).*(S_prim_0l_Rb.*K_ltau(:,1));

    %     d_d_0tau_11 = B.*d_r3_tau_11;
    %     d_d_0tau_21 = B.*d_r3_tau_21;

        H_d_0tau = [h_d_0tau_1;h_d_0tau_2];
    %     D_d_0tau = [d_d_0tau_11 d_r3_tau_12;...
    %                 d_d_0tau_21 d_r3_tau_22];

        coeff = D_r3_tau\H_d_0tau;      % -> ATTENTION MODIFICATION 16/05/2013

        a_0i = coeff(1:Ni+1);
        c_0l = coeff(Ni+2:end);

        b_0l = B.*besselj(0,lambda_0*Rp).*K_ltau(:,1) + K_ltau*a_0i;    % -> ATTENTION MODIFICATION 16/05/2013
        d_0n = M_ntau*c_0l;
        f_0j = B.*besselj(0,lambda_0*Rp).*L_jtau(:,1) + L_jtau*a_0i;    % -> ATTENTION MODIFICATION 16/05/2013
    
        tmp1 = sqrt(2)*Rb*((-1).^n.*besseli(1,alpha_n.*Rb)./(alpha_n.*besseli(0,alpha_n.*Rb)))'*d_0n(2:end);
        fz_1 = 2*pi*a0*rho*( .5*Rb^2*d_0n(1) + tmp1);

        tmp1 = .5*Rp^2*f_0j(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*f_0j(2:end);
        tmp2 = .5*Rb^2*d_0n(1) + sqrt(2)*Rb*( besseli(1,alpha_n.*Rb)./(alpha_n.*besseli(0,alpha_n.*Rb)) )'*d_0n(2:end);
        tmp3 = sum(N_gamma_l.^(-.5).*S_0l_int.*b_0l);
        tmp4 = sum(N_gamma_l.^(-.5).*S_0l_tild_int.*c_0l);
        fz_2 = 2*pi*a0*rho*(tmp1 - tmp2 - (tmp3 + tmp4));

% % % % %         % *** HASKIND (GREEN's Second theorem on Sb)
% % % % %         if options.Haskind
% % % % % 
% % % % %             step = .001;
% % % % %             V_r = Rb:step:Rp;
% % % % % 
% % % % %             C = 2*pi*rho*g*A*lambda_0/cosh(lambda_0*h);
% % % % % 
% % % % %             tmp1 = ( C*cosh(lambda_0*(h-b))/lambda_0^2 )*Rb*besselj(1,lambda_0*Rb);
% % % % % 
% % % % %             cst1 = C*sinh(lambda_0*(h-b));
% % % % %             tmp21 = (e1-b)*Rb*besselj(1,lambda_0*Rb)/(2*lambda_0);
% % % % %             tmp22 = (1/(4*(e1-b)*lambda_0))*( Rb^3*besselj(1,lambda_0*Rb) - 2*Rb^2*besselj(2,lambda_0*Rb)/lambda_0);
% % % % %             tmp23 = D_r3_k1(1)*Rb*besselj(1,lambda_0*Rb)/lambda_0;
% % % % %             tmp24 = (lambda_0*Rb.*besseli(0,alpha_n.*Rb).*besselj(1,lambda_0*Rb)...
% % % % %                         + alpha_n.*Rb.*besseli(1,alpha_n.*Rb).*besselj(0,lambda_0*Rb))./ (alpha_n.^2 + lambda_0^2);
% % % % %             tmp25 = sqrt(2)*sum((-1).^n.*D_r3_k1(2:end).*tmp24./besseli(0,alpha_n.*Rb));
% % % % %             tmp2 = cst1*(tmp21 - tmp22 + tmp23 + tmp25);
% % % % % 
% % % % %             cst2 = C*sinh(lambda_0*(h-e1));
% % % % %             tmp35 = sqrt(2)*sum(D_r3_k1(2:end).*tmp24./besseli(0,alpha_n.*Rb));
% % % % %             tmp3 = cst2*(-tmp22 + tmp23 + tmp35);
% % % % % 
% % % % %             tmp41 = besselh(0,gamma_0.*Rp).*besselh(0,2,gamma_0.*Rb) - besselh(0,gamma_0.*Rb).*besselh(0,2,gamma_0.*Rp);
% % % % %             tmp42 = besseli(0,gamma_l.*Rp).*besselk(0,gamma_l.*Rb) - besseli(0,gamma_l.*Rb).*besselk(0,gamma_l.*Rp);
% % % % %             tmp431 = besselh(0,2,gamma_0*Rb)*sum(besselh(0,gamma_0.*V_r).*besselj(0,lambda_0.*V_r).*V_r)*step...
% % % % %                         - besselh(0,gamma_0*Rb)*sum(besselh(0,2,gamma_0.*V_r).*besselj(0,lambda_0.*V_r).*V_r)*step;
% % % % %             tmp432 = besselh(0,gamma_0*Rp)*sum(besselh(0,2,gamma_0.*V_r).*besselj(0,lambda_0.*V_r).*V_r)*step...
% % % % %                         - besselh(0,2,gamma_0*Rp)*sum(besselh(0,gamma_0.*V_r).*besselj(0,lambda_0.*V_r).*V_r)*step;
% % % % %             tmp43 = N_gamma_l(1)^(-.5)*(B_r3_k1(1)*tmp431 + C_r3_k1(1)*tmp432)/tmp41;
% % % % %             tmp441 = besselk(0,gamma_l.*Rb).*sum(besseli(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step...
% % % % %                         - besseli(0,gamma_l.*Rb).*sum(besselk(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step;
% % % % %             tmp442 = besseli(0,gamma_l.*Rp).*sum(besselk(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step...
% % % % %                         - besselk(0,gamma_l.*Rp).*sum(besseli(0,gamma_l*V_r).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r).*V_r)),2).*step;
% % % % %             tmp44 = sum(N_gamma_l(2:end).^(-.5).*(B_r3_k1(2:end).*tmp441 + C_r3_k1(2:end).*tmp442)./tmp42);
% % % % %             tmp4 = cst2*(tmp43 + tmp44);
% % % % % 
% % % % %             cst4 = C*sinh(lambda_0*(h-e2));
% % % % %             tmp51 = F_r3_k1(1)*Rp*besselj(1,lambda_0*Rp)/lambda_0;
% % % % %             tmp52 = (lambda_0*Rp.*besseli(0,beta_j.*Rp).*besselj(1,lambda_0*Rp)...
% % % % %                         + beta_j.*Rp.*besseli(1,beta_j.*Rp).*besselj(0,lambda_0*Rp))./ (beta_j.^2 + lambda_0^2);
% % % % %             tmp53 = sqrt(2)*sum((-1).^j.*F_r3_k1(2:end).*tmp52./besseli(0,beta_j.*Rp));
% % % % %             tmp5 = cst4*(tmp51 + tmp53);
% % % % % 
% % % % %             tmp6 = C*Rb*besselj(1,lambda_0*Rb)*(N_gamma_l'.^(-.5)*(C_r3_k1.*integraleInterface( Gamma_l, Lambda_i(1),b, 0, e1, h, epsilon )));
% % % % %             tmp7 = C*Rp*besselj(1,lambda_0*Rp)*(N_lambda_i.^(-.5)*(A_r3_k1.*integraleInterface( Lambda_i, Lambda_i(1),e2, e1, h, h, epsilon )));
% % % % %             fz_1_Haskind = 1*( tmp1 - (tmp2 - (tmp3 + tmp4) + tmp5 + tmp6 + tmp7 ) ) ;
% % % % %         else
% % % % %             fz_1_Haskind = 0;
% % % % %         end% *** END OF HASKIND
    else
        fz_1 = 0;
        fz_2 = 0;
        
        fz_1_Haskind = 0;
    end
    
    if ~strcmp(options.dof,'heave')
        Chi_1 = 2*1i*(lambda_0*Rp*besselj(0,lambda_0*Rp)-besselj(1,lambda_0*Rp));

        h_d_1tau_1 = B.*( [Chi_1;zeros(Ni,1)] -  (2*1i*besselj(1,lambda_0*Rp)).*...
                    ( (e1/h).*(K_ltau'*(S_prim_1l_Rp.*K_ltau(:,1)))...
                    + ((h-e2)/h).*(L_jtau'*(Gamma_1j.*L_jtau(:,1))) ) );
        h_d_1tau_2 = B*(2*1i*besselj(1,lambda_0*Rp)).*(S_prim_1l_Rb.*K_ltau(:,1));

    %     d_d_1tau_11 = B.*d_r1_tau_11;
    %     d_d_1tau_21 = B.*d_r1_tau_21;

        H_d_1tau = [h_d_1tau_1;h_d_1tau_2];
    %     D_d_1tau = [d_d_1tau_11 d_r1_tau_12;...
    %                 d_d_1tau_21 d_r1_tau_22];

        coeff = D_r1_tau\H_d_1tau;      % -> ATTENTION MODIFICATION 16/05/2013

        a_1i = coeff(1:Ni+1);
        c_1l = coeff(Ni+2:end);

        b_1l = B*2*1i*besselj(1,lambda_0*Rp).*K_ltau(:,1) + K_ltau*a_1i;    % -> ATTENTION MODIFICATION 16/05/2013
        d_1n = M_ntau*c_1l;
        f_1j = B*2*1i*besselj(1,lambda_0*Rp).*L_jtau(:,1) + L_jtau*a_1i;    % -> ATTENTION MODIFICATION 16/05/2013
    end
    
    if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
        fx_1 = pi*a0*rho*Rb*( I1'*c_1l );
        fx_2 = -pi*a0*rho*Rp*( 2*1i*B*besselj(1,lambda_0*Rp)*I2(1) + I2*a_1i );
    else
        fx_1 = 0;
        fx_2 = 0;
    end

    if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
        tmp1 = .25*d_1n(1)*Rb^3 + sqrt(2)*Rb^2*(((-1).^n.*besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*d_1n(2:end);
        ty_1 = - pi*a0*rho*( Rb*I5'*c_1l + tmp1);

        tmp1 = .25*f_1j(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*f_1j(2:end);
        tmp2 = sum(N_gamma_l.^(-.5).*S_1l_int.*b_1l);
        tmp3 = sum(N_gamma_l.^(-.5).*S_1l_tild_int.*c_1l);
        tmp4 = .25*d_1n(1)*Rb^3 + sqrt(2)*Rb^2*((besseli(2,alpha_n.*Rb))./(alpha_n.*besseli(1,alpha_n.*Rb)))'*d_1n(2:end);   
        ty_2 = - pi*a0*rho*( Rp*(B*2*1i*besselj(1,lambda_0*Rp)*I6(1) + I6*a_1i) + tmp1 - (tmp2 + tmp3 + tmp4));
    else
        ty_1 = 0;
        ty_2 = 0;
    end

    
% %%% ***********************************************************************

% z1 = -h:0.001:0;
% z2 = -e1:.001:0;
% z3 = -e1:0.001:-b;
% z4 = -h:0.001:-e2;
% 
% % % % Condition sur la dérivée du potentiel de vitesse dans Omega I.
% Z_i_I = zeros(Ni+1,length(z1));
% for i = 1 : Ni
%     Z_i_I(i,:) = N_lambda_i(i)^(-.5).*cos(Lambda_i(i).*(z1+h));
% end
% 
% Z_l_II = zeros(Nl+1,length(z2));
% for i = 1 : Nl
%     Z_l_II(i,:) = N_kappa_l(i)^(-.5).*cos(Kappa_l(i).*(z2+e1));
% end
% 
% z_0_III = ones(1,length(z3));
% z_n_III = zeros(Nn,length(z3));
% for i = 1 : Nn
%     z_n_III(i,:) = sqrt(2).*cos(alpha_n(i).*(z3+e1));
% end
% Z_n_III = [z_0_III;z_n_III];
% 
% z_0_IV = ones(1,length(z4));
% z_j_IV = zeros(Nj,length(z4));
% for i = 1 : Nj
%     z_j_IV(i,:) = sqrt(2).*cos(beta_j(i).*(z4+h));
% end
% Z_j_IV = [z_0_IV;z_j_IV];
% 
% 
% phi_1_Rp_r5_k2 = conj(A_r5_k2)'*Z_i_I;
% d_phi_1_Rp_r5_k2 =  - (conj(A_r5_k2)'.*Delta_1i)*Z_i_I;
% 
% phi_2_Rp_r5_k2 = conj(B_r5_k2)'*Z_l_II - Rp.*(z2 + (g/omega^2));
% d_phi_2_Rp_r5_k2 = conj(B_r5_k2.*S_prim_1l_Rp + C_r5_k2.*S_tild_prim_1l_Rp)'*Z_l_II - Rp.*(z2 + (g/omega^2));
% 
% phi_4_Rp_r5_k2 = conj(F_r5_k2)'*Z_j_IV + (1/(8*(h-e2))).*( Rp^3 - 4*Rp.*(z4+h).^2 );
% d_phi_4_Rp_r5_k2 = (conj(F_r5_k2).*Gamma_1j)'*Z_j_IV + (1/(8*(h-e2))).*( 3*Rp^3 - 4*Rp.*(z4+h).^2 );
% 
% phi_2_Rb_r5_k2 = conj(C_r5_k2)'*Z_l_II - Rb.*(z2 + (g/omega^2));
% d_phi_2_Rb_r5_k2 = conj(B_r5_k2.*S_prim_1l_Rb + C_r5_k2.*S_tild_prim_1l_Rb)'*Z_l_II - Rb.*(z2 + (g/omega^2));
% 
% phi_3_Rb_r5_k2 = conj(D_r5_k2)'*Z_n_III - (1/(8*(e1-b))).*( Rb^3 - 4*Rb.*(z3+b).^2 );
% d_phi_3_Rb_r5_k2 = conj(D_r5_k2.*Pi_1n)'*Z_n_III - (1/(8*(e1-b))).*( 3*Rb^3 - 4*Rb.*(z3+b).^2 );
% 
% 
%    figure(1); hold on;
%    subplot(2,1,1), hold on;
%    plot(z4,abs(phi_4_Rp_r5_k2),'r','LineWidth',1);
%    plot(z2,abs(phi_2_Rp_r5_k2),'r','LineWidth',1);
%    plot(z1,abs(phi_1_Rp_r5_k2),'--');
%    subplot(2,1,2), hold on;
%    plot(z4,abs(d_phi_4_Rp_r5_k2)./Rp,'r','LineWidth',1);
%    plot(z2,abs(d_phi_2_Rp_r5_k2)./Rp,'r','LineWidth',1);
%    plot(z1,abs(d_phi_1_Rp_r5_k2)./Rp,'--');
% 
%    
%    figure(2); hold on;
%    subplot(2,1,1), hold on;
%    plot(z3,abs(phi_3_Rb_r5_k2),'r','LineWidth',1);
%    plot(z2,abs(phi_2_Rb_r5_k2),'--');
%    subplot(2,1,2), hold on;
%    plot(z3,abs(d_phi_3_Rb_r5_k2)./Rb,'r','LineWidth',1);
%    plot(z2,abs(d_phi_2_Rb_r5_k2)./Rb,'--');

    % *** ATTENTION -> ici on reprend une dépendance en temps positive
    % exp(-iwt) -> exp(iwt)
    out.fz_1 = fz_1';
    out.fz_2 = fz_2';
    
    out.fz_1_Haskind = fz_1_Haskind';
    out.fz_2_Haskind = fz_2_Haskind';

    out.fz_1_FK = 0;%fz_1_FK';
    out.fz_2_FK = 0;%fz_2_FK';
    
    out.fx_1 = -fx_1';
    out.fx_2 = fx_2';
    
    out.fx_1_Haskind = fx_1_Haskind';
    out.fx_2_Haskind = fx_2_Haskind';
    
    out.ty_1 = ty_1';
    out.ty_2 = ty_2';
    
    out.ty_1_Haskind = ty_1_Haskind';
    out.ty_2_Haskind = ty_2_Haskind';
    
    % *** F_r(w) = - Z_r(w) * u(w) avec Z_r = B(w) + iw M_a(w)
    out.z_r33_11 = imag(z_r33_11)*omega + 1i*omega*real(z_r33_11);
    out.z_r33_12 = imag(z_r33_12)*omega + 1i*omega*real(z_r33_12);
    out.z_r33_22 = imag(z_r33_22)*omega + 1i*omega*real(z_r33_22);
    out.z_r33_21 = imag(z_r33_21)*omega + 1i*omega*real(z_r33_21);
    
    out.z_r11_11 = imag(z_r11_11)*omega + 1i*omega*real(z_r11_11);
    out.z_r11_12 = imag(z_r11_12)*omega + 1i*omega*real(z_r11_12);
    out.z_r11_22 = imag(z_r11_22)*omega + 1i*omega*real(z_r11_22);
    out.z_r11_21 = imag(z_r11_21)*omega + 1i*omega*real(z_r11_21);
    
    out.z_r15_11 = imag(z_r15_11)*omega + 1i*omega*real(z_r15_11);
    out.z_r15_12 = imag(z_r15_12)*omega + 1i*omega*real(z_r15_12);
    out.z_r15_22 = imag(z_r15_22)*omega + 1i*omega*real(z_r15_22);
    out.z_r15_21 = imag(z_r15_21)*omega + 1i*omega*real(z_r15_21);
    
    out.z_r55_11 = imag(z_r55_11)*omega + 1i*omega*real(z_r55_11);
    out.z_r55_12 = imag(z_r55_12)*omega + 1i*omega*real(z_r55_12);
    out.z_r55_22 = imag(z_r55_22)*omega + 1i*omega*real(z_r55_22);
    out.z_r55_21 = imag(z_r55_21)*omega + 1i*omega*real(z_r55_21);
    
    out.z_r51_11 = imag(z_r51_11)*omega + 1i*omega*real(z_r51_11);
    out.z_r51_12 = imag(z_r51_12)*omega + 1i*omega*real(z_r51_12);
    out.z_r51_22 = imag(z_r51_22)*omega + 1i*omega*real(z_r51_22);
    out.z_r51_21 = imag(z_r51_21)*omega + 1i*omega*real(z_r51_21);
end