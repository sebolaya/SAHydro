function out = twoCylinder_Rp_inf_a_Rb_withColumn(WECstructure,omega,Lambda_i,options)


    g = 9.81;
    rho = 1025;    
    h = 150;
    
    Rb = WECstructure.Rbo;
    Rc = WECstructure.Rc;
    Rp = WECstructure.Rp;
    
    b = WECstructure.b;
    e1 = WECstructure.e1;
    e2 = WECstructure.e2;
    
    Ni = options.Truncate.Ni;
    Nn = options.Truncate.Nn;
    Nj = options.Truncate.Nj;
    Nl = options.Truncate.Nl;
    
    epsilon = 1e-7;
    
    Zc = options.Zc;

    n = (1:1:Nn)';
    j = (1:1:Nj)';
    l = (1:1:Nl)';
    
    beta_j = (pi/(h-e2)).*j;
    alpha_n = (pi/(e1-b)).*n;   
    gamma_l = (pi/(h-b)).*l;
    
    lambda_0 = -imag(Lambda_i(1));
    lambda_i = Lambda_i(2:Ni+1);
        
    if options.Haskind
        
        step = 1e-5;

        V_z1 = (-e2:step:-e1);
        V_z2 = (-b:step:0);
        V_z3 = (-e1:step:-b);

        V_r1 = (Rc:step:Rp);
        V_r2 = (Rp:step:Rb);
        V_r3 = (0:step:Rp);

        C1 = -2*pi*rho*g/cosh(lambda_0*h);
        C2 = -2*1i*pi*rho*g/cosh(lambda_0*h);
        
    end
     
    Lambda_i = [-1i*lambda_0 lambda_i];
    
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
    
    delta_00 = lambda_0*Rb*(besselh(1,lambda_0*Rb) / besselh(0,lambda_0*Rb));
    delta_0i = lambda_i.*Rb.*(besselk(1,lambda_i.*Rb) ./ besselk(0,lambda_i.*Rb));   

    delta_10 = 1 - lambda_0*Rb*(besselh(0,lambda_0*Rb) / besselh(1,lambda_0*Rb));
    delta_1i = 1 + lambda_i.*Rb.*( besselk(0,lambda_i.*Rb) ./ besselk(1,lambda_i.*Rb) );
    
    Gamma_0j = [gamma_00;gamma_0j];
    Gamma_1j = [gamma_10;gamma_1j];
    
    Delta_0i = [delta_00 delta_0i];
    Delta_1i = [delta_10 delta_1i];

    [S_0l_prim_Rb,S_0l_tild_prim_Rb] = Tp_func(0,Rb,Rb,Rp,gamma_l);
    [S_0l_prim_Rp,S_0l_tild_prim_Rp] = Tp_func(0,Rp,Rb,Rp,gamma_l);
    
    [S_1l_prim_Rb,S_1l_tild_prim_Rb] = Tp_func(1,Rb,Rb,Rp,gamma_l);
    [S_1l_prim_Rp,S_1l_tild_prim_Rp] = Tp_func(1,Rp,Rb,Rp,gamma_l);
    
    [S_0l_int,S_0l_tild_int] = Ti_func(0,Rp,Rb,gamma_l);
    [S_1l_int,S_1l_tild_int] = Ti_func(1,Rp,Rb,gamma_l);

    [T_0n_prim_Rp,T_0n_tild_prim_Rp] = Tp_func(0,Rp,Rp,Rc,alpha_n);
    [T_0n_prim_Rc,T_0n_tild_prim_Rc] = Tp_func(0,Rc,Rp,Rc,alpha_n);

    [T_1n_prim_Rp,T_1n_tild_prim_Rp] = Tp_func(1,Rp,Rp,Rc,alpha_n);
    [T_1n_prim_Rc,T_1n_tild_prim_Rc] = Tp_func(1,Rc,Rp,Rc,alpha_n);   
    
    [T_0n_int,T_0n_tild_int] = Ti_func(0,Rc,Rp,alpha_n);
    [T_1n_int,T_1n_tild_int] = Ti_func(1,Rc,Rp,alpha_n);
    

    I1 = N_lambda_i.^(-.5).*( sin(Lambda_i.*h)-sin(Lambda_i.*(h-b)) ) ./ Lambda_i;              % % Int_(-b)^(0){ Z_i^I }
    I2 = [(e2-e1);sqrt(2).*(sin(gamma_l.*(h-e1)) - sin(gamma_l.*(h-e2)))./gamma_l];             % % Int_(-e2)^(-e1){ Z_l^II }
    I3 = N_lambda_i.^(-.5).*f2(Lambda_i,-b,0,0,h);                                                 % % Int_(-b)^(0){ z*Z_i^I }
    I4 = [.5*(e1^2-e2^2);sqrt(2).*f2(gamma_l,-e2,-e1,0,h)];                                         % % Int_(-e2)^(-e1){ z*Z_l^II }
    I5_1 = [(e1-b)^3/3;sqrt(2).*f1(gamma_l,e1,b,e1,h)];                                         % % Int_{-e1}^{-b}{(z+e1)^2*Z_l^II}
    I5_2 = [(e1-b)^3/3;sqrt(2).*f1(gamma_l,e1,b,b,h)];                                          % % Int_{-e1}^{-b}{(z+e1)^2*Z_l^II}
    I6 = [(h-b)^3/3;sqrt(2).*f1(gamma_l,h,b,h,h)];                                              % % Int_{-h}^{-b}{(z+h)^2*Z_l^II}
    I7 = N_lambda_i.^(-.5).*f1(Lambda_i,h,b,h,h);                                               % % Int_(-h)^(-b){ (z+h)^2*Z_i^I }
    I8 = [(h-e2)^3/3;sqrt(2).*f1(gamma_l,h,e2,h,h)];                                            % % Int_{-h}^{-e2}{(z+h)^2*Z_l^II}
    I9 = [.5*(b^2-e1^2);sqrt(2).*f2(alpha_n,-e1,-b,0,e1)];                                          % % Int_(-e1)^(-b){ z*Z_n^III }
    I10_1 = [(e1-b)^3/3;sqrt(2).*f1(alpha_n,e1,b,e1,e1)];                                        % % Int_(-e1)^(-b){ (z+e1)^2*Z_n^III }
    I10_2 = [(e1-b)^3/3;sqrt(2).*f1(alpha_n,e1,b,b,e1)];                                        % % Int_(-e1)^(-b){ (z+b)^2*Z_n^III }
    
%%% ***********************************************************************
%%% ***         Définition des matrices intrinsèques à la structure
%%% ***********************************************************************
% *** m = 0 -> z direction
    % % Interface r = Rb
    d_r3_tau_11 = diag(Delta_0i) + ((h-b)/h).*(K_ltau'*(diag(S_0l_prim_Rb)*K_ltau));
    d_r3_tau_12 = ((h-b)/h).*(diag(S_0l_tild_prim_Rb)*K_ltau)';
    d_r3_tau_13 = zeros(Ni+1,Nn+1);

    % % Interface r = Rp
    d_r3_tau_21 = - diag(S_0l_prim_Rp)*K_ltau ;
    d_r3_tau_22 = ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_0j)*L_jtau))...
                 + ((e1-b)/(h-b)).*(M_ntau'*(diag(T_0n_prim_Rp)*M_ntau)) - diag(S_0l_tild_prim_Rp);
    d_r3_tau_23 = ((e1-b)/(h-b)) .* (diag(T_0n_tild_prim_Rp)*M_ntau)';
    
    % % Interface r = Rc
    d_r3_tau_31 = zeros(Nn+1,Ni+1);
    d_r3_tau_32 = diag(T_0n_prim_Rc)*M_ntau;
    d_r3_tau_33 = diag(T_0n_tild_prim_Rc);
    
    D_r3_tau = [d_r3_tau_11 d_r3_tau_12 d_r3_tau_13;...
                d_r3_tau_21 d_r3_tau_22 d_r3_tau_23;...
                d_r3_tau_31 d_r3_tau_32 d_r3_tau_33];

% *** m = 1 -> x and beta direction        
    d_r1_tau_11 = diag(Delta_1i) + ((h-b)/h).*(K_ltau'*(diag(S_1l_prim_Rb)*K_ltau));
    d_r1_tau_12 = ((h-b)/h).*(diag(S_1l_tild_prim_Rb)*K_ltau)';
    d_r1_tau_13 = zeros(Ni+1,Nn+1);

    % % Interface r = Rp
    d_r1_tau_21 = - diag(S_1l_prim_Rp)*K_ltau ;
    d_r1_tau_22 = ((h-e2)/(h-b)).*(L_jtau'*(diag(Gamma_1j)*L_jtau))...
                 + ((e1-b)/(h-b)).*(M_ntau'*(diag(T_1n_prim_Rp)*M_ntau)) - diag(S_1l_tild_prim_Rp);
    d_r1_tau_23 = ((e1-b)/(h-b)) .* (diag(T_1n_tild_prim_Rp)*M_ntau)';
    
    % % Interface r = Rc
    d_r1_tau_31 = zeros(Nn+1,Ni+1);
    d_r1_tau_32 = diag(T_1n_prim_Rc)*M_ntau;
    d_r1_tau_33 = diag(T_1n_tild_prim_Rc);

    D_r1_tau = [d_r1_tau_11 d_r1_tau_12 d_r1_tau_13;...
                d_r1_tau_21 d_r1_tau_22 d_r1_tau_23;...
                d_r1_tau_31 d_r1_tau_32 d_r1_tau_33];

%%% ***********************************************************************
%%% ***********************************************************************
%%%                     Problème de RADIATION
%%% ***********************************************************************
%%% ***********************************************************************
%%% ***********************************************************************

%%%                         HEAVE MODE
if strcmp(options.dof,'all') || strcmp(options.dof,'heave')

%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, column fixed
%%% *********************************************************************** 
   % % solution particulière pour r = R
    p_03 = ( (e1-b) / 12 ) * ( 2 - 3*(Rp^2 / (e1-b)^2) );
    p_n3 = sqrt(2).*f1(alpha_n,e1,b,e1,e1)./(2*(e1-b)^2);
    P_n3 = [p_03;p_n3];
    
    r1_03 = ( 2*((h-b)^3-(h-e1)^3) - 3*Rp^2*(e1-b) ) / (12*(h-b)*(e1-b)) ;
    r_n3 = sqrt(2).*f1(alpha_n,e1,b,h,e1)./(2*(h-b)*(e1-b));
    R_n3 = [r1_03;r_n3];
    
    r2_03 = ( (h-b) / 12 ) * ( 2 - 3*(Rb^2 / (h-b)^2) );
    r_l3 = sqrt(2).*f1(gamma_l,h,b,h,h)./(2*(h-b)^2);
    R_l3 = [r2_03;r_l3];

    r3_03 = ( (h-e2)^2 / (12*(h-b)) ) * ( 2 - 3*(Rp^2 / (h-e2)^2 ));
    r_j3 = sqrt(2).*f1(beta_j,h,e2,h,h)./(2*(h-b)*(h-e2));
    R_j3 = [r3_03;r_j3];

    % % Interface r = Rb
    h_r3_tau_1_k1 = ((h-b)/h) .* K_ltau'*(S_0l_prim_Rb.*R_l3) +  (.5*Rb^2/h).*K_ltau(1,:)';
    % % Interface r = Rp
    h_r3_tau_2_k1 = (.5*Rp^2/(h-b)).*(M_ntau(1,:)'-[1;zeros(Nl,1)]) - S_0l_prim_Rp.*R_l3...
                    + ((e1-b)/(h-b)).*(M_ntau'*(T_0n_prim_Rp.*(P_n3-R_n3)))...
                    - ((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*R_j3));
    % % Interface r = Rc
    h_r3_tau_3_k1 = T_0n_prim_Rc.*(P_n3-R_n3) + [(.5*Rc^2/(e1-b));zeros(Nn,1)];
    
    H_r3_tau_k1 = [h_r3_tau_1_k1;h_r3_tau_2_k1;h_r3_tau_3_k1];

    coeff = D_r3_tau\H_r3_tau_k1;

    A_r3_k1 = coeff(1:Ni+1);
    C_r3_k1 = coeff(Ni+2:Ni+Nl+2);
    E_r3_k1 = coeff(Ni+Nl+3:end);
    
    B_r3_k1 = K_ltau*A_r3_k1 - R_l3;
    D_r3_k1 = M_ntau*C_r3_k1 - (P_n3 - R_n3);
    F_r3_k1 = L_jtau*C_r3_k1 + R_j3;

    tmp1 = S_0l_int(1)*B_r3_k1(1) + sqrt(2)*sum((-1).^l.*S_0l_int(2:end).*B_r3_k1(2:end));
    tmp2 = S_0l_tild_int(1)*C_r3_k1(1) + sqrt(2)*sum((-1).^l.*S_0l_tild_int(2:end).*C_r3_k1(2:end));
	tmp3 = .5*(h-b)*(.5*(Rb^2 - Rp^2) - (1/8)*( (Rb^4-Rp^4)/(h-b)^2 ));
    tmp4 = T_0n_int(1)*D_r3_k1(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*D_r3_k1(2:end));
    tmp5 = T_0n_tild_int(1)*E_r3_k1(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*E_r3_k1(2:end));
    tmp6 = .5*(e1-b)*(.5*(Rp^2-Rc^2) - (1/8)*( (Rp^4-Rc^4)/(e1-b)^2 ));
    z_r33_11 = 2*pi*rho*( tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6);

    tmp1 = .5*Rp^2*F_r3_k1(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k1(2:end);
    tmp2 = T_0n_int(1)*D_r3_k1(1) + sqrt(2)*sum(T_0n_int(2:end).*D_r3_k1(2:end));
    tmp3 = T_0n_tild_int(1)*E_r3_k1(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*E_r3_k1(2:end));
    tmp4 = - (Rp^4-Rc^4)/(16*(e1-b));
    z_r33_12 = 2*pi*rho*( tmp1 - (tmp2 + tmp3 + tmp4) );
            
    
%%% ***********************************************************************    
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************     
    % % solution particulière pour r = Rc
    p_03 = ( (e1-b) / 12 ) * ( 2 - 3*( Rp^2/(e1-b)^2) );
    p_n3 = (sqrt(2)*(e1-b)) ./ (n.^2 .* pi^2);
    P_n3 = -[p_03;p_n3];
    
    q_03 = ( (h-e2) / 12 ) * ( 2 - 3*(Rp/(h-e2))^2);
    q_j3 = (sqrt(2)*(h-e2).*(-1).^j) ./ (j.^2 .* pi^2);
    Q_j3 = [q_03;q_j3];
    
    h_r3_tau_1_k2 = zeros(Ni+1,1);
    h_r3_tau_2_k2 = (.5*Rp^2/(h-b)).*(L_jtau(1,:)'-M_ntau(1,:)')...
                    + ((h-e2)/(h-b)).*(L_jtau'*(Gamma_0j.*Q_j3))...
                    + ((e1-b)/(h-b)).*(M_ntau'*(T_0n_prim_Rp.*P_n3));
    h_r3_tau_3_k2 = T_0n_prim_Rc.*P_n3 - [(.5*Rc^2/(e1-b));zeros(Nn,1)];
    
    H_r3_tau_k2 = [h_r3_tau_1_k2;h_r3_tau_2_k2;h_r3_tau_3_k2];
    
    % % Calcul des coefficients A_mi, C_mn & E_mn
    coeff = D_r3_tau\H_r3_tau_k2;

    A_r3_k2 = coeff(1:Ni+1);
    C_r3_k2 = coeff(Ni+2:Ni+Nl+2);
    E_r3_k2 = coeff(Ni+Nl+3:end);
   
    B_r3_k2 = K_ltau*A_r3_k2;
    D_r3_k2 = M_ntau*C_r3_k2 - P_n3;
    F_r3_k2 = L_jtau*C_r3_k2 - Q_j3;
    
    tmp1 = S_0l_int(1)*B_r3_k2(1) + sqrt(2)*sum((-1).^l.*S_0l_int(2:end).*B_r3_k2(2:end));
    tmp2 = S_0l_tild_int(1)*C_r3_k2(1) + sqrt(2)*sum((-1).^l.*S_0l_tild_int(2:end).*C_r3_k2(2:end));
    tmp3 = T_0n_int(1)*D_r3_k2(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*D_r3_k2(2:end));
    tmp4 = T_0n_tild_int(1)*E_r3_k2(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*E_r3_k2(2:end));
    tmp5 = (Rp^4-Rc^4)/(16*(e1-b));
    z_r33_21 = 2*pi*rho*( tmp1 + tmp2 + tmp3 + tmp4 + tmp5);

    tmp1 = .5*Rp^2*F_r3_k2(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*F_r3_k2(2:end);
    tmp2 = .5*(h-e2)*(.5*Rp^2 - (1/8)*( Rp^4/(h-e2)^2 ));
    tmp3 = T_0n_int(1)*D_r3_k2(1) + sqrt(2)*sum(T_0n_int(2:end).*D_r3_k2(2:end));
    tmp4 = T_0n_tild_int(1)*E_r3_k2(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*E_r3_k2(2:end));
    tmp5 = - .5*(e1-b)*(.5*(Rp^2-Rc^2) - (1/8)*( (Rp^4-Rc^4)/(e1-b)^2 ));
    z_r33_22 = 2*pi*rho*( tmp1 + tmp2 - (tmp3 + tmp4 + tmp5));
    

    if options.Haskind

        tmp01 = -cosh(lambda_0*(h-b))*sum( V_r1.*besselj(0,lambda_0.*V_r1) )*step;
        tmp02 = -cosh(lambda_0*(h-b))*sum( V_r2.*besselj(0,lambda_0.*V_r2) )*step;
        tmp0 = tmp01 + tmp02;

        % *** Vertical Axis

        tmp1_i = sum(cos(conj(Lambda_i')*(V_z2+h)) .* (ones(Ni+1,1)*cosh(lambda_0.*(V_z2+h))),2)*step;
        tmp1 = (-Rb*lambda_0*besselj(1,lambda_0*Rb))*(N_lambda_i.^(-.5)*(A_r3_k1.*tmp1_i));
        
        tmp21_l = sum(cos(gamma_l*(V_z1+h)) .* (ones(Nl,1)*cosh(lambda_0.*(V_z1+h))),2)*step;
        tmp21 = sqrt(2).*sum(C_r3_k1(2:end).*tmp21_l);
        tmp22 = C_r3_k1(1)*sum( cosh(lambda_0.*(V_z1+h)) )*step;
        tmp23 = (.5/(h-b))*sum( ((V_z1+h).^2 - .5*Rp^2) .* cosh(lambda_0.*(V_z1+h)) )*step;
        tmp2 = (-Rp*lambda_0*besselj(1,lambda_0*Rp))*(tmp21 + tmp22 + tmp23);
        
        tmp31_n = sum(cos(alpha_n*(V_z3+e1)) .* (ones(Nn,1)*cosh(lambda_0.*(V_z3+h))),2)*step;
        tmp31 = sqrt(2).*sum(E_r3_k1(2:end).*tmp31_n);
        tmp32 = E_r3_k1(1)*sum( cosh(lambda_0.*(V_z3+h)) )*step;
        tmp33 = (.5/(e1-b))*sum( ((V_z3+e1).^2 - .5*Rc^2) .* cosh(lambda_0.*(V_z3+h)) )*step;
        tmp3 = (-Rc*lambda_0*besselj(1,lambda_0*Rc))*(tmp31 + tmp32 + tmp33);
        
        % *** Horizontal Axis
        
        cst1 = besseli(0,alpha_n.*Rp).*besselk(0,alpha_n.*Rc) - besseli(0,alpha_n.*Rc).*besselk(0,alpha_n.*Rp);
        tmp411 = besselk(0,alpha_n.*Rc).*sum(besseli(0,alpha_n*V_r1).*(ones(Nn,1)*(besselj(0,lambda_0.*V_r1).*V_r1)),2).*step...
                    - besseli(0,alpha_n.*Rc).*sum(besselk(0,alpha_n*V_r1).*(ones(Nn,1)*(besselj(0,lambda_0.*V_r1).*V_r1)),2).*step;
        tmp412 = besseli(0,alpha_n.*Rp).*sum(besselk(0,alpha_n*V_r1).*(ones(Nn,1)*(besselj(0,lambda_0.*V_r1).*V_r1)),2).*step...
                    - besselk(0,alpha_n.*Rp).*sum(besseli(0,alpha_n*V_r1).*(ones(Nn,1)*(besselj(0,lambda_0.*V_r1).*V_r1)),2).*step;
        tmp41 = sqrt(2)*sum((-1).^n.*(D_r3_k1(2:end).*tmp411 + E_r3_k1(2:end).*tmp412)./cst1);

        tmp421 = sum( log(V_r1./Rc).*besselj(0,lambda_0.*V_r1).*V_r1 )*step/log(Rp/Rc);
        tmp422 = sum(log(Rp./V_r1).*besselj(0,lambda_0.*V_r1).*V_r1)*step/log(Rp/Rc);
        tmp42 = ( D_r3_k1(1)*tmp421 + E_r3_k1(1)*tmp422 );
        
        tmp43 = (.5/(e1-b))*sum( ((e1-b)^2 - .5.*V_r1.^2).*besselj(0,lambda_0.*V_r1).*V_r1 )*step;                  % solution particulière

        tmp4 = lambda_0*sinh(lambda_0*(h-b)) * (tmp41 + tmp42 + tmp43);
        
        tmp51 = sqrt(2)*sum( (D_r3_k1(2:end).*tmp411 + E_r3_k1(2:end).*tmp412)./cst1 );
        tmp53 = (.5/(e1-b))*sum( (-.5.*V_r1.^2).*besselj(0,lambda_0.*V_r1).*V_r1 )*step;                  % solution particulière
        tmp5 = lambda_0*sinh(lambda_0*(h-e1)) * (tmp51 + tmp42 + tmp53);
        
        cst2 = besseli(0,gamma_l.*Rb).*besselk(0,gamma_l.*Rp) - besseli(0,gamma_l.*Rp).*besselk(0,gamma_l.*Rb);
        tmp611 = besselk(0,gamma_l.*Rp).*sum(besseli(0,gamma_l*V_r2).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r2).*V_r2)),2).*step...
                    - besseli(0,gamma_l.*Rp).*sum(besselk(0,gamma_l*V_r2).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r2).*V_r2)),2).*step;
        tmp612 = besseli(0,gamma_l.*Rb).*sum(besselk(0,gamma_l*V_r2).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r2).*V_r2)),2).*step...
                    - besselk(0,gamma_l.*Rb).*sum(besseli(0,gamma_l*V_r2).*(ones(Nl,1)*(besselj(0,lambda_0.*V_r2).*V_r2)),2).*step;
        tmp61 = sqrt(2)*sum((-1).^l.*(B_r3_k1(2:end).*tmp611 + C_r3_k1(2:end).*tmp612)./cst2);
        
        tmp621 = sum( log(V_r2./Rp).*besselj(0,lambda_0.*V_r2).*V_r2 )*step/log(Rb/Rp);
        tmp622 = sum(log(Rb./V_r2).*besselj(0,lambda_0.*V_r2).*V_r2)*step/log(Rb/Rp);
        tmp62 = ( B_r3_k1(1)*tmp621 + C_r3_k1(1)*tmp622 );
        tmp63 = (.5/(h-b))*sum( ((h-b)^2 - .5.*V_r2.^2).*besselj(0,lambda_0.*V_r2).*V_r2 )*step;                  % solution particulière
        
        tmp6 = lambda_0*sinh(lambda_0*(h-b)) * (tmp61 + tmp62 + tmp63);

        tmp71 = F_r3_k1(1)*sum( besselj(0,lambda_0.*V_r3).*V_r3 )*step;
        tmp72_j = sum(besseli(0,beta_j*V_r3) .* (ones(Nj,1)*besselj(0,lambda_0.*V_r3)) .* (ones(Nj,1)*V_r3),2)*step;
        tmp72 = sqrt(2)*sum( (-1).^j.*tmp72_j.*F_r3_k1(2:end)./besseli(0,beta_j.*Rp) );
        tmp7 = lambda_0*sinh(lambda_0*(h-e2))*(tmp71 + tmp72);

        fz_1_Haskind = C1*( tmp0 - ( (tmp1 + tmp2 + tmp3) - tmp4 + tmp5 - tmp6 - tmp7) ) ;
        
        tmp01 = cosh(lambda_0*(h-e1))*sum( V_r1.*besselj(0,lambda_0.*V_r1) )*step;
        tmp02 = -cosh(lambda_0*(h-e2))*sum( V_r3.*besselj(0,lambda_0.*V_r3) )*step;
        tmp0 = tmp01 + tmp02;

        % *** Vertical Axis
        tmp1 = (-Rb*lambda_0*besselj(1,lambda_0*Rb))*(N_lambda_i.^(-.5)*(A_r3_k2.*tmp1_i));
        
        tmp21 = sqrt(2).*sum(C_r3_k2(2:end).*tmp21_l);
        tmp22 = C_r3_k2(1)*sum( cosh(lambda_0.*(V_z1+h)) )*step;
        tmp2 = (-Rp*lambda_0*besselj(1,lambda_0*Rp))*(tmp21 + tmp22);
        
        % *** column
        tmp31 = sqrt(2).*sum(E_r3_k2(2:end).*tmp31_n);
        tmp32 = E_r3_k2(1)*sum( cosh(lambda_0.*(V_z3+h)) )*step;
        tmp33 = -(.5/(e1-b))*sum( ((V_z3 + b).^2 -.5*Rc^2) .* cosh(lambda_0.*(V_z3+h)) )*step;        % solution particulière
        tmp3 = (-Rc*lambda_0*besselj(1,lambda_0*Rc))*(tmp31 + tmp32 + tmp33);
        
        % *** Horizontal Axis
        % *** buoy
        tmp41 = sqrt(2)*sum((-1).^n.*(D_r3_k2(2:end).*tmp411 + E_r3_k2(2:end).*tmp412)./cst1);
        tmp42 = ( D_r3_k2(1)*tmp421 + E_r3_k2(1)*tmp422 );
        tmp43 = -(.5/(e1-b))*sum( (-.5.*V_r1.^2).*besselj(0,lambda_0.*V_r1).*V_r1 )*step;                  % solution particulière
        tmp4 = lambda_0*sinh(lambda_0*(h-b)) * (tmp41 + tmp42 + tmp43);
        
        % *** plate
        tmp51 = sqrt(2)*sum( (D_r3_k2(2:end).*tmp411 + E_r3_k2(2:end).*tmp412)./cst1 );
        tmp53 = -(.5/(e1-b))*sum( ( (b-e1)^2 -.5.*V_r1.^2).*besselj(0,lambda_0.*V_r1).*V_r1 )*step;                  % solution particulière
        tmp5 = lambda_0*sinh(lambda_0*(h-e1)) * (tmp51 + tmp42 + tmp53);
        
        % buoy
        tmp61 = sqrt(2)*sum((-1).^l.*(B_r3_k2(2:end).*tmp611 + C_r3_k2(2:end).*tmp612)./cst2);
        tmp62 = ( B_r3_k2(1)*tmp621 + C_r3_k2(1)*tmp622 );
        tmp6 = lambda_0*sinh(lambda_0*(h-b)) * (tmp61 + tmp62);

        tmp71 = F_r3_k2(1)*sum( besselj(0,lambda_0.*V_r3).*V_r3 )*step;
        tmp72 = sqrt(2)*sum( (-1).^j.*tmp72_j.*F_r3_k2(2:end)./besseli(0,beta_j.*Rp) );
        tmp73 = (.5/(h-e2))*sum( ((h-e2)^2 - .5.*V_r3.^2) .* besselj(0,lambda_0.*V_r3) .*V_r3)*step;
        tmp7 = lambda_0*sinh(lambda_0*(h-e2))*(tmp71 + tmp72 + tmp73);

        fz_2_Haskind = C1*( tmp0 - ( (tmp1 + tmp2 + tmp3) - tmp4 + tmp5 - tmp6 - tmp7) ) ;
    else
        fz_1_Haskind = 0;
        fz_2_Haskind = 0;
    end % *** END OF HASKIND
else
    z_r33_11 = 0;
    z_r33_12 = 0;
    z_r33_21 = 0;
    z_r33_22 = 0;
    
    fz_1_Haskind = 0;
    fz_2_Haskind = 0;
end% *** END OF HEAVE


%%% ***********************************************************************

%%%                         SURGE MODE
if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
%%% ***********************************************************************
%%%     Cas k = 1       Buoy oscillating, column fixed
%%% ***********************************************************************

    h_r1_tau_1_k1 = - (Rb/h).*I1';          % Interface r = Rb
    h_r1_tau_2_k1 = zeros(Nl+1,1);          % Interface r = Rp
    h_r1_tau_3_k1 = zeros(Nn+1,1);          % Interface r = Rc
    H_r1_tau_k1 = [h_r1_tau_1_k1;h_r1_tau_2_k1;h_r1_tau_3_k1];
    
    coeff = D_r1_tau\H_r1_tau_k1;

    A_r1_k1 = coeff(1:Ni+1);
    C_r1_k1 = coeff(Ni+2:Ni+Nl+2);
    E_r1_k1 = coeff(Ni+Nl+3:end);
    
    B_r1_k1 = K_ltau*A_r1_k1;
    D_r1_k1 = M_ntau*C_r1_k1;
    F_r1_k1 = L_jtau*C_r1_k1;
    
    z_r11_11 = - pi*rho*Rb*(I1*A_r1_k1);
    z_r11_12 = - pi*rho*(Rp*I2'*C_r1_k1 + Rc*E_r1_k1(1)*(e1-b));
    
    tmp1 = S_1l_int(1)*B_r1_k1(1) + sqrt(2)*sum((-1).^l.*S_1l_int(2:end).*B_r1_k1(2:end));
    tmp2 = S_1l_tild_int(1)*C_r1_k1(1) + sqrt(2)*sum((-1).^l.*S_1l_tild_int(2:end).*C_r1_k1(2:end));
    tmp3 = T_1n_int(1)*D_r1_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r1_k1(2:end));
    tmp4 = T_1n_tild_int(1)*E_r1_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r1_k1(2:end));
    z_r15_11 = - pi*rho*(Rb*I3*A_r1_k1 + tmp1 + tmp2 + tmp3 + tmp4);
    
    tmp1 = .25*F_r1_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k1(2:end);
    tmp2 = T_1n_int(1)*D_r1_k1(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r1_k1(2:end));
    tmp3 = T_1n_tild_int(1)*E_r1_k1(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r1_k1(2:end));
    z_r15_12 = - pi*rho*( Rp*I4'*C_r1_k1 + tmp1 - (tmp2 + tmp3) + Rc*(I9'*E_r1_k1) );

%%% ***********************************************************************
%%%     Cas k = 2       Buoy fixed , column oscillating 
%%% ***********************************************************************


    h_r1_tau_1_k2 = zeros(Ni+1,1);          % Interface r = Rb
    h_r1_tau_2_k2 = -(Rp/(h-b)).*I2;        % Interface r = Rp
    h_r1_tau_3_k2 = [Rc;zeros(Nn,1)];       % Interface r = Rc
    H_r1_tau_k2 = [h_r1_tau_1_k2;h_r1_tau_2_k2;h_r1_tau_3_k2];
    
    coeff = D_r1_tau\H_r1_tau_k2;

    A_r1_k2 = coeff(1:Ni+1);
    C_r1_k2 = coeff(Ni+2:Ni+Nl+2);
    E_r1_k2 = coeff(Ni+Nl+3:end);
    
    B_r1_k2 = K_ltau*A_r1_k2;
    D_r1_k2 = M_ntau*C_r1_k2;
    F_r1_k2 = L_jtau*C_r1_k2;

    z_r11_21 = - pi*rho*Rb*(I1*A_r1_k2);
    z_r11_22 = - pi*rho*(Rp*I2'*C_r1_k2 + Rc*E_r1_k2(1)*(e1-b));
    
    tmp1 = S_1l_int(1)*B_r1_k2(1) + sqrt(2)*sum((-1).^l.*S_1l_int(2:end).*B_r1_k2(2:end));
    tmp2 = S_1l_tild_int(1)*C_r1_k2(1) + sqrt(2)*sum((-1).^l.*S_1l_tild_int(2:end).*C_r1_k2(2:end));
    tmp3 = T_1n_int(1)*D_r1_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r1_k2(2:end));
    tmp4 = T_1n_tild_int(1)*E_r1_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r1_k2(2:end));
    z_r15_21 = - pi*rho*(Rb*I3*A_r1_k2 + tmp1 + tmp2 + tmp3 + tmp4);
    
    tmp1 = .25*F_r1_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r1_k2(2:end);
    tmp2 = T_1n_int(1)*D_r1_k2(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r1_k2(2:end));
    tmp3 = T_1n_tild_int(1)*E_r1_k2(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r1_k2(2:end));
    z_r15_22 = - pi*rho*( Rp*I4'*C_r1_k2 + tmp1 - (tmp2+tmp3) + Rc*(I9'*E_r1_k2) );
    
    if options.Haskind
        tmp0 = Rb*besselj(1,lambda_0*Rb)*sum( cosh(lambda_0.*(V_z2+h)) )*step;

        % *** Vertical Axis
        % *** Buoy
        tmp1_i = sum(cos(conj(Lambda_i')*(V_z2+h)) .* (ones(Ni+1,1)*cosh(lambda_0.*(V_z2+h))),2)*step;
        tmp1 = (Rb*lambda_0*besselj(0,lambda_0*Rb) - besselj(1,lambda_0*Rb))*(N_lambda_i.^(-.5)*(A_r1_k1.*tmp1_i));
        
        % *** Plate
        tmp21_l = sum(cos(gamma_l*(V_z1+h)) .* (ones(Nl,1)*cosh(lambda_0.*(V_z1+h))),2)*step;
        tmp21 = sqrt(2).*sum(C_r1_k1(2:end).*tmp21_l);
        tmp22 = C_r1_k1(1)*sum( cosh(lambda_0.*(V_z1+h)) )*step;
        tmp2 = (Rp*lambda_0*besselj(0,lambda_0*Rp) - besselj(1,lambda_0*Rp))*(tmp21 + tmp22);
        
        % *** Column
        tmp31_n = sum(cos(alpha_n*(V_z3+e1)) .* (ones(Nn,1)*cosh(lambda_0.*(V_z3+h))),2)*step;
        tmp31 = sqrt(2).*sum(E_r1_k1(2:end).*tmp31_n);
        tmp32 = E_r1_k1(1)*sum( cosh(lambda_0.*(V_z3+h)) )*step;
        tmp3 = (Rc*lambda_0*besselj(0,lambda_0*Rc) - besselj(1,lambda_0*Rc))*(tmp31 + tmp32);
        
        % *** Horizontal Axis 
        % *** Buoy
        cst1 = besseli(1,alpha_n.*Rp).*besselk(1,alpha_n.*Rc) - besseli(1,alpha_n.*Rc).*besselk(1,alpha_n.*Rp);
        tmp411 = besselk(1,alpha_n.*Rc).*sum(besseli(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step...
                    - besseli(1,alpha_n.*Rc).*sum(besselk(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step;
        tmp412 = besseli(1,alpha_n.*Rp).*sum(besselk(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step...
                    - besselk(1,alpha_n.*Rp).*sum(besseli(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step;
        tmp41 = sqrt(2)*sum((-1).^n.*(D_r1_k1(2:end).*tmp411 + E_r1_k1(2:end).*tmp412)./cst1);

        tmp421 = sum( (V_r1./Rc - Rc./V_r1).*besselj(1,lambda_0.*V_r1).*V_r1 )*step/(Rp/Rc - Rc/Rp);
        tmp422 = sum( (Rp./V_r1 - V_r1./Rp).*besselj(1,lambda_0.*V_r1).*V_r1)*step/(Rp/Rc - Rc/Rp);
        tmp42 = ( D_r1_k1(1)*tmp421 + E_r1_k1(1)*tmp422 );
        tmp4 = lambda_0*sinh(lambda_0*(h-b)) * (tmp41 + tmp42);
        
        tmp51 = sqrt(2)*sum( (D_r1_k1(2:end).*tmp411 + E_r1_k1(2:end).*tmp412)./cst1 );
        tmp5 = lambda_0*sinh(lambda_0*(h-e1)) * (tmp51 + tmp42);
        
        
        cst2 = besseli(1,gamma_l.*Rb).*besselk(1,gamma_l.*Rp) - besseli(1,gamma_l.*Rp).*besselk(1,gamma_l.*Rb);
        tmp611 = besselk(1,gamma_l.*Rp).*sum(besseli(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step...
                    - besseli(1,gamma_l.*Rp).*sum(besselk(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step;
        tmp612 = besseli(1,gamma_l.*Rb).*sum(besselk(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step...
                    - besselk(1,gamma_l.*Rb).*sum(besseli(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step;
        tmp61 = sqrt(2)*sum((-1).^l.*(B_r1_k1(2:end).*tmp611 + C_r1_k1(2:end).*tmp612)./cst2);
        
        tmp621 = sum( (V_r2./Rp - Rp./V_r2).*besselj(1,lambda_0.*V_r2).*V_r2 )*step/(Rb/Rp - Rp/Rb);
        tmp622 = sum( (Rb./V_r2 - V_r2./Rb).*besselj(1,lambda_0.*V_r2).*V_r2)*step/(Rb/Rp - Rp/Rb);
        tmp62 = ( B_r1_k1(1)*tmp621 + C_r1_k1(1)*tmp622 );
        tmp6 = lambda_0*sinh(lambda_0*(h-b)) * (tmp61 + tmp62);

        tmp71 = (F_r1_k1(1)/Rp)*sum( V_r3 .* besselj(1,lambda_0.*V_r3) .* V_r3 )*step;
        tmp72_j = sum(besseli(1,beta_j*V_r3) .* (ones(Nj,1)*besselj(1,lambda_0.*V_r3)) .* (ones(Nj,1)*V_r3),2)*step;
        tmp72 = sqrt(2)*sum( (-1).^j.*tmp72_j.*F_r1_k1(2:end)./besseli(1,beta_j.*Rp) );
        tmp7 = lambda_0*sinh(lambda_0*(h-e2))*(tmp71 + tmp72);

        fx_1_Haskind = C2*( tmp0 - ( (tmp1 + tmp2 + tmp3) - tmp4 + tmp5 - tmp6 - tmp7) ) ;

        tmp01 = Rp*besselj(1,lambda_0*Rp)*sum( cosh(lambda_0.*(V_z1+h)) )*step;
        tmp02 = Rc*besselj(1,lambda_0*Rc)*sum( cosh(lambda_0.*(V_z3+h)) )*step;
        tmp0 = tmp01 + tmp02;

        % *** Vertical Axis
        % *** Buoy
        tmp1 = (Rb*lambda_0*besselj(0,lambda_0*Rb) - besselj(1,lambda_0*Rb))*(N_lambda_i.^(-.5)*(A_r1_k2.*tmp1_i));
        
        % *** Plate
        tmp21 = sqrt(2).*sum(C_r1_k2(2:end).*tmp21_l);
        tmp22 = C_r1_k2(1)*sum( cosh(lambda_0.*(V_z1+h)) )*step;
        tmp2 = (Rp*lambda_0*besselj(0,lambda_0*Rp) - besselj(1,lambda_0*Rp))*(tmp21 + tmp22);
        
        % *** Column
        tmp31 = sqrt(2).*sum(E_r1_k2(2:end).*tmp31_n);
        tmp32 = E_r1_k2(1)*sum( cosh(lambda_0.*(V_z3+h)) )*step;
        tmp3 = (Rc*lambda_0*besselj(0,lambda_0*Rc) - besselj(1,lambda_0*Rc))*(tmp31 + tmp32);
        
        % *** Horizontal Axis 
        % *** Buoy
        tmp41 = sqrt(2)*sum((-1).^n.*(D_r1_k2(2:end).*tmp411 + E_r1_k2(2:end).*tmp412)./cst1);
        tmp42 = ( D_r1_k2(1)*tmp421 + E_r1_k2(1)*tmp422 );
        tmp4 = lambda_0*sinh(lambda_0*(h-b)) * (tmp41 + tmp42);
        
        tmp51 = sqrt(2)*sum( (D_r1_k2(2:end).*tmp411 + E_r1_k2(2:end).*tmp412)./cst1 );
        tmp5 = lambda_0*sinh(lambda_0*(h-e1)) * (tmp51 + tmp42);

        tmp61 = sqrt(2)*sum((-1).^l.*(B_r1_k2(2:end).*tmp611 + C_r1_k2(2:end).*tmp612)./cst2);
        tmp62 = ( B_r1_k2(1)*tmp621 + C_r1_k2(1)*tmp622 );
        tmp6 = lambda_0*sinh(lambda_0*(h-b)) * (tmp61 + tmp62);

        tmp71 = (F_r1_k2(1)/Rp)*sum( V_r3 .* besselj(1,lambda_0.*V_r3) .* V_r3 )*step;
        tmp72 = sqrt(2)*sum( (-1).^j.*tmp72_j.*F_r1_k2(2:end)./besseli(1,beta_j.*Rp) );
        tmp7 = lambda_0*sinh(lambda_0*(h-e2))*(tmp71 + tmp72);

        fx_2_Haskind = C2*( tmp0 - ( (tmp1 + tmp2 + tmp3) - tmp4 + tmp5 - tmp6 - tmp7) ) ;

    else
        fx_1_Haskind = 0;
        fx_2_Haskind = 0;
    end% *** END OF HASKIND 
else
    z_r11_11 = 0;
    z_r11_12 = 0;
    z_r15_11 = 0;
    z_r15_12 = 0;
    
    z_r11_21 = 0;
    z_r11_22 = 0;
    z_r15_21 = 0;
    z_r15_22 = 0;

    fx_1_Haskind = 0;
    fx_2_Haskind = 0;

end% *** END OF HEAVE

%%% ***********************************************************************

%%%                         PITCH MODE

if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
%%% ***********************************************************************
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
    
    r_05 = (Rb^3/(24*(h-b))) * ( 3 - 4*((h-b)/Rb)^2 );
    r_l5 = -((sqrt(2)*Rb)/(2*(h-b)^2)).*f1(gamma_l,h,b,h,h);
    R_l5 = [r_05;r_l5];

    q_05 = (Rp^3/(24*(h-b))) * ( 3 - 4*((h-e2)/Rp)^2);
    q_j5 = -((sqrt(2)*Rp)/(2*(h-b)*(h-e2))).*f1(beta_j,h,e2,h,h);
    Q_j5 = [q_05;q_j5];

    h_r5_tau_1_k1 = ((h-b)/h).*K_ltau'*(S_1l_prim_Rb.*R_l5) - ((3*Rb^3)/(8*h)).*K_ltau(1,:)' + (Rb/(2*h*(h-b))).*I7' - (Rb/h).*I3'; 
    h_r5_tau_2_k1 = - S_1l_prim_Rp.*R_l5 + [3*Rp^3/(8*(h-b));zeros(Nl,1)] - (.5*Rp/(h-b)^2).*I6...
            + ((e1-b)/(h-b)).*(M_ntau'*(T_1n_prim_Rp.*(P_n5-S_n5))) - ((3*Rp^3)/(8*(h-b))).*M_ntau(1,:)' + (Rp/(2*(e1-b)*(h-b))).*I5_1...
            - ((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*Q_j5)); 
    h_r5_tau_3_k1 = (T_1n_prim_Rc.*(P_n5-S_n5)) - [(3*Rc^3)/(8*(e1-b));zeros(Nn,1)] + (Rc/(2*(e1-b)^2)).*I10_1;
    
    H_r5_tau_k1 = [h_r5_tau_1_k1;h_r5_tau_2_k1;h_r5_tau_3_k1];
    
    coeff = D_r1_tau\H_r5_tau_k1;

    A_r5_k1 = coeff(1:Ni+1);
    C_r5_k1 = coeff(Ni+2:Ni+Nl+2);
    E_r5_k1 = coeff(Ni+Nl+3:end);
    
    B_r5_k1 = K_ltau*A_r5_k1 - R_l5;
    D_r5_k1 = M_ntau*C_r5_k1 - (P_n5 - S_n5);
    F_r5_k1 = L_jtau*C_r5_k1 + Q_j5;
    
    tmp1 = S_1l_int(1)*B_r5_k1(1) + sqrt(2)*sum((-1).^l.*S_1l_int(2:end).*B_r5_k1(2:end));
    tmp2 = S_1l_tild_int(1)*C_r5_k1(1) + sqrt(2)*sum((-1).^l.*S_1l_tild_int(2:end).*C_r5_k1(2:end));
    tmp3 = ( (Rb^6-Rp^6)/6 - (Rb^4-Rp^4)*(h-b)^2 ) / (8*(h-b));
    tmp4 = T_1n_int(1)*D_r5_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r5_k1(2:end));
    tmp5 = T_1n_tild_int(1)*E_r5_k1(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r5_k1(2:end));
    tmp6 = ((Rp^6-Rc^6)/6 - (Rp^4-Rc^4)*(e1-b)^2) / (8*(e1-b));
    z_r55_11 = - pi*rho*(Rb*I3*A_r5_k1 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5 + tmp6);
    
    tmp1 = (.5*Rp^3*(e1^2-e2^2) - 4*Rp*( e2*(h-e2)^3 - e1*(h-e1)^3 - ((h-e1)^4-(h-e2)^4)/4 )/3 ) / (8*(h-b));
    tmp2 = .25*F_r5_k1(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k1(2:end);
    tmp3 = T_1n_int(1)*D_r5_k1(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r5_k1(2:end));
    tmp4 = T_1n_tild_int(1)*E_r5_k1(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r5_k1(2:end));
    tmp5 = (Rp^6-Rc^6)/(48*(e1-b));
    tmp6 = ( .125/(e1-b) )*(.5*Rc^3*(b^2-e1^2) - 4*Rc*( -(1/12)*(e1-b)^4 - (1/3)*b*(e1-b)^3 ));
    z_r55_12 = - pi*rho*( Rp*(I4'*C_r5_k1 + tmp1) + tmp2 - (tmp3 + tmp4 + tmp5) + Rc*(I9'*E_r5_k1 + tmp6) );
    
    z_r51_11 =  - pi*rho*Rb*(I1*A_r5_k1);
    
    tmp1 = (Rp^3*(e2-e1) - 4*Rp*((h-e1)^3-(h-e2)^3)/3)/(8*(h-b));
    tmp2 = Rc^3/8 - Rc*(e1-b)^2/6;
    z_r51_12 = - pi*rho*( Rp*(I2'*C_r5_k1 + tmp1) + Rc*(E_r5_k1(1)*(e1-b) + tmp2) );
    


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
    h_r5_tau_2_k2 = ((e1-b)/(h-b)).*(M_ntau'*(T_1n_prim_Rp.*P_n5)) + ((3*Rp^3)/(8*(h-b))).*M_ntau(1,:)' - (Rp/(2*(e1-b)*(h-b))).*I5_2...
               + ((h-e2)/(h-b)).*(L_jtau'*(Gamma_1j.*Q_j5)) - ((3*Rp^3)/(8*(h-b))).*L_jtau(1,:)' + (Rp/(2*(h-e2)*(h-b))).*I8 - (Rp/(h-b)).*I4;
    h_r5_tau_3_k2 = (T_1n_prim_Rc.*P_n5) + [(3*Rc^3)/(8*(e1-b));zeros(Nn,1)] - (Rc/(2*(e1-b)^2)).*I10_2  + (Rc/(e1-b)).*I9;

    H_r5_tau_k2 = [h_r5_tau_1_k2;h_r5_tau_2_k2;h_r5_tau_3_k2];
    
    coeff = D_r1_tau\H_r5_tau_k2;
    
    A_r5_k2 = coeff(1:Ni+1);
    C_r5_k2 = coeff(Ni+2:Ni+Nl+2);
    E_r5_k2 = coeff(Ni+Nl+3:end);
    
    B_r5_k2 = K_ltau*A_r5_k2;
    D_r5_k2 = M_ntau*C_r5_k2 - P_n5 ;
    F_r5_k2 = L_jtau*C_r5_k2 - Q_j5;
    
    
    tmp1 = S_1l_int(1)*B_r5_k2(1) + sqrt(2)*sum((-1).^l.*S_1l_int(2:end).*B_r5_k2(2:end));
    tmp2 = S_1l_tild_int(1)*C_r5_k2(1) + sqrt(2)*sum((-1).^l.*S_1l_tild_int(2:end).*C_r5_k2(2:end));
    tmp3 = T_1n_int(1)*D_r5_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*D_r5_k2(2:end));
    tmp4 = T_1n_tild_int(1)*E_r5_k2(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*E_r5_k2(2:end));
    tmp5 = -((Rp^6-Rc^6)/6) / (8*(e1-b));
    z_r55_21 = - pi*rho*(Rb*I3*A_r5_k2 + tmp1 + tmp2 + tmp3 + tmp4 + tmp5);
    
    tmp1 = .25*F_r5_k2(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*F_r5_k2(2:end);
    tmp2 = (Rp^6/6 - Rp^4*(h-e2)^2) / (8*(h-e2));
    tmp3 = T_1n_int(1)*D_r5_k2(1) + sqrt(2)*sum(T_1n_int(2:end).*D_r5_k2(2:end));
    tmp4 = T_1n_tild_int(1)*E_r5_k2(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*E_r5_k2(2:end));
    tmp5 = -( (Rp^6-Rc^6)/6 - (Rp^4-Rc^4)*(e1-b)^2) / (8*(e1-b));
    tmp6 = - ( .125/(e1-b) )*(.5*Rc^3*(b^2-e1^2) - 4*Rc*( (1/12)*(e1-b)^4 - (1/3)*e1*(e1-b)^3 ));
    z_r55_22 = - pi*rho*( Rp*I4'*C_r5_k2 + tmp1 + tmp2 - (tmp3 + tmp4 + tmp5) + Rc*(I9'*E_r5_k2 + tmp6) );

    z_r51_21 = - pi*rho*Rb*(I1*A_r5_k2);
    
    tmp1 = - Rc^3/8 + Rc*(e1-b)^2/6;
    z_r51_22 = - pi*rho*(Rp*I2'*C_r5_k2 + Rc*((e1-b)*E_r5_k2(1) + tmp1));
    
    % *** HASKIND (GREEN's Second theorem on Sb) *** %
    if options.Haskind
      
        tmp01 = Rb*besselj(1,lambda_0*Rb)*sum( (V_z2 - Zc(1)).*cosh(lambda_0.*(V_z2+h)) )*step;
        tmp02 = cosh(lambda_0*(h-b))*sum( V_r1.^2.*besselj(1,lambda_0.*V_r1) )*step;
        tmp03 = cosh(lambda_0*(h-b))*sum( V_r2.^2.*besselj(1,lambda_0.*V_r2) )*step;
        tmp0 = tmp01 + tmp02 + tmp03;

        % *** Vertical Axis
        % *** Buoy
        tmp1_i = sum(cos(conj(Lambda_i')*(V_z2+h)) .* (ones(Ni+1,1)*cosh(lambda_0.*(V_z2+h))),2)*step;
        tmp1 = (Rb*lambda_0*besselj(0,lambda_0*Rb) - besselj(1,lambda_0*Rb))*(N_lambda_i.^(-.5)*(A_r5_k1.*tmp1_i));
        
        % *** Plate
        tmp21_l = sum(cos(gamma_l*(V_z1+h)) .* (ones(Nl,1)*cosh(lambda_0.*(V_z1+h))),2)*step;
        tmp21 = sqrt(2).*sum(C_r5_k1(2:end).*tmp21_l);
        tmp22 = C_r5_k1(1)*sum( cosh(lambda_0.*(V_z1+h)) )*step;
        tmp23 = (.125/(h-b))*sum((Rp^3 - 4*Rp.*(V_z1+h).^2).*cosh(lambda_0.*(V_z1+h)))*step;        % solution particulière
        tmp2 = (Rp*lambda_0*besselj(0,lambda_0*Rp) - besselj(1,lambda_0*Rp))*(tmp21 + tmp22 + tmp23);
        
        % *** Column
        tmp31_n = sum(cos(alpha_n*(V_z3+e1)) .* (ones(Nn,1)*cosh(lambda_0.*(V_z3+h))),2)*step;
        tmp31 = sqrt(2).*sum(E_r5_k1(2:end).*tmp31_n);
        tmp32 = E_r5_k1(1)*sum( cosh(lambda_0.*(V_z3+h)) )*step;
        tmp33 = (.125/(e1-b))*sum((Rc^3 - 4*Rc.*(V_z3+e1).^2).*cosh(lambda_0.*(V_z3+h)))*step;
        tmp3 = (Rc*lambda_0*besselj(0,lambda_0*Rc) - besselj(1,lambda_0*Rc))*(tmp31 + tmp32 + tmp33);
        
        % *** Horizontal Axis 
        % *** Buoy
        cst1 = besseli(1,alpha_n.*Rp).*besselk(1,alpha_n.*Rc) - besseli(1,alpha_n.*Rc).*besselk(1,alpha_n.*Rp);
        tmp411 = besselk(1,alpha_n.*Rc).*sum(besseli(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step...
                    - besseli(1,alpha_n.*Rc).*sum(besselk(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step;
        tmp412 = besseli(1,alpha_n.*Rp).*sum(besselk(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step...
                    - besselk(1,alpha_n.*Rp).*sum(besseli(1,alpha_n*V_r1).*(ones(Nn,1)*(besselj(1,lambda_0.*V_r1).*V_r1)),2).*step;
        tmp41 = sqrt(2)*sum((-1).^n.*(D_r5_k1(2:end).*tmp411 + E_r5_k1(2:end).*tmp412)./cst1);
        tmp421 = sum( (V_r1./Rc - Rc./V_r1).*besselj(1,lambda_0.*V_r1).*V_r1 )*step/(Rp/Rc - Rc/Rp);
        tmp422 = sum( (Rp./V_r1 - V_r1./Rp).*besselj(1,lambda_0.*V_r1).*V_r1)*step/(Rp/Rc - Rc/Rp);
        tmp42 = ( D_r5_k1(1)*tmp421 + E_r5_k1(1)*tmp422 );
        tmp43 = (.125/(e1-b))*sum((V_r1.^3 - 4.*V_r1.*(e1-b)^2).*besselj(1,lambda_0.*V_r1).*V_r1)*step;
        tmp4 = lambda_0*sinh(lambda_0*(h-b)) * (tmp41 + tmp42 + tmp43);
        
        tmp51 = sqrt(2)*sum( (D_r5_k1(2:end).*tmp411 + E_r5_k1(2:end).*tmp412)./cst1 );
        tmp53 = (.125/(e1-b))*sum(V_r1.^3.*besselj(1,lambda_0.*V_r1).*V_r1)*step;
        tmp5 = lambda_0*sinh(lambda_0*(h-e1)) * (tmp51 + tmp42 + tmp53);
        
        cst2 = besseli(1,gamma_l.*Rb).*besselk(1,gamma_l.*Rp) - besseli(1,gamma_l.*Rp).*besselk(1,gamma_l.*Rb);
        tmp611 = besselk(1,gamma_l.*Rp).*sum(besseli(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step...
                    - besseli(1,gamma_l.*Rp).*sum(besselk(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step;
        tmp612 = besseli(1,gamma_l.*Rb).*sum(besselk(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step...
                    - besselk(1,gamma_l.*Rb).*sum(besseli(1,gamma_l*V_r2).*(ones(Nl,1)*(besselj(1,lambda_0.*V_r2).*V_r2)),2).*step;
        tmp61 = sqrt(2)*sum((-1).^l.*(B_r5_k1(2:end).*tmp611 + C_r5_k1(2:end).*tmp612)./cst2);
        
        tmp621 = sum( (V_r2./Rp - Rp./V_r2).*besselj(1,lambda_0.*V_r2).*V_r2 )*step/(Rb/Rp - Rp/Rb);
        tmp622 = sum( (Rb./V_r2 - V_r2./Rb).*besselj(1,lambda_0.*V_r2).*V_r2)*step/(Rb/Rp - Rp/Rb);
        tmp62 = ( B_r5_k1(1)*tmp621 + C_r5_k1(1)*tmp622 );
        tmp63 = (.125/(h-b))*sum((V_r2.^3 - 4.*V_r2.*(h-b)^2).*besselj(1,lambda_0.*V_r2).*V_r2)*step;
        tmp6 = lambda_0*sinh(lambda_0*(h-b)) * (tmp61 + tmp62 + tmp63);

        tmp71 = (F_r5_k1(1)/Rp)*sum( V_r3 .* besselj(1,lambda_0.*V_r3) .* V_r3 )*step;
        tmp72_j = sum(besseli(1,beta_j*V_r3) .* (ones(Nj,1)*besselj(1,lambda_0.*V_r3)) .* (ones(Nj,1)*V_r3),2)*step;
        tmp72 = sqrt(2)*sum( (-1).^j.*tmp72_j.*F_r5_k1(2:end)./besseli(1,beta_j.*Rp) );
        tmp7 = lambda_0*sinh(lambda_0*(h-e2))*(tmp71 + tmp72);

        ty_1_Haskind = C2*( tmp0 - ( (tmp1 + tmp2 + tmp3) - tmp4 + tmp5 - tmp6 - tmp7) ) ;

        tmp01 = Rp*besselj(1,lambda_0*Rp)*sum( (V_z1 - Zc(2)).*cosh(lambda_0.*(V_z1+h)) )*step;
        tmp02 = Rc*besselj(1,lambda_0*Rc)*sum( (V_z3 - Zc(2)).*cosh(lambda_0.*(V_z3+h)) )*step;        
        tmp03 = cosh(lambda_0*(h-e1))*sum( V_r1.^2.*besselj(1,lambda_0.*V_r1) )*step;
        tmp05 = cosh(lambda_0*(h-e2))*sum( V_r3.^2.*besselj(1,lambda_0.*V_r3) )*step;
        
        tmp0 = tmp01 + tmp02 - tmp03 + tmp05;

        % *** Vertical Axis
        % *** Buoy
        tmp1 = (Rb*lambda_0*besselj(0,lambda_0*Rb) - besselj(1,lambda_0*Rb))*(N_lambda_i.^(-.5)*(A_r5_k2.*tmp1_i));
        
        % *** Plate
        tmp21 = sqrt(2).*sum(C_r5_k2(2:end).*tmp21_l);
        tmp22 = C_r5_k2(1)*sum( cosh(lambda_0.*(V_z1+h)) )*step;
        tmp2 = (Rp*lambda_0*besselj(0,lambda_0*Rp) - besselj(1,lambda_0*Rp))*(tmp21 + tmp22);
        
        % *** Column
        tmp31 = sqrt(2).*sum(E_r5_k2(2:end).*tmp31_n);
        tmp32 = E_r5_k2(1)*sum( cosh(lambda_0.*(V_z3+h)) )*step;
        tmp33 = -(.125/(e1-b))*sum((Rc^3 - 4*Rc.*(V_z3+b).^2).*cosh(lambda_0.*(V_z3+h)))*step;
        tmp3 = (Rc*lambda_0*besselj(0,lambda_0*Rc) - besselj(1,lambda_0*Rc))*(tmp31 + tmp32 + tmp33);
        
        % *** Horizontal Axis 
        % *** Buoy
        tmp41 = sqrt(2)*sum((-1).^n.*(D_r5_k2(2:end).*tmp411 + E_r5_k2(2:end).*tmp412)./cst1);
        tmp42 = ( D_r5_k2(1)*tmp421 + E_r5_k2(1)*tmp422 );
        tmp43 = -(.125/(e1-b))*sum(V_r1.^3.*besselj(1,lambda_0.*V_r1).*V_r1)*step;
        tmp4 = lambda_0*sinh(lambda_0*(h-b)) * (tmp41 + tmp42 + tmp43);
        
        % *** plate
        tmp51 = sqrt(2)*sum( (D_r5_k2(2:end).*tmp411 + E_r5_k2(2:end).*tmp412)./cst1 );
        tmp53 = -(.125/(e1-b))*sum((V_r1.^3 - 4.*V_r1.*(b-e1)^2).*besselj(1,lambda_0.*V_r1).*V_r1)*step;
        tmp5 = lambda_0*sinh(lambda_0*(h-e1)) * (tmp51 + tmp42 + tmp53);
        
        tmp61 = sqrt(2)*sum((-1).^l.*(B_r5_k2(2:end).*tmp611 + C_r5_k2(2:end).*tmp612)./cst2);
        tmp62 = ( B_r5_k2(1)*tmp621 + C_r5_k2(1)*tmp622 );
        tmp6 = lambda_0*sinh(lambda_0*(h-b)) * (tmp61 + tmp62);

        tmp71 = (F_r5_k2(1)/Rp)*sum(V_r3 .* besselj(1,lambda_0.*V_r3).* V_r3 )*step;
        tmp72 = sqrt(2)*sum( (-1).^j.*tmp72_j.*F_r5_k2(2:end)./besseli(1,beta_j.*Rp) );
        tmp73 = (.125/(h-e2))*sum((V_r3.^3 - 4.*V_r3.*(h-e2)^2).*besselj(1,lambda_0.*V_r3).*V_r3)*step;
        tmp7 = lambda_0*sinh(lambda_0*(h-e2))*(tmp71 + tmp72 + tmp73);

        ty_2_Haskind = C2*( tmp0 - ( (tmp1 + tmp2 + tmp3) - tmp4 + tmp5 - tmp6 - tmp7) ) ;
        
    else

        ty_1_Haskind = 0;
        ty_2_Haskind = 0;

    end% *** END OF HASKIND
else
    z_r55_11 = 0;
    z_r55_12 = 0;
    z_r51_11 = 0;
    z_r51_12 = 0;
    
    z_r55_21 = 0;
    z_r55_22 = 0;
    z_r51_21 = 0;
    z_r51_22 = 0;
    
    ty_1_Haskind = 0;
    ty_2_Haskind = 0;
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
        Chi_0 = -lambda_0*Rb*besselj(1,lambda_0*Rb);
        
        h_d_0tau_1 = B.*( [Chi_0;zeros(Ni,1)] -  besselj(0,lambda_0*Rb)*((h-b)/h).*(K_ltau'*(S_0l_prim_Rb.*K_ltau(:,1))) );
        h_d_0tau_2 = B*besselj(0,lambda_0*Rb).*(S_0l_prim_Rp.*K_ltau(:,1));
        h_d_0tau_3 = zeros(Nn+1,1);

        H_d_0tau = [h_d_0tau_1;h_d_0tau_2;h_d_0tau_3];

        coeff = D_r3_tau\H_d_0tau;

        a_0i = coeff(1:Ni+1);
        c_0l = coeff(Ni+2:Ni+Nl+2);
        e_0n = coeff(Ni+Nl+3:end);

        b_0l = B.*besselj(0,lambda_0*Rb).*K_ltau(:,1) + K_ltau*a_0i;
        d_0n = M_ntau*c_0l;
        f_0j = L_jtau*c_0l;    

        tmp1 = S_0l_int(1)*b_0l(1) + sqrt(2)*sum((-1).^l.*S_0l_int(2:end).*b_0l(2:end));
        tmp2 = S_0l_tild_int(1)*c_0l(1) + sqrt(2)*sum((-1).^l.*S_0l_tild_int(2:end).*c_0l(2:end));
        tmp3 = T_0n_int(1)*d_0n(1) + sqrt(2)*sum((-1).^n.*T_0n_int(2:end).*d_0n(2:end));
        tmp4 = T_0n_tild_int(1)*e_0n(1) + sqrt(2)*sum((-1).^n.*T_0n_tild_int(2:end).*e_0n(2:end));
        fz_1 = 2*pi*a0*rho*(tmp1 + tmp2 + tmp3 + tmp4);

        tmp1 = .5*Rp^2*f_0j(1) + sqrt(2)*Rp*((-1).^j.*besseli(1,beta_j.*Rp)./(beta_j.*besseli(0,beta_j.*Rp)))'*f_0j(2:end);
        tmp2 = T_0n_int(1)*d_0n(1) + sqrt(2)*sum(T_0n_int(2:end).*d_0n(2:end));
        tmp3 = T_0n_tild_int(1)*e_0n(1) + sqrt(2)*sum(T_0n_tild_int(2:end).*e_0n(2:end));
        fz_2 = 2*pi*a0*rho*( tmp1 - (tmp2 + tmp3) );
    else
        fz_1 = 0;
        fz_2 = 0;
    end
    
    
    if ~strcmp(options.dof,'heave')
        Chi_1 = 2*1i*(lambda_0*Rb*besselj(0,lambda_0*Rb)-besselj(1,lambda_0*Rb));
        h_d_1tau_1 = B.*( [Chi_1;zeros(Ni,1)] -  (2*1i*besselj(1,lambda_0*Rb))*((h-b)/h).*(K_ltau'*(S_1l_prim_Rb.*K_ltau(:,1))) );
        h_d_1tau_2 = B*(2*1i*besselj(1,lambda_0*Rb)).*(S_1l_prim_Rp.*K_ltau(:,1));
        h_d_1tau_3 = zeros(Nn+1,1);

        H_d_1tau = [h_d_1tau_1;h_d_1tau_2;h_d_1tau_3];

        coeff = D_r1_tau\H_d_1tau;

        a_1i = coeff(1:Ni+1);
        c_1l = coeff(Ni+2:Ni+Nl+2);
        e_1n = coeff(Ni+Nl+3:end);

        b_1l = B*2*1i*besselj(1,lambda_0*Rb).*K_ltau(:,1) + K_ltau*a_1i;
        d_1n = M_ntau*c_1l;
        f_1j = L_jtau*c_1l;
    end

    if strcmp(options.dof,'all') || strcmp(options.dof,'surge')
        fx_1 = -pi*a0*rho*Rb*( B*2*1i*besselj(1,lambda_0*Rb)*I1(1) + I1*a_1i );
        fx_2 = -pi*a0*rho*( Rp*I2'*c_1l + Rc*e_1n(1)*(e1-b) );
    else
        fx_1 = 0;
        fx_2 = 0;
    end

    if strcmp(options.dof,'all') || strcmp(options.dof,'pitch')
        tmp1 = S_1l_int(1)*b_1l(1) + sqrt(2)*sum((-1).^l.*S_1l_int(2:end).*b_1l(2:end));
        tmp2 = S_1l_tild_int(1)*c_1l(1) + sqrt(2)*sum((-1).^l.*S_1l_tild_int(2:end).*c_1l(2:end));
        tmp3 = T_1n_int(1)*d_1n(1) + sqrt(2)*sum((-1).^n.*T_1n_int(2:end).*d_1n(2:end));
        tmp4 = T_1n_tild_int(1)*e_1n(1) + sqrt(2)*sum((-1).^n.*T_1n_tild_int(2:end).*e_1n(2:end));
        ty_1 = - pi*a0*rho*( Rb*(B*2*1i*besselj(1,lambda_0*Rb)*I3(1) + I3*a_1i) + tmp1 + tmp2 + tmp3 + tmp4 );

        tmp1 = .25*f_1j(1)*Rp^3 + sqrt(2)*Rp^2*(((-1).^j.*besseli(2,beta_j.*Rp))./(beta_j.*besseli(1,beta_j.*Rp)))'*f_1j(2:end);
        tmp2 = T_1n_int(1)*d_1n(1) + sqrt(2)*sum(T_1n_int(2:end).*d_1n(2:end));
        tmp3 = T_1n_tild_int(1)*e_1n(1) + sqrt(2)*sum(T_1n_tild_int(2:end).*e_1n(2:end));  
        ty_2 = - pi*a0*rho*( Rp*I4'*c_1l + tmp1 - (tmp2 + tmp3) + Rc*(I9'*e_1n) );
    else
        ty_1 = 0;
        ty_2 = 0;
    end

% % % % %%% ***********************************************************************
if ((options.matchingConditions == 1) && abs(omega - options.matchingFrequency)<=epsilon)

    step = 1e-3;

    V_z1 = -h:step:0;
    V_z2 = -h:step:-b;
    V_z3 = -e1:step:-b;
    V_z4 = -h:step:-e2;

    % *** Construction des solutions en z
    Z_i_1 = zeros(Ni+1,length(V_z1));
    for i = 1 : Ni+1
        Z_i_1(i,:) = N_lambda_i(i)^(-.5).*cos(Lambda_i(i).*(V_z1+h));
    end
    
    z_0_2 = ones(1,length(V_z2));
    z_l_2 = zeros(Nl,length(V_z2));
    for i = 1 : Nl
        z_l_2(i,:) = sqrt(2).*cos(gamma_l(i).*(V_z2+h));
    end
    Z_l_2 = [z_0_2;z_l_2];
    
    z_0_3 = ones(1,length(V_z3));
    z_n_3 = zeros(Nn,length(V_z3));
    for i = 1 : Nn
        z_n_3(i,:) = sqrt(2).*cos(alpha_n(i).*(V_z3+e1));
    end
    Z_n_3 = [z_0_3;z_n_3];

    z_0_4 = ones(1,length(V_z4));
    z_j_4 = zeros(Nj,length(V_z4));
    for i = 1 : Nj
        z_j_4(i,:) = sqrt(2).*cos(beta_j(i).*(V_z4+h));
    end
    Z_j_4 = [z_0_4;z_j_4];

    % *** *** Z direction - Scattering Problem
    phi_1_d_0 = B.*besselj(0,lambda_0*Rb).*Z_i_1(1,:) + conj(a_0i)'*Z_i_1;
    gradPhi_1_d_0 = B.*Chi_0.*Z_i_1(1,:)  - ((conj(a_0i)'.*Delta_0i)*Z_i_1);

    phi_2_d_Rb_0 = conj(b_0l)'*Z_l_2;
    gradPhi_2_d_Rb_0 = conj(S_0l_prim_Rb.*b_0l + S_0l_tild_prim_Rb.*c_0l)'*Z_l_2;
    
    phi_2_d_Rp_0 = conj(c_0l)'*Z_l_2;
    gradPhi_2_d_Rp_0 = conj(S_0l_prim_Rp.*b_0l + S_0l_tild_prim_Rp.*c_0l)'*Z_l_2;
 
    phi_3_d_0 = conj(d_0n)'*Z_n_3;
    gradPhi_3_d_0 = conj(T_0n_prim_Rp.*d_0n + T_0n_tild_prim_Rp.*e_0n)'*Z_n_3;
    
    phi_4_d_0 = conj(f_0j)'*Z_j_4;
    gradPhi_4_d_0 = conj(Gamma_0j.*f_0j)'*Z_j_4;
    
    % *** *** Heaving Mode
    % *** k = 1
    phi_1_r3_k1 = conj(A_r3_k1)'*Z_i_1;
    gradPhi_1_r3_k1 = - ((conj(A_r3_k1)'.*Delta_0i)*Z_i_1);
    
    phi_2_r3_Rb_k1 = conj(B_r3_k1)'*Z_l_2 + .5.*((V_z2+h).^2 - .5*Rb^2)./(h-b);
    gradPhi_2_r3_Rb_k1 = conj(S_0l_prim_Rb.*B_r3_k1 + S_0l_tild_prim_Rb.*C_r3_k1)'*Z_l_2 - .5*Rb^2/(h-b);
    
    phi_2_r3_Rp_k1 = conj(C_r3_k1)'*Z_l_2 + .5.*((V_z2+h).^2 - .5*Rp^2)./(h-b);
    gradPhi_2_r3_Rp_k1 = conj(S_0l_prim_Rp.*B_r3_k1 + S_0l_tild_prim_Rp.*C_r3_k1)'*Z_l_2 - .5*Rp^2/(h-b);

    phi_3_r3_k1 = conj(D_r3_k1)'*Z_n_3 + .5.*((V_z3+e1).^2 - .5*Rp^2)./(e1-b);
    gradPhi_3_r3_k1 = conj(T_0n_prim_Rp.*D_r3_k1 + T_0n_tild_prim_Rp.*E_r3_k1)'*Z_n_3 - .5*Rp^2/(e1-b);
    
    phi_4_r3_k1 = conj(F_r3_k1)'*Z_j_4;
    gradPhi_4_r3_k1 = conj(Gamma_0j.*F_r3_k1)'*Z_j_4;

     % *** k = 2
    phi_1_r3_k2 = conj(A_r3_k2)'*Z_i_1;
    gradPhi_1_r3_k2 = - ((conj(A_r3_k2)'.*Delta_0i)*Z_i_1);

    phi_2_r3_Rb_k2 = conj(B_r3_k2)'*Z_l_2;
    gradPhi_2_r3_Rb_k2 = conj(S_0l_prim_Rb.*B_r3_k2 + S_0l_tild_prim_Rb.*C_r3_k2)'*Z_l_2;
    
    phi_2_r3_Rp_k2 = conj(C_r3_k2)'*Z_l_2;
    gradPhi_2_r3_Rp_k2 = conj(S_0l_prim_Rp.*B_r3_k2 + S_0l_tild_prim_Rp.*C_r3_k2)'*Z_l_2;

    phi_3_r3_k2 = conj(D_r3_k2)'*Z_n_3 - .5.*((V_z3+b).^2 - .5*Rp^2)./(e1-b);
    gradPhi_3_r3_k2 = conj(T_0n_prim_Rp.*D_r3_k2 + T_0n_tild_prim_Rp.*E_r3_k2)'*Z_n_3 + .5*Rp^2/(e1-b);
    
    phi_4_r3_k2 = conj(F_r3_k2)'*Z_j_4 + .5.*((V_z4+h).^2 - .5*Rp^2)./(h-e2);
    gradPhi_4_r3_k2 = conj(Gamma_0j.*F_r3_k2)'*Z_j_4 - .5*Rp^2/(h-e2);
 
%        figure(1); hold on;
%        subplot(2,1,1), hold on;
%        plot(V_z2,abs(phi_2_d_Rb_0),'r','LineWidth',1);
%        plot(V_z1,abs(phi_1_d_0),'--');
%        subplot(2,1,2), hold on;
%        plot(V_z2,abs(gradPhi_2_d_Rb_0)./Rb,'r','LineWidth',1);
%        plot(V_z1,abs(gradPhi_1_d_0)./Rb,'--');
%        
%        figure(2); hold on;
%        subplot(2,1,1), hold on;
%        plot(V_z4,abs(phi_4_d_0),'r','LineWidth',1);
%        plot(V_z3,abs(phi_3_d_0),'r','LineWidth',1);
%        plot(V_z2,abs(phi_2_d_Rp_0),'--');
%        subplot(2,1,2), hold on;
%        plot(V_z4,abs(gradPhi_4_d_0)./Rp,'r','LineWidth',1);
%        plot(V_z3,abs(gradPhi_3_d_0)./Rp,'r','LineWidth',1);
%        plot(V_z2,abs(gradPhi_2_d_Rp_0)./Rp,'--');
% 
%        figure(3); hold on;
%        subplot(2,1,1), hold on;
%        plot(V_z2,abs(phi_2_r3_Rb_k1),'r','LineWidth',1);
%        plot(V_z1,abs(phi_1_r3_k1),'--');
%        subplot(2,1,2), hold on;
%        plot(V_z2,abs(gradPhi_2_r3_Rb_k1)./Rb,'r','LineWidth',1);
%        plot(V_z1,abs(gradPhi_1_r3_k1)./Rb,'--');
%        
%        figure(4); hold on;
%        subplot(2,1,1), hold on;
%        plot(V_z4,abs(phi_4_r3_k1),'r','LineWidth',1);
%        plot(V_z3,abs(phi_3_r3_k1),'r','LineWidth',1);
%        plot(V_z2,abs(phi_2_r3_Rp_k1),'--');
%        subplot(2,1,2), hold on;
%        plot(V_z4,abs(gradPhi_4_r3_k1)./Rp,'r','LineWidth',1);
%        plot(V_z3,abs(gradPhi_3_r3_k1)./Rp,'r','LineWidth',1);
%        plot(V_z2,abs(gradPhi_2_r3_Rp_k1)./Rp,'--');
% 
%        figure(5); hold on;
%        subplot(2,1,1), hold on;
%        plot(V_z2,abs(phi_2_r3_Rb_k2),'r','LineWidth',1);
%        plot(V_z1,abs(phi_1_r3_k2),'--');
%        subplot(2,1,2), hold on;
%        plot(V_z2,abs(gradPhi_2_r3_Rb_k2)./Rb,'r','LineWidth',1);
%        plot(V_z1,abs(gradPhi_1_r3_k2)./Rb,'--');
%        
%        figure(6); hold on;
%        subplot(2,1,1), hold on;
%        plot(V_z4,abs(phi_4_r3_k2),'r','LineWidth',1);
%        plot(V_z3,abs(phi_3_r3_k2),'r','LineWidth',1);
%        plot(V_z2,abs(phi_2_r3_Rp_k2),'--');
%        subplot(2,1,2), hold on;
%        plot(V_z4,abs(gradPhi_4_r3_k2)./Rp,'r','LineWidth',1);
%        plot(V_z3,abs(gradPhi_3_r3_k2)./Rp,'r','LineWidth',1);
%        plot(V_z2,abs(gradPhi_2_r3_Rp_k2)./Rp,'--');

    REP = 'matchingConditions/Rp_inf_Rb/';
%     REP = 'matchingConditions/Chau2010plate/';
    
    if ~isdir([options.DIR REP])
        mkdir(options.DIR,REP);
    end
    
    text1 = '# *** Matching conditions in heaving mode for Region %d *** #\n';
    text2 = '# For omega = %.2f, h = %.2f, Rb = %.2f, Rp = %.2f, Rc = %.2f, b = %.2f, e1 = %.2f, e2 = %.2f\n';
    text3 = 'z phi_k1 phi_k2 gradPhi_k1 gradPhi_k2\n\n';
    
    f_id = fopen([options.DIR REP 'Rp_' num2str(Rp/Rb) 'Rb_Heave_R1' '.dat'],'w');
    % print a title, followed by a blank line
    fprintf(f_id,text1,1);
    fprintf(f_id,text2,omega,h,Rb,Rp,Rc,b,e1,e2);
    fprintf(f_id,text3);
    for i = 1 : length(V_z1)
        fprintf(f_id,'%.3f %e %e %e %e\n',...
            V_z1(i),abs(phi_1_r3_k1(i)),abs(phi_1_r3_k2(i)),...
            abs(gradPhi_1_r3_k1(i)),abs(gradPhi_1_r3_k2(i)));
    end
    fclose(f_id);
    
    f_id = fopen([options.DIR REP 'Rp_' num2str(Rp/Rb) 'Rb_Heave_R2' '.dat'],'w');
    % print a title, followed by a blank line
    fprintf(f_id,text1,2);
    fprintf(f_id,text2,omega,h,Rb,Rp,Rc,b,e1,e2);
    fprintf(f_id,'z phi_k1_Rb phi_k2_Rb gradPhi_k1_Rb gradPhi_k2_Rb phi_k1_Rp phi_k2_Rp gradPhi_k1_Rp gradPhi_k2_Rp\n\n');
    for i = 1 : length(V_z2)
        fprintf(f_id,'%.3f %e %e %e %e %e %e %e %e\n',...
            V_z2(i),abs(phi_2_r3_Rb_k1(i)),abs(phi_2_r3_Rb_k2(i)),...
            abs(gradPhi_2_r3_Rb_k1(i)),abs(gradPhi_2_r3_Rb_k2(i)),...
            abs(phi_2_r3_Rp_k1(i)),abs(phi_2_r3_Rp_k2(i)),...
            abs(gradPhi_2_r3_Rp_k1(i)),abs(gradPhi_2_r3_Rp_k2(i)));
    end
    fclose(f_id);
    
    f_id = fopen([options.DIR REP 'Rp_' num2str(Rp/Rb) 'Rb_Heave_R3' '.dat'],'w');
    % print a title, followed by a blank line
    fprintf(f_id,text1,3);
    fprintf(f_id,text2,omega,h,Rb,Rp,Rc,b,e1,e2);
    fprintf(f_id,text3);
    for i = 1 : length(V_z3)
        fprintf(f_id,'%.3f %e %e %e %e\n',...
            V_z3(i),abs(phi_3_r3_k1(i)),abs(phi_3_r3_k2(i)),...
            abs(gradPhi_3_r3_k1(i)),abs(gradPhi_3_r3_k2(i)));
    end
    fclose(f_id);
    
    f_id = fopen([options.DIR REP 'Rp_' num2str(Rp/Rb) 'Rb_Heave_R4' '.dat'],'w');
    % print a title, followed by a blank line
    fprintf(f_id,text1,4);
    fprintf(f_id,text2,omega,h,Rb,Rp,Rc,b,e1,e2);
    fprintf(f_id,text3);
    for i = 1 : length(V_z4)
        fprintf(f_id,'%.3f %e %e %e %e\n',...
            V_z4(i),abs(phi_4_r3_k1(i)),abs(phi_4_r3_k2(i)),...
            abs(gradPhi_4_r3_k1(i)),abs(gradPhi_4_r3_k2(i)));
    end
    fclose(f_id);
       
end
   
%     out.fz_1 = fz_1;
%     out.fz_2 = fz_2;
%     
%     out.fz_1_Haskind = fz_1_Haskind;
%     out.fz_2_Haskind = fz_2_Haskind;
%     
%     out.fx_1 = fx_1;
%     out.fx_2 = fx_2;
%     
%     out.fx_1_Haskind = fx_1_Haskind;
%     out.fx_2_Haskind = fx_2_Haskind;
%     
%     out.ty_1 = ty_1;
%     out.ty_2 = ty_2;
%     
%     out.ty_1_Haskind = ty_1_Haskind;
%     out.ty_2_Haskind = ty_2_Haskind;
%     
%     out.z_r33_11 = z_r33_11;
%     out.z_r33_12 = z_r33_12;
%     out.z_r33_22 = z_r33_22;
%     out.z_r33_21 = z_r33_21;
%     
%     out.z_r11_11 = z_r11_11;
%     out.z_r11_12 = z_r11_12;
%     out.z_r11_22 = z_r11_22;
%     out.z_r11_21 = z_r11_21;
%     
%     out.z_r15_11 = z_r15_11;
%     out.z_r15_12 = z_r15_12;
%     out.z_r15_22 = z_r15_22;
%     out.z_r15_21 = z_r15_21;
%     
%     out.z_r55_11 = z_r55_11;
%     out.z_r55_12 = z_r55_12;
%     out.z_r55_22 = z_r55_22;
%     out.z_r55_21 = z_r55_21;
%     
%     out.z_r51_11 = z_r51_11;
%     out.z_r51_12 = z_r51_12;
%     out.z_r51_22 = z_r51_22;
%     out.z_r51_21 = z_r51_21;
%     
%%
    % *** ATTENTION -> ici on reprend une dépendance en temps positive
    % exp(-iwt) -> exp(iwt)
    out.fz_1 = fz_1';
    out.fz_2 = fz_2';
    
    out.fz_1_Haskind = fz_1_Haskind';
    out.fz_2_Haskind = fz_2_Haskind';
    
    out.fx_1 = fx_1';
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