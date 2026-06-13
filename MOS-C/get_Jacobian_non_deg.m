
den_Rn=(tau_p/tau1)*((N1/N2)*n+sqrt(N1/N2)*ni/N_ni)+(tau_n/tau1)*(p+sqrt(N1/N2)*ni/N_ni);
% Rn=(n.*p-(ni/N_ni)^2)./den_Rn;
% dRn_by_dn=(p-(tau_p/tau1)*Rn*N1/N2)./den_Rn;
% dRn_by_dp= (n-(tau_n/tau1)*Rn*N1/N2)./den_Rn;

dRn_by_dn=(p)./den_Rn;  %equilibrium case
dRn_by_dp= (n)./den_Rn;  %equilibrium case

% Rp=Rn*(tau2/tau1)*(N1/N2);
dRp_by_dn=dRn_by_dn*(tau2/tau1)*(N1/N2);
dRp_by_dp=dRn_by_dp*(tau2/tau1)*(N1/N2);


%calculate electric field dependent mobility. Derivative of mobility is required for Jacobian
E=-diff([V_L,V,V_R]*VN)/delta_x; 
mu_n = fun_get_mobility(mu0_n, 2, v_sat, E);   Dn=mu_n*V_thermal; 
mu_p = fun_get_mobility(mu0_p, 1, v_sat, E);  Dp=mu_p*V_thermal; 

mu_n_dash = fun_get_mobility_der(mu0_n, 2, v_sat, E);
mu_p_dash = fun_get_mobility_der(mu0_p, 1, v_sat, E);

U_diff=diff([V_L,V,V_R])*(VN/V_thermal); 
B_U_diff=B(U_diff);
B_minusU_diff=B_U_diff.*exp(U_diff);
Bdash_U_diff=Bdash(U_diff);
Bdash_minusU_diff=Bdash(-U_diff);


%%% for electrons 

array1=Dn.*B_U_diff; array2=Dn.*B_minusU_diff;
array1=array1/D0_n; array2=array2/D0_n;
a=-(array1(1:end-1)+array2(2:end)); b=array1(2:end-1); c=array2(2:end-1);

M1=(diag(a)+(diag(b,1))+(diag(c,-1)))*tn;

J2_11=-diag(dRn_by_dn);

J2_12=-diag(dRn_by_dp);

% hn=Dn.*([n_L,n].*Bdash_minusU_diff + [n,n_R].*Bdash_U_diff); hn=hn/D0_n; hn=hn*VN/V_thermal;
% hn=hn-(mu_n_dash.*Jn./mu_n)*VN/(N1*q*D0_n);
% hn=hn*tn;

hn=-tn*Dn.*B_U_diff.*[n,n_R]/V_thermal*(VN/D0_n); %equilibrium case

J2_13=-diag(hn(1:end-1)+hn(2:end))+(diag(hn(2:end-1),1))+(diag(hn(2:end-1),-1));

%nbdy=[n_L*array2(1),zeros(1,Np-2),n_R*array1(end)]*tn;

%% for holes

array1=Dp.*B_minusU_diff; array2=Dp.*B_U_diff;
array1=array1/D0_p; array2=array2/D0_p;
a=-(array1(1:end-1)+array2(2:end)); b=array1(2:end-1); c=array2(2:end-1);

M2=(diag(a)+(diag(b,1))+(diag(c,-1)))*tp;
J2_22=-diag(dRp_by_dp);

J2_21=-diag(dRp_by_dn);

%hp=Dp.*([p_L,p].*Bdash_U_diff + [p,p_R].*Bdash_minusU_diff); hp=hp/D0_p; hp=hp*VN/V_thermal;
%hp=hp-(mu_p_dash.*Jp./mu_p)*VN/(N2*q*D0_p);
%hp=-tp*hp;

hp=tp*Dp.*B_minusU_diff.*[p,p_R]/V_thermal*(VN/D0_p); %equilibrium case

J2_23=-diag(hp(1:end-1)+hp(2:end))+(diag(hp(2:end-1),1))+(diag(hp(2:end-1),-1));

%pbdy=[p_L*array2(1),zeros(1,Np-2),p_R*array1(end)]*tp;

%%%% Poisson
J2_31=-(N1/N)*eye(Np); 
J2_32=(N2/N)*eye(Np);
J2_33=zeros(Np,Np);
rho_tilda=doping_profile-n*(N1/N)+p*(N2/N);
%%%%%%%%%%%%%%%%%

J1=blkdiag(M1,M2,M3);
J2=[J2_11,J2_12,J2_13; J2_21,J2_22,J2_23; J2_31,J2_32,J2_33];
% Vbdy=t_poi*[V_L,zeros(1,Np-2),V_R];
% bdy=[nbdy,pbdy,Vbdy];
% f=J1*[n,p,V]'-[Rn,Rp,-rho_tilda]'+bdy';
