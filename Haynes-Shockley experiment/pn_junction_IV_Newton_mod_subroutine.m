

E_old=E;
V_new=V;
V_L=0/VN;  V_R=(-VA)/VN; %V_L and V_R are normalized to VN

nbdy_add=n*2*tau1/(delta_t);  %see eqn 6.8a
pbdy_add=p*2*tau2/(delta_t); %see eqn 6.8b
bdy_add=[nbdy_add,pbdy_add,zeros(1,Np)]+f;  %see eqn 6.8

Vbdy=t_poi*[V_L,zeros(1,Np-2),V_R];

change=1; no_iterations=0;

while change>1e-6
V=V_new;
get_Jacobian_mod;

%n_p_V_new=(J1+J2)\(J2*[n,p,V]'+[Rn,Rp,-rho_tilda]'-bdy');
n_p_V_new=(J1_dash+J2)\(J2*[n,p,V]'+[Rn,Rp,-rho_tilda]'-(bdy+bdy_add)');

n=n_p_V_new(1:Np); p=n_p_V_new(Np+1:2*Np); V_new=n_p_V_new(2*Np+1:end);

n=n'; p=p'; V_new=V_new';


change=(max((abs(V_new-V))));

no_iterations=no_iterations+1;
end
V=V_new;

%calculate electric field dependent mobility.
E=-diff([V_L,V,V_R]*VN)/delta_x;
mu_n = fun_get_mobility(mu0_n, 2, v_sat, E);   Dn=mu_n*V_thermal; 
mu_p = fun_get_mobility(mu0_p, 1, v_sat, E);  Dp=mu_p*V_thermal; 


U_diff=diff([V_L,V,V_R])*(VN/V_thermal); %U=potential/kT
B_U_diff=B(U_diff);
B_minusU_diff=B_U_diff.*exp(U_diff);

n_all=[n_L,n,n_R]; p_all=[p_L,p,p_R];
Jn=N1*q/(delta_x)*Dn.*(-n_all(1:end-1).*B_minusU_diff + n_all(2:end).*B_U_diff);
Jp=-N2*q/(delta_x)*Dp.*(-p_all(1:end-1).*B_U_diff + p_all(2:end).*B_minusU_diff);
JD=epsilon0*Ks*(E-E_old)/delta_t;
J_total=Jn+Jp+JD;

% get drift and diffusion current density
Jn_drift=N1*q/(delta_x)*Dn.*B_U_diff.*(n_all(2:end)+n_all(1:end-1).*exp(U_diff/2)).*(1-exp(U_diff/2));
Jn_diff=N1*q/(delta_x)*Dn.*B_U_diff.*(n_all(2:end)-n_all(1:end-1)).*exp(U_diff/2);

Jp_drift=-N2*q/(delta_x)*Dp.*B_minusU_diff.*(p_all(2:end)+p_all(1:end-1).*exp(-U_diff/2)).*(1-exp(-U_diff/2));
Jp_diff=N2*q/(delta_x)*Dp.*B_minusU_diff.*(p_all(1:end-1)-p_all(2:end)).*exp(-U_diff/2);

% J1=blkdiag(M1,M2,M3); J1_dash=blkdiag(M1_dash,M2_dash,M3);
% J2=[J2_11,J2_12,J2_13; J2_21,J2_22,J2_23; J2_31,J2_32,zeros(Np,Np)];
% Vbdy=t_poi*[V_L,zeros(1,Np-2),V_R];
% bdy=[nbdy,pbdy,Vbdy];
f=J1*[n,p,V]'-[Rn,Rp,-rho_tilda]'+bdy'; f=f';
