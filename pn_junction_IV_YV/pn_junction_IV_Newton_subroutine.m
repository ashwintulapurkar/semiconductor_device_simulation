
xn=sqrt(2*K1*K2*epsilon0*NA*abs(Vbi-VA)/(q*ND*(K1*NA+K2*ND)));
xp=sqrt(2*K1*K2*epsilon0*ND*abs(Vbi-VA)/(q*NA*(K1*NA+K2*ND)));

V_new=V;
V_L=0/VN;  V_R=(Vbi-VA)/VN; %V_L and V_R are normalized to VN

Vbdy=t_poi*[V_L,zeros(1,Np-2),V_R];  %eqn. 3.13 and 3.14
change=1; no_iterations=0;

while change>1e-8
V=V_new;
get_Jacobian;

n_p_V_new=(J1+J2)\(J2*[n,p,V]'+[Rn,Rp,-rho_tilda]'-bdy');    %eqn 3.17
n=n_p_V_new(1:Np); p=n_p_V_new(Np+1:2*Np); V_new=n_p_V_new(2*Np+1:end);

n=n'; p=p'; V_new=V_new'; %column to row


change=(max((abs(V_new-V))));

no_iterations=no_iterations+1;
end
V=V_new;

%calculate electric field dependent mobility.
E=-diff([V_L,V,V_R]*VN)/delta_x;
mu_n = fun_get_mobility(mu0_n, 2, v_sat, E);   Dn=mu_n*V_thermal;  %Einstein reln
mu_p = fun_get_mobility(mu0_p, 1, v_sat, E);  Dp=mu_p*V_thermal;  %Einstein reln
%end calculate electric field dependent mobility.

U_diff=diff([V_L,V,V_R])*(VN/V_thermal);  %U=potential/kT, eqn 3.7c
B_U_diff=B(U_diff);
B_minusU_diff=B_U_diff.*exp(U_diff);

%calculate current from converged solution
n_all=[n_L,n,n_R]; p_all=[p_L,p,p_R];
Jn=N1*q/(delta_x)*Dn.*(-n_all(1:end-1).*B_minusU_diff + n_all(2:end).*B_U_diff);
Jp=-N2*q/(delta_x)*Dp.*(-p_all(1:end-1).*B_U_diff + p_all(2:end).*B_minusU_diff);
J_total=Jn+Jp;
%end calculate current from converged solution

% get drift and diffusion current density
Jn_drift=N1*q/(delta_x)*Dn.*B_U_diff.*(n_all(2:end)+n_all(1:end-1).*exp(U_diff/2)).*(1-exp(U_diff/2)); %appendix
Jn_diff=N1*q/(delta_x)*Dn.*B_U_diff.*(n_all(2:end)-n_all(1:end-1)).*exp(U_diff/2);

Jp_drift=-N2*q/(delta_x)*Dp.*B_minusU_diff.*(p_all(2:end)+p_all(1:end-1).*exp(-U_diff/2)).*(1-exp(-U_diff/2));
Jp_diff=N2*q/(delta_x)*Dp.*B_minusU_diff.*(p_all(1:end-1)-p_all(2:end)).*exp(-U_diff/2);
% end get drift and diffusion current density

%index = find(x_array == 0);
xn0=sqrt(2*K1*K2*epsilon0*NA*abs(Vbi)/(q*ND*(K1*NA+K2*ND)));
xp0=sqrt(2*K1*K2*epsilon0*ND*abs(Vbi)/(q*NA*(K1*NA+K2*ND)));

index = find(x_array > 0);
%index = find(x_array > (xn0-xp0)/2);

Efn_minus_Efp=kT*log(n(index(1))*p(index(1))*N1*N2/ni^2);

