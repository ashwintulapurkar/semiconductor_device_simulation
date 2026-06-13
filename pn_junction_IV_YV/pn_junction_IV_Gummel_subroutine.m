xn=sqrt(2*K1*K2*epsilon0*NA*abs(Vbi-VA)/(q*ND*(K1*NA+K2*ND)));
xp=sqrt(2*K1*K2*epsilon0*ND*abs(Vbi-VA)/(q*NA*(K1*NA+K2*ND)));

V_new=V;
V_L=0/VN;  V_R=(Vbi-VA)/VN; %V_L and V_R are normalized to VN
%V_L=(VA/2)/VN;  V_R=(Vbi-(VA/2))/VN; %V_L and V_R are normalized to VN


Vbdy=t_poi*[V_L,zeros(1,Np-2),V_R];  %eqn. 3.13 and 3.14
change=1; no_iterations=0;

while change>1e-6
V=V_new;

% R is defined by eqn 3.1c, the normalized versions Rn and Rp are defined in eqn 3.9
den_Rn=(tau_p/tau1)*((N1/N2)*n+sqrt(N1/N2)*ni/N_ni)+(tau_n/tau1)*(p+sqrt(N1/N2)*ni/N_ni);
Rn=(n.*p-(ni/N_ni)^2)./den_Rn;
dRn_by_dn=(p-(tau_p/tau1)*Rn*N1/N2)./den_Rn;
dRn_by_dp= (n-(tau_n/tau1)*Rn*N1/N2)./den_Rn;

% if we assume tau_n and tau_p->inf in the depletion region.-->Rn=0 and dRn/dn=0, dRn/dp=0
% uncomment following 3 lines if you want to switch off recmbination-generation in depletion region.
% Rn=Rn.*(heaviside(x_array-1*xn)+(1-heaviside(x_array+1*xp))); % use to turn off R-G in the depletion region
% dRn_by_dn=dRn_by_dn.*(heaviside(x_array-1*xn)+(1-heaviside(x_array+1*xp))); % use to turn off R-G in the depletion region
% dRn_by_dp=dRn_by_dp.*(heaviside(x_array-1*xn)+(1-heaviside(x_array+1*xp))); % use to turn off R-G in the depletion region


Rp=Rn*(tau2/tau1)*(N1/N2);
dRp_by_dn=dRn_by_dn*(tau2/tau1)*(N1/N2);
dRp_by_dp=dRn_by_dp*(tau2/tau1)*(N1/N2);


%calculate electric field dependent mobility
E=-diff([V_L,V,V_R]*VN)/delta_x; 
mu_n = fun_get_mobility(mu0_n, 2, v_sat, E);   Dn=mu_n*V_thermal; %Einstein reln
mu_p = fun_get_mobility(mu0_p, 1, v_sat, E);  Dp=mu_p*V_thermal;  %Einstein reln
%  E, mu_n, mu_p defined at mid points
%end calculate electric field dependent mobility

U_diff=diff([V_L,V,V_R])*(VN/V_thermal); %U=potential/kT, eqn 3.7c
B_U_diff=B(U_diff); %B=Brillouin function
B_minusU_diff=B_U_diff.*exp(U_diff); %B_minusU_diff=B(-U_diff);


%%% for electrons %electron continuity eqn
array1=Dn.*B_U_diff; array2=Dn.*B_minusU_diff;
array1=array1/D0_n; array2=array2/D0_n;
a=-(array1(1:end-1)+array2(2:end)); b=array1(2:end-1); c=array2(2:end-1);

M1=(diag(a)+(diag(b,1))+(diag(c,-1)))*tn; %eqn 3.13

nbdy=[n_L*array2(1),zeros(1,Np-2),n_R*array1(end)]*tn;   %eqn. 3.13 and 3.14

J2n=-diag(dRn_by_dn); %see eqn 3.20

n_new=(M1+J2n)\(J2n*n'+Rn'-nbdy');  %see eqn 3.20

%% for holes
array1=Dp.*B_minusU_diff; array2=Dp.*B_U_diff;
array1=array1/D0_p; array2=array2/D0_p;
a=-(array1(1:end-1)+array2(2:end)); b=array1(2:end-1); c=array2(2:end-1);

M2=(diag(a)+(diag(b,1))+(diag(c,-1)))*tp; %eqn 3.13

pbdy=[p_L*array2(1),zeros(1,Np-2),p_R*array1(end)]*tp;    %eqn. 3.13 and 3.14

J2p=-diag(dRn_by_dp);  %see eqn 3.21

p_new=(M2+J2p)\(J2p*p'+Rn'-pbdy'); %eqn 3.21

%update n and p now
n=n_new';
p=p_new';

%%%%%%%%%%%%%%%
%Poisson eqn

J2poi=-(VN/V_thermal)*diag(((N2/N)*p+(N1/N)*n)); % second term in eqn 3.19, see eqn 3.22

rho_tilda=doping_profile-n*(N1/N)+p*(N2/N); %eqn. 3.11 a

V_new=(M3+J2poi)\(J2poi*V'-rho_tilda'-Vbdy'); %M3 is defined in get_equi_U_n_p_norm,  % solution of eqn. 3.22
V_new=V_new'; %column to row

change=(max((abs(V_new-V))));

no_iterations=no_iterations+1;


end


%calculate electric field dependent mobility
E=-diff([V_L,V,V_R]*VN)/delta_x; 
mu_n = fun_get_mobility(mu0_n, 2, v_sat, E);   Dn=mu_n*V_thermal;  %Einstein reln
mu_p = fun_get_mobility(mu0_p, 1, v_sat, E);  Dp=mu_p*V_thermal;  %Einstein reln
%end calculate electric field dependent mobility

U_diff=diff([V_L,V,V_R])*(VN/V_thermal);  %U=potential/kT, eqn 3.7c
B_U_diff=B(U_diff);
B_minusU_diff=B_U_diff.*exp(U_diff);

%calculate current from converged solution
n_all=[n_L,n,n_R]; p_all=[p_L,p,p_R];
Jn=N1*q/(delta_x)*Dn.*(-n_all(1:end-1).*B_minusU_diff + n_all(2:end).*B_U_diff); %eqn 3.23a
Jp=-N2*q/(delta_x)*Dp.*(-p_all(1:end-1).*B_U_diff + p_all(2:end).*B_minusU_diff); %eqn 3.23b
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


