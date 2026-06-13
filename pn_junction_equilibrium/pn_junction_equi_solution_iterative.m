
clear
hbar=1.06e-34;q=1.6e-19; epsilon0=8.85e-12; m0=9.1e-31; 
kB=1.38e-23; T=300; V_thermal=kB*T/q; kT=V_thermal;

input_for_pn_junction_Si;
%input_for_pn_junction_GaAs;

Vbi_deg=(EF_R_deg-EF_L_deg);
Vbi_non_deg=(EF_R_non_deg-EF_L_non_deg);

delta_x=20e-9;  
x_L=-1e-6; x_R=2e-6; 


x_array_L=[-delta_x:-delta_x:x_L]; x_array_R=[delta_x:delta_x:x_R];
no_L=length(x_array_L); no_R=length(x_array_R);
x_array=[flip(x_array_L),0,x_array_R]; Np=length(x_array);

doping_profile_L=(ND_L-NA_L)*ones(1,length(x_array_L));
doping_profile_R=(ND_R-NA_R)*ones(1,length(x_array_R));
doping_profile=[doping_profile_L,((ND_L-NA_L)+(ND_R-NA_R))/2,doping_profile_R];

x_array_for_E=x_array(1:end-1)+delta_x/2; %derivatives e.g. electric field are evaluated at mid-points
x_array_R=[0,x_array_R];
%get depletion apprximation solution (can be used as initial guess for itertaive solution)
[V_depl,E_depl, rho_depl, wL,wR] = fun_get_V_E_rho_from_depletion_approx(x_array_L,x_array_R,Ks_L,Ks_R, ND_L-NA_L, ND_R-NA_R, 0, Vbi_non_deg);
%end get depletion apprximation solution


t1=epsilon0*Ks_L/delta_x^2; t2=epsilon0*Ks_R/delta_x^2;  %eqn 2.6

array1=[-2*t1*ones(1,no_L),-t1-t2,-2*t2*ones(1,no_R)];
array2=[t1*ones(1,no_L),t2*ones(1,no_R)];
J1=((diag(array1))+(diag(array2,1))+(diag(array2,-1)));  %eqn 2.7a

Vbdy=[0,zeros(1,Np-2),Vbi_non_deg*t2];  %eqn 2.7d with phi_L=0, phi_R=Vbi.


%begin non-degenerate soln
%phi_array=1e-4*ones(1,Np);  %initial pot
phi_array=V_depl;  %initial pot
no_iterations=0; 
change=1;

while change>1e-6

n_L=(n0_L_non_deg)*exp(phi_array(1:no_L)/kT); p_L=(p0_L_non_deg)*exp(-phi_array(1:no_L)/kT);  % eqn 1.6
n_R=(n0_R_non_deg)*exp((phi_array(no_L+2:end)-Vbi_non_deg)/kT); p_R=(p0_R_non_deg)*exp((-phi_array(no_L+2:end)+Vbi_non_deg)/kT); % eqn 1.6
n=[n_L,(n_L(end)+n_R(1))/2,n_R]; p=[p_L,(p_L(end)+p_R(1))/2,p_R];

J2=(-q/kT)*diag((p+n)); % eqn 2.7b and 2.8
%J2=J2*0; %see discussion in the 2nd para below eqn 2.11

RHS=J2*phi_array'-Vbdy'-q*(p-n+doping_profile)'; %RHS of eqn 2.7c
phi_array_new=(J1+J2)\RHS; %solution of eqn 2.7c
phi_array_new=phi_array_new'; %change from column to row

change=(max((abs(phi_array_new-phi_array))));
if mod(no_iterations,100)==1, disp(change);  end
no_iterations=no_iterations+1;

phi_array=phi_array_new;

%beta=.1; phi_array=phi_array+beta*(phi_array_new-phi_array); % with J2=0, see discussion in the 2nd para below eqn 2.11

end

n_L=(n0_L_non_deg)*exp(phi_array(1:no_L)/kT); p_L=(p0_L_non_deg)*exp(-phi_array(1:no_L)/kT);
n_R=(n0_R_non_deg)*exp((phi_array(no_L+2:end)-Vbi_non_deg)/kT); p_R=(p0_R_non_deg)*exp((-phi_array(no_L+2:end)+Vbi_non_deg)/kT);
n=[n_L,(n_L(end)+n_R(1))/2,n_R]; p=[p_L,(p_L(end)+p_R(1))/2,p_R];

phi_array_non_deg=phi_array;
n_array_non_deg=n; p_array_non_deg=p;
rho_array_non_deg=(doping_profile-n+p);
E_array_non_deg=-diff(phi_array_non_deg)./diff(x_array);  %E=-dphi/dx
%end non-degenerate soln


% begin degenerate  soln
Vbdy=[0,zeros(1,Np-2),Vbi_deg*t2];
%phi_array=1e-4*ones(1,Np);  %initial pot
%phi_array=phi_array_non_deg;  %initial pot
phi_array=V_depl;  %initial pot
no_iterations=0; 
change=1;

while change>1e-6

zz=fun_get_charge_density(Nc_L, Nv_L, Ec_L, Ev_L, Ed_L, Ea_L, ND_L, NA_L, EF_L_deg, V_thermal, gD_L, gA_L, phi_array(1:no_L));
rho_L=zz.rho; n_L=zz.n; p_L=zz.p; NDplus_L=zz.Ndplus; NAminus_L=zz.Naminus;

zz=fun_get_charge_density(Nc_R, Nv_R, Ec_R, Ev_R, Ed_R, Ea_R, ND_R, NA_R, EF_R_deg, V_thermal, gD_R, gA_R, phi_array(no_L+2:end)-Vbi_deg);
rho_R=zz.rho; n_R=zz.n; p_R=zz.p;  NDplus_R=zz.Ndplus; NAminus_R=zz.Naminus;

n=[n_L,(n_L(end)+n_R(1))/2,n_R]; p=[p_L,(p_L(end)+p_R(1))/2,p_R];

rho=[rho_L,(rho_L(end)+rho_R(1))/2,rho_R];

drho_by_dphi_L=fun_get_drho_by_dphi(Nc_L, Nv_L, Ec_L, Ev_L, Ed_L, Ea_L, ND_L, NA_L, EF_L_deg, V_thermal, gD_L, gA_L, phi_array(1:no_L));
drho_by_dphi_R=fun_get_drho_by_dphi(Nc_R, Nv_R, Ec_R, Ev_R, Ed_R, Ea_R, ND_R, NA_R, EF_R_deg, V_thermal, gD_R, gA_R, phi_array(no_L+2:end)-Vbi_deg);
drho_by_dphi=[drho_by_dphi_L,(drho_by_dphi_L(end)+drho_by_dphi_R(1))/2,drho_by_dphi_R];

J2=q*diag(drho_by_dphi); %eqn 2.7b
%J2=J2*0; %see discussion in the 2nd para below eqn 2.11


RHS=J2*phi_array'-Vbdy'-q*rho'; %RHS of eqn 2.7c
phi_array_new=(J1+J2)\RHS;  %solution of eqn 2.7c
phi_array_new=phi_array_new'; %change from column to row


change=(max((abs(phi_array_new-phi_array))));
if mod(no_iterations,100)==1, disp(change);  end
no_iterations=no_iterations+1;

phi_array=phi_array_new;
%beta=.1; phi_array=phi_array+beta*(phi_array_new-phi_array); % with J2=0, see discussion in the 2nd para below eqn 2.11

end

zz=fun_get_charge_density(Nc_L, Nv_L, Ec_L, Ev_L, Ed_L, Ea_L, ND_L, NA_L, EF_L_deg, V_thermal, gD_L, gA_L, phi_array(1:no_L));
rho_array_L_deg=zz.rho; n_L=zz.n; p_L=zz.p; NDplus_L=zz.Ndplus; NAminus_L=zz.Naminus;

zz=fun_get_charge_density(Nc_R, Nv_R, Ec_R, Ev_R, Ed_R, Ea_R, ND_R, NA_R, EF_R_deg, V_thermal, gD_R, gA_R, phi_array(no_L+2:end)-Vbi_deg);
rho_array_R_deg=zz.rho; n_R=zz.n; p_R_deg=zz.p;  NDplus_R=zz.Ndplus; NAminus_R=zz.Naminus;

n=[n_L,(n_L(end)+n_R(1))/2,n_R]; p=[p_L,(p_L(end)+p_R(1))/2,p_R];

ionized_impurity_profile_L=NDplus_L-NAminus_L; ionized_impurity_profile_R=NDplus_R-NAminus_R;
ionized_impurity_profile=[ionized_impurity_profile_L,(ionized_impurity_profile_L(end)+ionized_impurity_profile_R(1))/2,ionized_impurity_profile_R];

phi_array_deg=phi_array;
n_array_deg=n; p_array_deg=p;
rho_array_deg=(ionized_impurity_profile-n+p);
E_array_deg=-diff(phi_array_deg)./diff(x_array);    %E=-dphi/dx
%end degenerate soln

plot_p_n_junction;
