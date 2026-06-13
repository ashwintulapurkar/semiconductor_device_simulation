
%get equilibrium n,p and phi using non-iterative method
[X_array, phi_array1]=ode45(@(z,phi) fun_dphi_by_dz_deg(z,phi),X_array,phi_s_deg,opts);
phi_array1=real(phi_array1');
  
zz=fun_get_charge_density(Nc, Nv, Ec, Ev, Ed, Ea, ND, NA, EF_deg, V_thermal, gD, gA, phi_array1);

N11=max(zz.n); N22=max(zz.p); N33=sqrt(N11*N22);
N1=max(max(NA,ND),N33); N2=max(max(NA,ND),N33); N_ni=sqrt(N1*N2); % we are setting normalization factors. If you don't want it, comment it.

n=zz.n(2:end-1)/N1; p=zz.p(2:end-1)/N2; % p and n are normalized   
ND_plus=zz.Ndplus(2:end-1); NA_minus=zz.Naminus(2:end-1);

V=phi_array1(2:end-1)/VN; %V is normalized now
V_L=phi_s_deg/VN; V_R=phi_array1(end)/N2; %V_R=0;


n_L=zz.n(1)/N1; p_L=zz.p(1)/N2;
n_R=zz.n(end)/N1; p_R=zz.p(end)/N2;

ND_plus_L=zz.Ndplus(1); NA_minus_L=zz.Naminus(1);
ND_plus_R=zz.Ndplus(end); NA_minus_R=zz.Naminus(end);

%end get equilibrium n,p and phi using non-iterative method

%return

% if the semiconductor width is small, you can use the following method to
% get equilibrium n,p and V

phi_array=phi_array1(2:end-1);  %initial pot from non-iterative method

t0=epsilon0*Ks/delta_x^2;

J1=(-2*t0*diag(ones(1,Np)))+(t0*diag(ones(1,Np-1),1))+(t0*diag(ones(1,Np-1),-1));

no_iterations=0; 
change=1;
doping_profile=(ND-NA)*ones(1,length(x_array));
Vbdy=[phi_s_deg,zeros(1,Np-1)]*t0;
while change>1e-8
    
zz=fun_get_charge_density(Nc, Nv, Ec, Ev, Ed, Ea, ND, NA, EF_deg, V_thermal, gD, gA, phi_array);
rho_L=zz.rho; n_L=zz.n; p_L=zz.p; NDplus_L=zz.Ndplus; NAminus_L=zz.Naminus;
rho=zz.rho; n=zz.n; p=zz.p; NDplus=zz.Ndplus; NAminus=zz.Naminus;


drho_by_dphi=fun_get_drho_by_dphi(Nc, Nv, Ec, Ev, Ed, Ea, ND, NA, EF_deg, V_thermal, gD, gA, phi_array);

J2=q*diag(drho_by_dphi); %J2=J2*0;


RHS=J2*phi_array'-Vbdy'-q*rho';

phi_array_new=(J1+J2)\RHS; phi_array_new=phi_array_new';

 
change=(max((abs(phi_array_new-phi_array))));

no_iterations=no_iterations+1;

beta=1; phi_array=phi_array+beta*(phi_array_new-phi_array);

end

zz=fun_get_charge_density(Nc, Nv, Ec, Ev, Ed, Ea, ND, NA, EF_deg, V_thermal, gD, gA, phi_array);
n=zz.n/N1; p=zz.p/N2; % p and n are normalized   
ND_plus=zz.Ndplus; NA_minus=zz.Naminus;

V=phi_array/VN; %V is normalized now
V_L=phi_s_deg/VN; V_R=0; 

zz=fun_get_charge_density(Nc, Nv, Ec, Ev, Ed, Ea, ND, NA, EF_deg, V_thermal, gD, gA, phi_s_deg);
n_L=zz.n/N1; p_L=zz.p/N2;
ND_plus_L=zz.Ndplus; NA_minus_L=zz.Naminus;

zz=fun_get_charge_density(Nc, Nv, Ec, Ev, Ed, Ea, ND, NA, EF_deg, V_thermal, gD, gA, 0);
n_R=zz.n/N1; p_R=zz.p/N2;
ND_plus_R=zz.Ndplus; NA_minus_R=zz.Naminus;

% end get equilibrium n,p and phi using iterative method

% Now calculate the electric field at left end and find out what VG should be. Use this value of VG for plotting etc. 



