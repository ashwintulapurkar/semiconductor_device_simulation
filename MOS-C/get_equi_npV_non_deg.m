%get equilibrium n,p and phi using non-iterative method
[X_array, phi_array1]=ode45(@(z,phi) fun_dphi_by_dz_non_deg(z,phi),X_array,phi_s_non_deg,opts);
phi_array1=real(phi_array1');
     

n=n0*exp(phi_array1(2:end-1)/kT)/N1; p=p0*exp(-phi_array1(2:end-1)/kT)/N2; % p and n are normalized   
V=phi_array1(2:end-1)/VN; %V is normalized now
V_L=phi_s_non_deg/VN; V_R=phi_array1(end)/N2; %V_R=0;
n_L=n0*exp(phi_s_non_deg/kT)/N1; p_L=p0*exp(-phi_s_non_deg/kT)/N2; 
n_R=n0*exp(phi_array1(end)/kT)/N1; p_R=p0*exp(-phi_array1(end)/kT)/N2; %n_R=n0/N1; p_R=p0/N2;

%end get equilibrium n,p and phi using non-iterative method
return
% if the semiconductor width is small, you can use the following method to
% get equilibrium n,p and V.
phi_array=phi_array1(2:end-1);  %initial pot from non-iterative method

t0=epsilon0*Ks/delta_x^2;

J1=(-2*t0*diag(ones(1,Np)))+(t0*diag(ones(1,Np-1),1))+(t0*diag(ones(1,Np-1),-1));

no_iterations=0; 
change=1;
doping_profile=(ND-NA)*ones(1,length(x_array));
Vbdy=[phi_s_non_deg,zeros(1,Np-1)]*t0;
while change>1e-8
 n=n0*exp(phi_array/kT); p=p0*exp(-phi_array/kT);
 
J2=(-q/kT)*diag((p+n)); %J2=J2*0;

RHS=J2*phi_array'-Vbdy'-q*(p-n+doping_profile)';
phi_array_new=(J1+J2)\RHS; phi_array_new=phi_array_new';

change=(max((abs(phi_array_new-phi_array))));

no_iterations=no_iterations+1;

beta=1; phi_array=phi_array+beta*(phi_array_new-phi_array);

end

n=n0*exp(phi_array/kT)/N1; p=p0*exp(-phi_array/kT)/N2; % p and n are normalized   
V=phi_array/VN; %V is normalized now
V_L=phi_s_non_deg/VN; V_R=0;
n_L=n0*exp(phi_s_non_deg/kT)/N1; p_L=p0*exp(-phi_s_non_deg/kT)/N2; 
n_R=n0/N1; p_R=p0/N2;
% end get equilibrium n,p and phi using iterative method

% Now calculate the electric field at left end and find out what VG should be. Use this value of VG for plotting etc. 

