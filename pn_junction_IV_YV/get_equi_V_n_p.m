
n0_R=Nc*exp((EF_R-Ec)/V_thermal); p0_R=Nv*exp((Ev-EF_R)/V_thermal);  %eqn 1.7
n0_L=Nc*exp((EF_L-Ec)/V_thermal); p0_L=Nv*exp((Ev-EF_L)/V_thermal); %eqn 1.7

n_R=n0_R/N1; p_R=p0_R/N2; %normlaized values  eqn. 3.7a
n_L=n0_L/N1; p_L=p0_L/N2; %normlaized values   eqn. 3.7a


V_L=0/VN;  V_R=(Vbi)/VN; %V_L and V_R are normalized to VN,   eqn. 3.7a


K1=Ks; K2=Ks; xn=sqrt(2*K1*K2*epsilon0*NA*abs(Vbi)/(q*ND*(K1*NA+K2*ND)));
xp=sqrt(2*K1*K2*epsilon0*ND*abs(Vbi)/(q*NA*(K1*NA+K2*ND)));

VA=0; 
[V_depletion_approx_array,E_depletion_approx_array, rho_depletion_approx_array,wL,wR]=fun_get_V_E_rho_from_depletion_approx(x_array,Ks,Ks, NA,ND, 0, Vbi);


doping_profile=-NA*(1-heaviside(x_array))+ND*heaviside(x_array);
doping_profile=doping_profile/N; %doping prifile normalized by N, eqn 3.11 a
LD=sqrt((Ks*epsilon0*VN/(N*q))); %eqn 3.11 a
t_poi=(LD/delta_x)^2;  %eqn 3.11 a
Vbdy=t_poi*[V_L,zeros(1,Np-2),V_R]; %eqn. 3.13 and 3.14

array1=-2*t_poi*ones(1,Np);
array2=t_poi*ones(1,Np-1);
M3=((diag(array1))+(diag(array2,1))+(diag(array2,-1)));  %eqn. 3.13

no_iterations=0; 

change=1;
V=V_depletion_approx_array/VN; %initial guess


while change>1e-6

n=n_L*exp(VN*V/V_thermal); p=p_L*exp(-VN*V/V_thermal); %n_L and p_L are scaled by N1 and N2--> n and p are scaled by N1 and N2; eqn. 1.6

J2=-(VN/V_thermal)*diag(((N2/N)*p+(N1/N)*n));  % second term in eqn 3.19, see eqn 3.22

rho_tilda=doping_profile-n*(N1/N)+p*(N2/N); %eqn. 3.11 a
RHS=(J2*V')'-rho_tilda-Vbdy; % RHS of eqn. 3.22

V_new=(M3+J2)\RHS'; % solution of eqn. 3.22
V_new=V_new'; %column to row

change=(max((abs(V_new-V))));
no_iterations=no_iterations+1;

V=V_new; 
end

n=n_L*exp(VN*V/V_thermal); p=p_L*exp(-VN*V/V_thermal); %eqn 1.6




