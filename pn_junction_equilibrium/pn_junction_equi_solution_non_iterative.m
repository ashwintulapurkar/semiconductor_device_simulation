clear all
clear
hbar=1.06e-34; kB=1.38e-23; q=1.6e-19; epsilon0=8.85e-12; m0=9.1e-31; T=300; V_thermal=kB*T/q; kT=V_thermal;

input_for_pn_junction_Si; 
%input_for_pn_junction_GaAs; %_HEMT;


Vbi_deg=(EF_R_deg-EF_L_deg); %built-in pot
Vbi_non_deg=(EF_R_non_deg-EF_L_non_deg); %built-in pot


delta_x=10e-9;  
x_L=-1e-6; x_R=1e-6; %widths of left and right semi-conductors

x_array_R=[0:delta_x:x_R]; x_array_L=[0:-delta_x:x_L];
x_array=[flip(x_array_L),x_array_R];

x_array_for_E=x_array; % x array where electric field is calculated.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define electric field as a function of pot for degenerate case for left hand side 
fun_f1_L_deg=@(phi)fermi(1.5,(Ev_L-EF_L_deg-phi)/kT)-fermi(1.5,(Ev_L-EF_L_deg)/kT); %eqn 1.9
fun_f2_L_deg=@(phi)fermi(1.5,(EF_L_deg-Ec_L+phi)/kT)-fermi(1.5,(EF_L_deg-Ec_L)/kT); %eqn 1.10
fun_f3_L_deg=@(phi) phi/kT+log( (1+gD_L*exp((EF_L_deg-Ed_L)/kT))./ (1+gD_L*exp((EF_L_deg-Ed_L+phi)/kT))); %eqn 1.13
fun_f4_L_deg=@(phi) phi/kT+log( (1+gA_L*exp((Ea_L-EF_L_deg-phi)/kT))./ (1+gA_L*exp((Ea_L-EF_L_deg)/kT))); %eqn 1.14
fun_integral_rho_L=@(phi) kT*( -Nv_L*fun_f1_L_deg(phi) -Nc_L*fun_f2_L_deg(phi) + ND_L*fun_f3_L_deg(phi) - NA_L*fun_f4_L_deg(phi)); %eqn 1.15

fun_E_L_deg=@(phi)  -sign(phi).*sqrt(-2/(epsilon0*Ks_L)*q*fun_integral_rho_L(phi)); %eqn 1.5
% end define electric field as a function of pot for degenerate case for left hand side 

% define electric field as a function of pot for degenerate case for right hand side 
fun_f1_R_deg=@(phi)fermi(1.5,(Ev_R-EF_R_deg-phi)/kT)-fermi(1.5,(Ev_R-EF_R_deg)/kT); %eqn 1.9
fun_f2_R_deg=@(phi)fermi(1.5,(EF_R_deg-Ec_R+phi)/kT)-fermi(1.5,(EF_R_deg-Ec_R)/kT); %eqn 1.10
fun_f3_R_deg=@(phi) phi/kT+log( (1+gD_R*exp((EF_R_deg-Ed_R)/kT))./ (1+gD_R*exp((EF_R_deg-Ed_R+phi)/kT))); %eqn 1.13
fun_f4_R_deg=@(phi) phi/kT+log( (1+gA_R*exp((Ea_R-EF_R_deg-phi)/kT))./ (1+gA_R*exp((Ea_R-EF_R_deg)/kT))); %eqn 1.14
fun_integral_rho_R=@(phi) kT*( -Nv_R*fun_f1_R_deg(phi) -Nc_R*fun_f2_R_deg(phi) + ND_R*fun_f3_R_deg(phi) - NA_R*fun_f4_R_deg(phi)); %eqn 1.15

fun_E_R_deg=@(phi)  sign(phi).*sqrt(-2/(epsilon0*Ks_R)*q*fun_integral_rho_R(phi)); %eqn 1.5
% end define electric field as a function of pot for degenerate case for right hand side 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% define electric field as a function of pot for non-degenerate case for left hand side 
f1=2*q*kT/(Ks_L*epsilon0);
fun_E_L_non_deg=@(phi) -sign(phi).*sqrt(f1*( p0_L_non_deg*(exp(-phi/kT)-1)+n0_L_non_deg*(exp(phi/kT)-1)+ (NA_L-ND_L)*phi/kT) ); %eqn 1.8
% end define electric field as a function of pot for non-degenerate case for left hand side 

% define electric field as a function of pot for non-degenerate case for right hand side 
f1=2*q*kT/(Ks_R*epsilon0);
fun_E_R_non_deg=@(phi) sign(phi).*sqrt(f1*( p0_R_non_deg*(exp(-phi/kT)-1)+n0_R_non_deg*(exp(phi/kT)-1)+ (NA_R-ND_R)*phi/kT) ); %eqn 1.8
% end define electric field as a function of pot for non-degenerate case for right hand side 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% get phi at the interface

phi_s_non_deg=fzero(@(phi) Ks_L*fun_E_L_non_deg(phi)-Ks_R*fun_E_R_non_deg(phi-Vbi_non_deg),Vbi_non_deg/2); % eqn 2.1
phi_s_deg=fzero(@(phi) Ks_L*fun_E_L_deg(phi)-Ks_R*fun_E_R_deg(phi-Vbi_deg),phi_s_non_deg);% eqn 2.1

phi_s_non_deg_exact=((ND_R-NA_R)*Vbi_non_deg+(V_thermal*(p0_L_non_deg+n0_L_non_deg-p0_R_non_deg-n0_R_non_deg)))/(NA_L-ND_L+ND_R-NA_R); %eqn 2.3a

%% end get phi at the interface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   

get_phi_E_rho_from_ode45;

phi_array_deg=[flip(phi_array_L_deg),phi_array_R_deg+Vbi_deg]; 
E_array_deg=[flip(E_array_L_deg),E_array_R_deg]; 
rho_array_deg=[flip(rho_array_L_deg),rho_array_R_deg];
n_array_deg=[flip(n_array_L_deg),n_array_R_deg];
p_array_deg=[flip(p_array_L_deg),p_array_R_deg];

phi_array_non_deg=[flip(phi_array_L_non_deg),phi_array_R_non_deg+Vbi_non_deg]; 
E_array_non_deg=[flip(E_array_L_non_deg),E_array_R_non_deg]; 
rho_array_non_deg=[flip(rho_array_L_non_deg),rho_array_R_non_deg];
n_array_non_deg=[flip(n_array_L_non_deg),n_array_R_non_deg];
p_array_non_deg=[flip(p_array_L_non_deg),p_array_R_non_deg];


plot_p_n_junction;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



