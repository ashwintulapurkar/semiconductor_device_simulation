%solve eqn d_phi/dx=-E(phi) using ODE45. Need phi_s as boundary condition. 

fun_dphi_by_dx_L_deg=@(x,phi) -fun_E_L_deg(phi);
fun_dphi_by_dx_R_deg=@(x,phi) -fun_E_R_deg(phi);

fun_dphi_by_dx_L_non_deg=@(x,phi) -fun_E_L_non_deg(phi);
fun_dphi_by_dx_R_non_deg=@(x,phi) -fun_E_R_non_deg(phi);


opts = odeset('RelTol',1e-6,'AbsTol',1e-6);
%%% deg soln
[z_tilda_array_R, phi_array_R_deg]=ode45(fun_dphi_by_dx_R_deg,x_array_R,phi_s_deg-Vbi_deg,opts); 
phi_array_R_deg=real(phi_array_R_deg');
E_array_R_deg=fun_E_R_deg(phi_array_R_deg); 
zz=fun_get_charge_density(Nc_R, Nv_R, Ec_R, Ev_R, Ed_R, Ea_R, ND_R, NA_R, EF_R_deg, V_thermal, gD_R, gA_R, phi_array_R_deg);
rho_array_R_deg=zz.rho; n_array_R_deg=zz.n; p_array_R_deg=zz.p;

[z_tilda_array_L, phi_array_L_deg]=ode45(fun_dphi_by_dx_L_deg,x_array_L,phi_s_deg,opts); 
phi_array_L_deg=real(phi_array_L_deg');
E_array_L_deg=fun_E_L_deg(phi_array_L_deg); 
zz=fun_get_charge_density(Nc_L, Nv_L, Ec_L, Ev_L, Ed_L, Ea_L, ND_L, NA_L, EF_L_deg, V_thermal, gD_L, gA_L, phi_array_L_deg);
rho_array_L_deg=zz.rho; n_array_L_deg=zz.n; p_array_L_deg=zz.p;
%%%%%%%%%%%%%%%%%
%non deg soln
[z_tilda_array_R, phi_array_R_non_deg]=ode45(fun_dphi_by_dx_R_non_deg,x_array_R,phi_s_non_deg-Vbi_non_deg,opts); 
phi_array_R_non_deg=real(phi_array_R_non_deg');
E_array_R_non_deg=fun_E_R_non_deg(phi_array_R_non_deg); 
n_array_R_non_deg=n0_R_non_deg*exp(phi_array_R_non_deg/kT);
p_array_R_non_deg=p0_R_non_deg*exp(-phi_array_R_non_deg/kT);
rho_array_R_non_deg=ND_R-NA_R+p_array_R_non_deg-n_array_R_non_deg;


[z_tilda_array_L, phi_array_L_non_deg]=ode45(fun_dphi_by_dx_L_non_deg,x_array_L,phi_s_non_deg,opts); 
phi_array_L_non_deg=real(phi_array_L_non_deg');
E_array_L_non_deg=fun_E_L_non_deg(phi_array_L_non_deg); 
n_array_L_non_deg=n0_L_non_deg*exp(phi_array_L_non_deg/kT);
p_array_L_non_deg=p0_L_non_deg*exp(-phi_array_L_non_deg/kT);
rho_array_L_non_deg=ND_L-NA_L+p_array_L_non_deg-n_array_L_non_deg;



