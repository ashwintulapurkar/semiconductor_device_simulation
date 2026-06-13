
%Al_0.7 Ga_0.3 As/GaAs junction (n+p junction. )
%put correct values of parameters

kB=1.38e-23; T=300; V_thermal=kB*T/q; 

Ks_L=12.90-2.84*.3; Eg_L=1.424+(1.247*.3); 
ND_L=5e17*1e6; NA_L=0*1e15*1e6; %doping in /m^3
xi_L=4.07-1.1*.3; %electron affinity in eV i.e. Evacuum-Ec
xi_L=4.07-.23; %put here correct value


Nc_L=4.21e17*1e6; Nv_L=9.52e18*1e6; ni_L=sqrt(Nc_L*Nv_L)*exp(-Eg_L/(2*V_thermal));

Ec_L=-xi_L; Ev_L=Ec_L-Eg_L; Ei_L=(Ec_L+Ev_L)/2 + V_thermal*log(Nv_L/Nc_L)/2;

Ed_L=Ec_L-0.0078; Ea_L=Ev_L+0.0078;   gD_L=2; gA_L=4; %put here correct values

fun =@(EF) fun_get_charge_density(Nc_L, Nv_L, Ec_L, Ev_L, Ed_L, Ea_L, ND_L, NA_L, EF, V_thermal, gD_L, gA_L, 0).rho;
EF_L_deg = fzero(fun,Ei_L); %deg soln

EF_L_non_deg=-asinh(0.5*(NA_L-ND_L)/ni_L)*V_thermal+Ei_L; %non-deg soln

n0_L_non_deg=Nc_L*exp((EF_L_non_deg-Ec_L)/V_thermal); p0_L_non_deg=Nv_L*exp((Ev_L-EF_L_non_deg)/V_thermal); 
LD_L=sqrt(Ks_L*epsilon0*V_thermal/((n0_L_non_deg+p0_L_non_deg)*q)); %Debye length

n0_L_deg=Nc_L*fermi(0.5,(EF_L_deg-Ec_L)/V_thermal); p0_L_deg=Nv_L*fermi(0.5,(Ev_L-EF_L_deg)/V_thermal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ks_R=12.9; Eg_R=1.424;
ND_R=0e11*1e6; NA_R=1e14*1e6;
xi_R=4.07;

Nc_R=4.21e17*1e6; Nv_R=9.52e18*1e6; ni_R=sqrt(Nc_R*Nv_R)*exp(-Eg_R/(2*V_thermal));
Ec_R=-xi_R; Ev_R=Ec_R-Eg_R; Ei_R=(Ec_R+Ev_R)/2 + V_thermal*log(Nv_R/Nc_R)/2;

Ed_R=Ec_R-0.0078; Ea_R=Ev_R+0.0078;   gD_R=2; gA_R=4;

fun =@(EF) fun_get_charge_density(Nc_R, Nv_R, Ec_R, Ev_R, Ed_R, Ea_R, ND_R, NA_R, EF, V_thermal, gD_R, gA_R, 0).rho;
EF_R_deg = fzero(fun,Ei_R); %deg soln

EF_R_non_deg=-asinh(0.5*(NA_R-ND_R)/ni_R)*V_thermal+Ei_R;  %non-deg soln

n0_R_non_deg=Nc_R*exp((EF_R_non_deg-Ec_R)/V_thermal); p0_R_non_deg=Nv_R*exp((Ev_R-EF_R_non_deg)/V_thermal); %reqd for non-deg case
LD_R=sqrt(Ks_R*epsilon0*V_thermal/((n0_R_non_deg+p0_R_non_deg)*q)); %Debye length

n0_R_deg=Nc_L*fermi(0.5,(EF_R_non_deg-Ec_R)/V_thermal); p0_R_deg=Nv_L*fermi(0.5,(Ev_R-EF_R_non_deg)/V_thermal);

%%%%%%%%%%%%
