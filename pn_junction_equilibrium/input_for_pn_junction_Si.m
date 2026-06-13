
%Si junction (pn junction. )
%hbar=1.06e-34; kB=1.38e-23; q=1.6e-19; epsilon0=8.85e-12; m0=9.1e-31; T=300; V_thermal=kB*T/q;

%%% left side semiconductor 
ND_L=0e16*1e6; NA_L=2e16*1e6; %doping in /m^3 of the left semiconductor

Ks_L=11.8; Eg_L=1.12; %dielectric constant and bandgap (Si at 300 K)   

xi_L=4.03; %electron affinity in eV 
%mn_L=1.18*m0; mp_L=0.81*m0; Nc_L=2*(mn_L*V_thermal*q/(2*pi*hbar^2))^1.5; Nv_L=2*(mp_L*V_thermal*q/(2*pi*hbar^2))^1.5; 
Nc_L=3.23e19*1e6; Nv_L=1.83e19*1e6; %effective dos of electrons and holes in /m^3
ni_L=sqrt(Nc_L*Nv_L)*exp(-Eg_L/(2*V_thermal)); %intrinsic carrier conc. eqn 1.2

Ec_L=-xi_L; Ev_L=Ec_L-Eg_L; Ei_L=(Ec_L+Ev_L)/2 + V_thermal*log(Nv_L/Nc_L)/2; %conduction band, valence band and intrinsic Fermilevel energies eqn 1.2

Ed_L=Ec_L-0.045; Ea_L=Ev_L+0.045; gD_L=2; gA_L=4; %donor and acceptor levels and their degenearcies

fun =@(EF) fun_get_charge_density(Nc_L, Nv_L, Ec_L, Ev_L, Ed_L, Ea_L, ND_L, NA_L, EF, V_thermal, gD_L, gA_L, 0).rho;
EF_L_deg = fzero(fun,Ei_L); %Fermi-level using degenerate equations. charge density=0 at correct value of EF.

EF_L_non_deg=-asinh(0.5*(NA_L-ND_L)/ni_L)*V_thermal+Ei_L; %Fermi-level using degenerate equations. eqn 1.2

n0_L_non_deg=Nc_L*exp((EF_L_non_deg-Ec_L)/V_thermal); p0_L_non_deg=Nv_L*exp((Ev_L-EF_L_non_deg)/V_thermal); %eqn 1.7
n0_L_deg=Nc_L*fermi(0.5,(EF_L_deg-Ec_L)/V_thermal); p0_L_deg=Nv_L*fermi(0.5,(Ev_L-EF_L_non_deg)/V_thermal); %eqn 1.1

LDebye_L=sqrt(Ks_L*epsilon0*V_thermal/((n0_L_non_deg+p0_L_non_deg)*q)); %Debye length, non-degenerate

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% right side semiconductor 
ND_R=1e16*1e6; NA_R=0e16*1e6; %doping in /m^3 of the right semiconductor

Ks_R=11.8; Eg_R=1.12;

xi_R=4.03;
%mn_R=1.18*m0; mp_R=0.81*m0; Nc_R=2*(mn_R*V_thermal*q/(2*pi*hbar^2))^1.5; Nv_R=2*(mp_R*V_thermal*q/(2*pi*hbar^2))^1.5;
Nc_R=3.23e19*1e6; Nv_R=1.83e19*1e6;
ni_R=sqrt(Nc_R*Nv_R)*exp(-Eg_R/(2*V_thermal));
Ec_R=-xi_R; Ev_R=Ec_R-Eg_R; Ei_R=(Ec_R+Ev_R)/2 + V_thermal*log(Nv_R/Nc_R)/2;

Ed_R=Ec_R-0.045; Ea_R=Ev_R+0.045;  gD_R=2; gA_R=4;

fun =@(EF) fun_get_charge_density(Nc_R, Nv_R, Ec_R, Ev_R, Ed_R, Ea_R, ND_R, NA_R, EF, V_thermal, gD_R, gA_R, 0).rho;
EF_R_deg = fzero(fun,Ei_R); %deg soln

EF_R_non_deg=-asinh(0.5*(NA_R-ND_R)/ni_R)*V_thermal+Ei_R;  %non-deg soln

n0_R_non_deg=Nc_R*exp((EF_R_non_deg-Ec_R)/V_thermal); p0_R_non_deg=Nv_R*exp((Ev_R-EF_R_non_deg)/V_thermal); 
n0_R_deg=Nc_R*fermi(0.5,(EF_R_deg-Ec_R)/V_thermal); p0_R_deg=Nv_R*fermi(0.5,(Ev_R-EF_R_deg)/V_thermal); 

LDebye_R=sqrt(Ks_R*epsilon0*V_thermal/((n0_R_non_deg+p0_R_non_deg)*q)); %Debye length

