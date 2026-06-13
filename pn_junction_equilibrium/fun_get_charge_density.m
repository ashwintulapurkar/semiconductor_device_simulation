function all_rho=fun_get_charge_density(Nc, Nv, Ec, Ev, Ed, Ea, Nd, Na, EF, V_thermal, gD, gA, phi)

eta_c=(EF-Ec+phi)/V_thermal; n=Nc*fermi(0.5,eta_c); 
eta_v=(Ev-phi-EF)/V_thermal; p=Nv*fermi(0.5,eta_v); 

%gD=2; 
Ndplus=Nd./(1+gD*exp((EF-Ed+phi)/V_thermal)); %eq 4.65 in Advanced semiconductor fundamentals by Pierret

%gA=4; 
Naminus=Na./(1+gA*exp((Ea-phi-EF)/V_thermal));

rho=p-n+Ndplus-Naminus;

all_rho.rho=rho;
all_rho.n=n;
all_rho.p=p;
all_rho.Ndplus=Ndplus;
all_rho.Naminus=Naminus;


end