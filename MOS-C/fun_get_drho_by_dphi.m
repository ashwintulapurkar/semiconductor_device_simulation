function [drho_by_dphi]=fun_get_drho_by_dphi(Nc, Nv, Ec, Ev, Ed, Ea, Nd, Na, EF, V_thermal, gD, gA, phi)

%[drho_by_dphi, dNdplus_by_dphi,dNaminus_by_dphi]
eta_c=(EF-Ec+phi)/V_thermal; dn_by_dphi=Nc*fermi(-0.5,eta_c)/V_thermal; 
eta_v=(Ev-phi-EF)/V_thermal; dp_by_dphi=-Nv*fermi(-0.5,eta_v)/V_thermal;  

r1=1./(1+gD*exp((EF-Ed+phi)/V_thermal)); 
Ndplus=Nd.*r1; 
dNdplus_by_dphi=-(Ndplus).*(1-r1)/V_thermal;

r2=1./(1+gA*exp((Ea-phi-EF)/V_thermal));
Naminus=Na.*r2;
dNaminus_by_dphi=(Naminus).*(1-r2)/V_thermal;

drho_by_dphi=dp_by_dphi-dn_by_dphi+dNdplus_by_dphi-dNaminus_by_dphi;


end