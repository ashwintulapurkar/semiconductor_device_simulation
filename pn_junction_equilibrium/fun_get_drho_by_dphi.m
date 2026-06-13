function [drho_by_dphi]=fun_get_drho_by_dphi(Nc, Nv, Ec, Ev, Ed, Ea, Nd, Na, EF, V_thermal, gD, gA, phi)

eta_c=(EF-Ec+phi)/V_thermal; dn_by_dphi=Nc*fermi(-0.5,eta_c)/V_thermal;  %eqn 2.10
eta_v=(Ev-phi-EF)/V_thermal; dp_by_dphi=-Nv*fermi(-0.5,eta_v)/V_thermal; %eqn 2.10

r1=1./(1+gD*exp((EF-Ed+phi)/V_thermal)); 
Ndplus=Nd.*r1; 
dNdplus_by_dphi=-(Ndplus).*(1-r1)/V_thermal;  %eqn 2.9

r2=1./(1+gA*exp((Ea-phi-EF)/V_thermal));
Naminus=Na.*r2;
dNaminus_by_dphi=(Naminus).*(1-r2)/V_thermal;  %eqn 2.9

drho_by_dphi=dp_by_dphi-dn_by_dphi+dNdplus_by_dphi-dNaminus_by_dphi;  %eqn 2.11


end