clear
Ks=11.9; epsilon0=8.85e-12;q=1.6e-19; kB=1.38e-23; T=300; kT=kB*T/q; V_thermal=kB*T/q; %V_thermal=.0259
hbar=1.06e-34;  m0=9.1e-31; 
Eg=1.12; xi=4.03; %electron affinity in eV 

Kox=3.9; tox=10e-9;
Cox=epsilon0*Kox/tox;

NA=1e15*1e6; ND=0e16*1e6; %choose p-type or n-type substrate

Nc=3.23e19*1e6; Nv=1.83e19*1e6; ni=sqrt(Nc*Nv)*exp(-Eg/(2*V_thermal)); %values at 300 K
Ec=-xi; Ev=Ec-Eg; Ei=(Ec+Ev)/2 + V_thermal*log(Nv/Nc)/2;
Ed=Ec-0.045; Ea=Ev+0.045; gD=2; gA=4;

fun =@(EF) fun_get_charge_density(Nc, Nv, Ec, Ev, Ed, Ea, ND, NA, EF, V_thermal, gD, gA, 0).rho;
EF_deg = fzero(fun,Ei);
EF_non_deg=-asinh(0.5*(NA-ND)/ni)*V_thermal+Ei; %non-deg soln

phi_F=Ei-EF_non_deg;


tau_n=.1e-6; tau_p=.1e-6; %unit=m^2/(V-sec)
mu0_n=1360*1e-4; mu0_p=460*1e-4; v_sat=1e5; %used for electric field dependent mobility

D0_n=mu0_n*V_thermal; D0_p=mu0_p*V_thermal; %Einstein reln

N=max(NA,ND);
N1=N; N2=N; tau1=tau_n; tau2=tau_p; %normalization factors. can be changed.

VN=V_thermal; %normalization factor. can be changed.

N_ni=sqrt(N1*N2); %normalization factor for ni

Ln=sqrt(D0_n*tau1); Lp=sqrt(D0_p*tau2); % eqn 3.9

x_R=2e-6; delta_x=1e-9;
X_array=[0:delta_x:x_R];
x_array=[delta_x:delta_x:x_R-delta_x]; Np=length(x_array);
x_array_E=X_array(1:end-1)+delta_x/2; %electric field is defined at these points

tn=(Ln/delta_x)^2; tp=(Lp/delta_x)^2; %eqn 3.9

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define functions (non-degenerate case)
p0=ni*exp((-EF_non_deg+Ei)/kT); n0=ni*exp((EF_non_deg-Ei)/kT);  %eqn 1.7
f1=2*q*kT/(Ks*epsilon0);

fun_E_sub_non_deg=@(phi) sign(phi).*sqrt(-f1*( p0*(1-exp(-phi/kT))+n0*(1-exp(phi/kT))+ (ND-NA)*phi/kT) ); %eqn 1.8, electric field in substrate
fun_get_VG_non_deg=@(phi) phi+ (Ks/Kox)*fun_E_sub_non_deg(phi)*tox;
fun_dphi_by_dz_non_deg=@(z,phi) -fun_E_sub_non_deg(phi);
%%%%%%%%%%%%%%%% end define functions (non-degenerate case)

% define functions (degenerate case)
fun_f1_R_deg=@(phi)fermi(1.5,(Ev-EF_deg-phi)/kT)-fermi(1.5,(Ev-EF_deg)/kT);  %eqn 1.9
fun_f2_R_deg=@(phi)fermi(1.5,(EF_deg-Ec+phi)/kT)-fermi(1.5,(EF_deg-Ec)/kT); %eqn 1.10
fun_f3_R_deg=@(phi) phi/kT+log( (1+gD*exp((EF_deg-Ed)/kT))./ (1+gD*exp((EF_deg-Ed+phi)/kT))); %eqn 1.13
fun_f4_R_deg=@(phi) phi/kT+log( (1+gA*exp((Ea-EF_deg-phi)/kT))./ (1+gA*exp((Ea-EF_deg)/kT))); %eqn 1.14
fun_integral_rho_R=@(phi) kT*( -Nv*fun_f1_R_deg(phi) -Nc*fun_f2_R_deg(phi) + ND*fun_f3_R_deg(phi) - NA*fun_f4_R_deg(phi)); %eqn 1.15

fun_E_sub_deg=@(phi)  sign(phi).*sqrt(-2/(epsilon0*Ks)*q*fun_integral_rho_R(phi));  %eqn 1.5, electric field in substrate

fun_get_VG_deg=@(phi) phi+ (Ks/Kox)*fun_E_sub_deg(phi)*tox;
fun_dphi_by_dz_deg=@(z,phi) -fun_E_sub_deg(phi);
% end define functions (degenerate case)

Vth_deg=fun_get_VG_deg(2*phi_F); %threshold voltage


%let's find out VG for which Ec(surface)-EF=3kT and EF-Ev(surface)=3kT
phi_s_lim_Ec_3kT=phi_F + Eg/2 -3*kT; VG_lim_Ec_3kT=fun_get_VG_deg(phi_s_lim_Ec_3kT);
phi_s_lim_Ev_3kT=phi_F - Eg/2 +3*kT; VG_lim_Ev_3kT=fun_get_VG_deg(phi_s_lim_Ev_3kT);
%end find out VG for which EC(surface)-EF=3kT and EF-EV(surface)=3kT

%let's find out VG for which Ec(surface)-EF=0 and EF-Ev(surface)=0
phi_s_lim_Ec_0kT=phi_F + Eg/2 -0*kT; VG_lim_Ec_0kT=fun_get_VG_deg(phi_s_lim_Ec_0kT);
phi_s_lim_Ev_0kT=phi_F - Eg/2 +0*kT; VG_lim_Ev_0kT=fun_get_VG_deg(phi_s_lim_Ev_0kT);
%end find out VG for which EC(surface)-EF=0 and EF-EV(surface)=0


doping_profile=(ND-NA)*ones(1,length(x_array));
doping_profile=doping_profile/N; %doping prifile normalized
LD=sqrt((Ks*epsilon0*VN/(N*q)));  %eqn 3.11 a
t_poi=(LD/delta_x)^2; %Vbdy=t_poi*[V_L,zeros(1,Np-2),V_R];

array1=-2*t_poi*ones(1,Np); array2=t_poi*ones(1,Np-1);
M3=((diag(array1))+(diag(array2,1))+(diag(array2,-1)));  %eqn. 3.13

VG_array=[linspace(-2,3,60)];

freq=100; 
omega=2*pi*freq;

 opts = odeset('RelTol',1e-8,'AbsTol',1e-8);
for iVG=1:length(VG_array)
    VG=VG_array(iVG)
    phi_s_deg=fzero(@(phi_s) Kox*(VG-phi_s)/tox- Ks*fun_E_sub_deg(phi_s), VG/2) ;
    phi_s_non_deg=fzero(@(phi_s) Kox*(VG-phi_s)/tox- Ks*fun_E_sub_non_deg(phi_s), VG/2) ;
    phi_s_deg_array(iVG)=phi_s_deg; phi_s_non_deg_array(iVG)=phi_s_non_deg;


get_equi_npV_non_deg; get_MOS_capacitance_non_deg_sub; 

%get_equi_npV_deg; get_MOS_capacitance_deg_sub; 
    

 cap_array(iVG)=capacitance;  cond_array(iVG)=conductance; 
end


close all;
figure(1); plot(VG_array,cap_array,'r','LineWidth',2); hold on; ylabel('C (F/m^2)');
xline(VG_lim_Ec_0kT); xline(VG_lim_Ev_0kT);  xline(Vth_deg);
set(gca,'fontsize', 14);

figure(2); plot(VG_array,cond_array,'b','LineWidth',2); hold on;  ylabel('G (S/m^2)');
xline(VG_lim_Ec_0kT); xline(VG_lim_Ev_0kT);  xline(Vth_deg);
set(gca,'fontsize', 14);

theta=angle(delta_VG); 
delta_Jn=delta_Jn*exp(-i*theta);
delta_Jp=delta_Jp*exp(-i*theta);
delta_JD=delta_JD*exp(-i*theta);
delta_J=delta_J*exp(-i*theta);
delta_n=delta_n*exp(-i*theta);
delta_p=delta_p*exp(-i*theta);
delta_E=delta_E*exp(-i*theta);


figure(3); plot(x_array_E,real(delta_Jn),'r',x_array_E,real(delta_Jp),'b',x_array_E,real(delta_JD),'g',x_array_E,real(delta_J),'k');
xlabel('x (m)'); ylabel('real \delta J (A/m^2)'); 
legend('real(\delta J_n)','real(\delta J_p)','real(\delta J_D)','real(\delta J)');
title(['V_G=',num2str(VG), ', f=',num2str(freq),' Hz']);

figure(4); plot(x_array_E,imag(delta_Jn),'ro',x_array_E,imag(delta_Jp),'bx',x_array_E,imag(delta_JD),'g',x_array_E,imag(delta_J),'k');
xlabel('x (m)'); ylabel('imag \delta J (A/m^2)');
legend('imag(\delta J_n)','imag(\delta J_p)','imag(\delta J_D)','imag(\delta J)');
title(['V_G=',num2str(VG), ', f=',num2str(freq),' Hz']);

