%n=density of electrons/N1;  p=density of holes/N2
% V=potential/V_N, symbol phi is used to denote potential in the document.
%doping=(ND-NA)/N
%rho =charge density/q
%U=potential/V_thermal
% values of normalization constants, N1, N2, N, tau1, tau2 can be set in the code

clear
Ks=11.9; epsilon0=8.85e-12;q=1.6e-19; kB=1.38e-23; T=300; V_thermal=kB*T/q; kT=V_thermal;
hbar=1.06e-34;  m0=9.1e-31; 
Eg=1.12; xi=4.03; %electron affinity in eV 

NA=1e18*1e6; ND=1e16*1e6; %doping in /m^3
Nc=3.23e19*1e6; Nv=1.83e19*1e6;  %effective dos of electrons and holes in /m^3

ni=sqrt(Nc*Nv)*exp(-Eg/(2*V_thermal));  %intrinsic carrier conc. eqn 1.2
Ec=-xi; Ev=Ec-Eg; Ei=(Ec+Ev)/2 + V_thermal*log(Nv/Nc)/2; %eqn 1.2

EF_L=-asinh(0.5*NA/ni)*V_thermal+Ei; %non-deg soln eqn 1.2
EF_R=asinh(0.5*ND/ni)*V_thermal+Ei; %non-deg soln eqn 1.2
Vbi=(EF_R-EF_L);

tau_n=0.1e-6; tau_p=0.1e-6; %unit=m^2/(V-sec)

mu0_n=1360*1e-4; mu0_p=460*1e-4; v_sat=1e5; %used for electric field dependent mobility

D0_n=mu0_n*V_thermal; D0_p=mu0_p*V_thermal; %Einstein reln

%N1=ND; N2=NA; %normalization factors. can be changed.
tau1=tau_n; tau2=tau_p; %normalization factors. can be changed.
N=max(NA,ND);N1=N; N2=N; %normalization factors. can be changed.

VN=V_thermal; %normalization factor. can be changed.

N_ni=sqrt(N1*N2); %normalization factor for ni

Ln=sqrt(D0_n*tau1); Lp=sqrt(D0_p*tau2); % eqn 3.9

LDebye_L=sqrt(Ks*epsilon0*V_thermal/(NA*q)); %Debye length  for n-type sc
LDebye_R=sqrt(Ks*epsilon0*V_thermal/(ND*q)); %Debye length  for n-type sc

x_L=-5e-6; x_R=15e-6; delta_x=20e-9; %widths of the left and right side of the junction.
% note: LDebye_L=4 nm for NA=1e18. But for p+n junction, n side is more
% important. So 20 nm is chosen as delta_x. You can reduce delta_x and
% check if there is any difference.

X_array=[x_L:delta_x:x_R];
x_array=[x_L+delta_x:delta_x:x_R-delta_x]; Np=length(x_array);
x_array_E=X_array(1:end-1)+delta_x/2; %electric field is defined at these points. Midway between X_array.

tn=(Ln/delta_x)^2; %eqn 3.9
tp=(Lp/delta_x)^2; %eqn 3.9

get_equi_V_n_p; %starting n,p and V=equilibrium values

delta_V=V_thermal/3; %increase voltage in small steps
VA_array=[0:delta_V:1.2]; VA_array=[VA_array,1.2];
%VA_array=[0:-delta_V:-1]; VA_array=[VA_array,-1];
Jn=0; Jp=0;
tic
for iVA=1:length(VA_array)
    VA=VA_array(iVA)
    pn_junction_IV_Gummel_subroutine; 
    %pn_junction_IV_Newton_subroutine; 
    
    J_array(iVA)=mean(J_total(4:end-4));
    Efn_minus_Efp_array(iVA)=Efn_minus_Efp; 
      
%     omega=1/tau_p; get_diode_admittance; % admittance is calculated at each voltage point. Set the value of omega you want.
%     G_array(iVA)=conductance; C_array(iVA)=capacitance;
%      
end
toc

close all;
plot_vs_voltage;  % plot various quantities as a function of voltage. e.g. current vs voltage
plot_for_last_voltage; % plot various quantities for the last voltage point e.g. current vs x, n vs x etc.
plot_diode_theory; %add depletion approx plots. Compare numerical results with depl approx.
return

% plot diode admittance. For admittance calculation,uncomment lines 69 and 70.
close all; %if you want to close the dc IV plots.
plot_diode_admittance_vs_VA; % small signal Y plots as a function of VA (for a fixed freq chosen on line 69)
plot_diode_admittance_vs_freq; % small signal Y plots as a function of freq, for the last VA. Can choose freq range in the code.

%%% small signal Y plots for the last VA and one freq. To see that real and imaginary parts of current are conserved.
omega=10/tau_p;
plot_diode_admittance_last_voltage_one_freq;


