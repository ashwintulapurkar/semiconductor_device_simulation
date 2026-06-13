clear
Ks=11.9; epsilon0=8.85e-12;q=1.6e-19; kB=1.38e-23; T=300; kT=kB*T/q; V_thermal=kB*T/q; 
hbar=1.06e-34;  m0=9.1e-31; 
Eg=1.12; xi=4.03; %electron affinity in eV 

NA=0e16*1e6; ND=1e16*1e6; %choose p-type or n-type substrate

Nc=3.23e19*1e6; Nv=1.83e19*1e6; ni=sqrt(Nc*Nv)*exp(-Eg/(2*V_thermal)); %values at 300 K
Ec=-xi; Ev=Ec-Eg; Ei=(Ec+Ev)/2 + V_thermal*log(Nv/Nc)/2;

EF=-asinh(0.5*(NA-ND)/ni)*V_thermal+Ei; %non-deg soln eqn 1.2

tau_n=100e-9; tau_p=100e-9; %unit=m^2/(V-sec)
mu0_n=1360*1e-4; mu0_p=460*1e-4; v_sat=1e5; %used for electric field dependent mobility

D0_n=mu0_n*V_thermal; D0_p=mu0_p*V_thermal; %Einstein reln

tdr=(epsilon0*Ks)/(q*ND*mu0_n); %dielectric relaxation time for n-type sc
LDebye=sqrt(Ks*epsilon0*V_thermal/(ND*q)); %Debye length  for n-type sc
t1=LDebye^2/(2*D0_p);


N=max(NA,ND); %N=ni;
N1=N; N2=N; tau1=tau_n; tau2=tau_p; %normalization factors. can be changed.

VN=V_thermal; %normalization factor. can be changed.

N_ni=sqrt(N1*N2); %normalization factor for ni

Ln=sqrt(D0_n*tau1); Lp=sqrt(D0_p*tau2);

x_R=10e-6; delta_x=20e-9; 
X_array=[0:delta_x:x_R];
x_array=[delta_x:delta_x:x_R-delta_x]; Np=length(x_array);
x_array_E=X_array(1:end-1)+delta_x/2; %electric field is defined at these points

tn=(Ln/delta_x)^2; tp=(Lp/delta_x)^2; %eqn 3.9
LD=sqrt((Ks*epsilon0*VN/(N*q)));  %eqn 3.9
t_poi=(LD/delta_x)^2;   %eqn 3.11 a
array1=-2*t_poi*ones(1,Np);
array2=t_poi*ones(1,Np-1);
M3=((diag(array1))+(diag(array2,1))+(diag(array2,-1)));  %eqn. 3.13

doping_profile=(ND-NA)*ones(1,length(x_array))/N; %doping prifile normalized by N, eqn 3.11 a

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
VA=.3;
% get initial conf
V_L=0; V_R=-VA/VN;
V=linspace(V_L,V_R,length(x_array)); 
n0=Nc*exp((EF-Ec)/V_thermal); p0=Nv*exp((Ev-EF)/V_thermal);   %eqn 1.7
n=(n0/N1)*ones(1,length(x_array)); p=(p0/N2)*ones(1,length(x_array)); %n and p normalized 
n_L=n0/N1; p_L=p0/N2; n_R=n0/N1; p_R=p0/N2;
L_Debye=sqrt(Ks*epsilon0*V_thermal/((n0+p0)*q));
Jn=q*mu0_n*n0*VA/x_R; Jp=q*mu0_p*p0*VA/x_R; 

pn_junction_IV_Newton_subroutine; %get equilibrium config
f=J1*[n,p,V]'-[Rn,Rp,-rho_tilda]'+bdy'; f=f';
Jn0=Jn; Jp0=Jp;

GL0=1e18*1e6; %e-h pair per m^3 per sec produced by light
x0=4e-6; sig0=200e-9; %light spot position and size
T_pulse=1e-12; %light pulse time width
GL=GL0*exp(-(x_array-x0).^2/(2*sig0^2));

ele_field=VA/x_R; vel=mu0_p*ele_field; %used mu0_p as holes are minority

delta_t_max=min(delta_x^2/(2*D0_p),tdr); 

delta_t=T_pulse/100;
time_array_pulse=[0:delta_t:T_pulse]; %time for which light pulse in on;


for j1=1:length(time_array_pulse)
    j1
    pn_junction_IV_Newton_mod_subroutine;
end

close all;
figure(1); plot(x_array,(n*N1-n0)*1e-6,'r'); hold on; xlabel('x (m)'); ylabel('\Delta n /cc');
figure(2); plot(x_array,(p*N2-p0)*1e-6,'ro-'); hold on;  xlabel('x (m)'); ylabel('\Delta p /cc');

delta_p_max=GL0*tau_p*(1-exp(-T_pulse/tau_p)); %if T_pulse very short
delta_p_expected_after_pulse=delta_p_max*exp(-(x_array-x0).^2/(2*sig0^2)); %if T_pulse very short

A=max(delta_p_expected_after_pulse)*sig0*sqrt(2*pi);
figure(2); plot(x_array,delta_p_expected_after_pulse*1e-6,'bo-'); hold on;


GL=0; %light switched off
sim_time=1e-9;
delta_t=delta_t_max/2;
time_array1=[0:delta_t:sim_time];
for j1=1:length(time_array1)
    j1
    pn_junction_IV_Newton_mod_subroutine;
end

figure(1); plot(x_array,(n*N1-n0)*1e-6,'b'); hold on; xlabel('x (m)'); ylabel('\Delta n /cc');
figure(2); plot(x_array,(p*N2-p0)*1e-6,'b'); hold on;  xlabel('x (m)'); ylabel('\Delta p /cc');

t=time_array1(end);
x_cen=x0+vel*t;
sig1=sqrt(sig0^2+2*D0_p*t);
delta_p_expected1=A*exp(-t/tau_n)*exp(-(x_array-x_cen).^2/(2*sig1^2));
delta_p_expected1=delta_p_expected1/(sig1*(sqrt(2*pi)));

figure(2); plot(x_array,delta_p_expected1*1e-6,'m'); hold on;  xlabel('x (m)'); ylabel('\Delta p /cc');


figure(3); plot(x_array_E,Jp-Jp0,'r',x_array_E,Jp_drift-Jp0,'b',x_array_E,Jp_diff,'g'); xlabel('x (m)'); ylabel('J (A/m^2)');
legend('\delta J_p', '\delta J_{p,drift}', '\delta J_{p,diff}');



