

xn=sqrt(2*Ks*epsilon0*NA*(Vbi-VA)/(q*ND*(NA+ND))); xp=sqrt(2*Ks*epsilon0*ND*(Vbi-VA)/(q*NA*(NA+ND))); w=xn+xp;
W1=-x_L-xp; W2=x_R-xn; %widths of the quasi-neutral sides on p and n sides

D_n=(kB*T/q)*mu0_n; D_p=(kB*T/q)*mu0_p; %L_n=sqrt(D_n*tau_n); L_p=sqrt(D_p*tau_p);
tau_n_eff=tau_n./(1-i*omega_array*tau_n); tau_p_eff=tau_p./(1-i*omega_array*tau_p);
L_n=sqrt(D_n*tau_n_eff); L_p=sqrt(D_p*tau_p_eff);

Js=q*((D_p./L_p).*coth(W2./L_p)*(ni^2/ND)+(D_n./L_n).*coth(W1./L_n)*(ni^2/NA)); %reverse saturation current density for narrow diode

Y_diffusion=(Js*q/(kB*T)).*(exp(q*VA/(kB*T))); %diffusion admittance %see eqns 7.26-7.28 in Pierret
G_diffusion=real(Y_diffusion); C_diffusion=-imag(Y_diffusion)./omega_array;
C_junction=epsilon0*Ks./w; 

tau=tau_n; f1=q*ni/(2*tau)*w; %for R.B.
f2=(exp(VA/V_thermal)-1);
f3=0.5*exp(0.5*VA/V_thermal).*((Vbi-VA))/V_thermal;
J_RG=(f1.*f2)./(1+f3);

wdash=-sqrt((2*Ks*epsilon0/q)*(NA+ND)/(NA*ND))/(2*sqrt(Vbi-VA)); f1dash=q*ni/(2*tau)*wdash;
f2dash=exp(VA/V_thermal)/V_thermal;
f3dash=0.5*exp(0.5*VA/V_thermal).*(-1)/V_thermal + 0.25*exp(0.5*VA/V_thermal).*((Vbi-VA))/V_thermal^2;

G_RG=f1dash.*f2./(1+f3)+f1.*f2dash./(1+f3)-f1.*f2.*f3dash./(1+f3).^2;  %low freq conductance from RG current

C_junction=C_junction*ones(1,length(omega_array));
G_RG=G_RG*ones(1,length(omega_array));

