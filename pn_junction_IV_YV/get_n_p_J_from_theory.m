%all in SI

D_n=(kB*T/q)*mu0_n; D_p=(kB*T/q)*mu0_p; L_n=sqrt(D_n*tau_n); L_p=sqrt(D_p*tau_p); %minority carrier diffusion lengths

 xn=sqrt(2*Ks*epsilon0*NA*(Vbi-VA)/(q*ND*(NA+ND))); xp=sqrt(2*Ks*epsilon0*ND*(Vbi-VA)/(q*NA*(NA+ND))); w=xn+xp;

[V_depletion_approx_array,E_depletion_approx_array, rho_depletion_approx_array,wL,wR]=fun_get_V_E_rho_from_depletion_approx(X_array,Ks,Ks, NA,ND, VA, Vbi);


W1=-X_array(1)-xp; W2=X_array(end)-xn; %widths of the quasi-neutral sides on p and n sides


 Ec_depl=Ec-V_depletion_approx_array; Ev_depl=Ev-V_depletion_approx_array; Ei_depl=Ei-V_depletion_approx_array; 
 
EF_L=(Ei_depl(1)-V_thermal*log(NA/ni)); EF_R=(Ei_depl(end)+V_thermal*log(ND/ni));
%%%%%%%%%%%%%%%%%%%%%%%%% calculate n and p
 delta_n_0=(ni^2/NA)*(exp(q*VA/(kB*T))-1); delta_p_0=(ni^2/ND)*(exp(q*VA/(kB*T))-1);
delta_n=delta_n_0*(sinh((-X_array(1)+X_array)/L_n)/sinh(W1/L_n)).*(1-heaviside(X_array+xp)); %for narrow diode 

delta_p=delta_p_0*(sinh((X_array(end)-X_array)/L_p)/sinh(W2/L_p)).*heaviside(X_array-xn); %for narrow diode 

n_depl=(ni^2/NA)*(1-heaviside(X_array+xp))+delta_n;
p_depl=(ni^2/ND)*heaviside(X_array-xn)+delta_p;

n_depl=n_depl+ND*heaviside(X_array-xn); %1st term=n on p-side, 2nd term=n on n-side
p_depl=p_depl+NA*(1-heaviside(X_array+xp));
%let's find out n and p in the depletion regions
FN_depl_reg=EF_R; 
FP_depl_reg=EF_L;

n_depl=n_depl+ni*exp((FN_depl_reg-Ei_depl)/V_thermal).*(1-heaviside(X_array-xn)).*heaviside(X_array+xp);
p_depl=p_depl+ni*exp((Ei_depl-FP_depl_reg)/V_thermal).*(1-heaviside(X_array-xn)).*heaviside(X_array+xp);
%%%%%%%%%%%%%%%%%%%%%%%%% end calculate n and p
%%%%%%%%%%%%%%%%%%%%%% calculate Jn Jp J

Js=q*((D_p/L_p)*coth(W2/L_p)*(ni^2/ND)+(D_n/L_n)*coth(W1/L_n)*(ni^2/NA)); %reverse saturation current density for narrow diode

Jn_0=q*((D_n/L_n)*(ni^2/NA))*(exp(q*VA/(kB*T))-1); Jp_0=q*((D_p/L_p)*(ni^2/ND))*(exp(q*VA/(kB*T))-1);  %for long diode
Jn_0=Jn_0*coth(W1/L_n); Jp_0=Jp_0*coth(W2/L_p); % for narrow diode


J_1=Js*(exp(q*VA/(kB*T))-1)*ones(1,length(X_array));
%avoid using J as variable name. used J_1
Jn_1=Jn_0*(cosh((-X_array(1)+X_array)/L_n)/cosh(W1/L_n)).*(1-heaviside(X_array+xp))+Jn_0*heaviside(X_array+xp).*(1-heaviside(X_array-xn)); %narrow diode
%avoid using Jn as variable name. used Jn_1
Jp_1=Jp_0*(cosh((X_array(end)-X_array)/L_p)/cosh(W2/L_p)).*(heaviside(X_array-xn))+Jp_0*heaviside(X_array+xp).*(1-heaviside(X_array-xn)); %narrow diode
%avoid using Jp as variable name.  used Jp_1
Jn_1=Jn_1+(J_1-Jp_1).*heaviside(X_array-xn);
Jp_1=Jp_1+(J_1-Jn_1).*(1-heaviside(X_array+xp));
%%%%%%%%%%%%%%%%%%%%%%%% end calculate Jn Jp J
%% calculate J as a function of VA
J_array1=Js*(exp(q*VA_array/(kB*T))-1);


%% add R-G in the depletion region
 xn_array=sqrt(2*Ks*epsilon0*NA*(Vbi-VA_array)/(q*ND*(NA+ND))); xp_array=sqrt(2*Ks*epsilon0*ND*(Vbi-VA_array)/(q*NA*(NA+ND))); w_array=xn_array+xp_array;
tau=(tau_n+tau_p)/2;
f1=q*ni/(2*tau)*w_array; %for R.B.
f2=(exp(VA_array/V_thermal)-1);
f3=0.5*exp(0.5*VA_array/V_thermal).*((Vbi-VA_array))/V_thermal;
J_RG_array=(f1.*f2)./(1+f3); %eqn 6.45 in Pierret

Jn_depl_approx=Jn_1; Jp_depl_approx=Jp_1; J_depl_approx=J_1;


FN_depl=Ei_depl+V_thermal*log(n_depl/ni);
FP_depl=Ei_depl-V_thermal*log(p_depl/ni);


%%%%
