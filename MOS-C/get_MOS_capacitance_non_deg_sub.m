
delta_n_R=0;  delta_p_R=0;
delta_n_L=0; delta_p_L=0; %delta_n_L and delta_p_L are arb
delta_phi_s=1e-3; %value of delta_phi_s doesnot matter while calculating capacitance

delta_V_L=delta_phi_s/VN;  %delta_V_L normalized
delta_V_R=0; %delta_V_R is normalized

get_Jacobian_non_deg;

M1(1,1)=-tn*(Dn(2)/D0_n)*B((V(1)-V(2))*VN/kT); %from bc that Jn=0 on LHS
J2_13(1,1)=-hn(2);  %from bc that Jn=0 on LHS

M2(1,1)=-tp*(Dp(2)/D0_p)*B((V(2)-V(1))*VN/kT); %from bc that Jp=0 on LHS
J2_23(1,1)=-hp(2); %from bc that Jp=0 on LHS

J1=blkdiag(M1,M2,M3);
J2=[J2_11,J2_12,J2_13; J2_21,J2_22,J2_23; J2_31,J2_32,zeros(Np,Np)];

J3=i*omega*blkdiag(tau1*eye(Np),tau2*eye(Np),zeros(Np,Np));



nbdy=[0,zeros(1,Np-2),delta_V_R*hn(end)];
pbdy=[0,zeros(1,Np-2),delta_V_R*hp(end)];
Vbdy=t_poi*[delta_V_L,zeros(1,Np-2),delta_V_R];

bdy=[nbdy,pbdy,Vbdy];


delta_n_p_V=(J1+J2+J3)\transpose(-bdy); 

delta_n_p_V=transpose(delta_n_p_V);
delta_n=delta_n_p_V(1:Np); delta_p=delta_n_p_V(Np+1:2*Np); delta_V=delta_n_p_V(2*Np+1:end);


delta_E=-diff([delta_V_L,delta_V,delta_V_R]*VN)./delta_x;
diff_delta_U=diff([delta_V_L,delta_V,delta_V_R]*(VN/V_thermal));
diff_delta_V=diff([delta_V_L,delta_V,delta_V_R]);


n_all=[delta_n_L,delta_n,delta_n_R]; p_all=[delta_p_L,delta_p,delta_p_R];
delta_Jn=N1*q/(delta_x)*Dn.*(-n_all(1:end-1).*B_minusU_diff + n_all(2:end).*B_U_diff);
delta_Jp=-N2*q/(delta_x)*Dp.*(-p_all(1:end-1).*B_U_diff + p_all(2:end).*B_minusU_diff);

delta_Jn=delta_Jn+(N1*q*delta_x/tau1)*diff_delta_V.*hn;
delta_Jp=delta_Jp-(N2*q*delta_x/tau2)*diff_delta_V.*hp;

delta_Jn(1)=0; delta_Jp(1)=0;

delta_JD=-Ks*epsilon0*i*omega*delta_E; %assume exp(-iwt) dependence

delta_J=delta_Jn+delta_Jp+delta_JD;


delta_J1=mean(delta_J(100:end-10)); %delta_J can be noisy on LHS

delta_E_Si=i*(delta_J1)/(Ks*epsilon0*omega); %at interface JD=J

delta_E_ox=delta_E_Si*Ks/Kox;
%delta_E_ox=delta_E(1)*Ks/Kox;
delta_VG=delta_E_ox*tox+delta_V_L*VN;   

Y=delta_J1/delta_VG;
conductance=real(Y); capacitance=-imag(Y)/omega;







