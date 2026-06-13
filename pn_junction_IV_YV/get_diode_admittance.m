
delta_n_L=0; delta_n_R=0; delta_p_L=0; delta_p_R=0; %eqn 4.3
delta_VA=10e-3; %value of delta_VA doesnot matter while calculating admittance
%delta_V_L=(delta_VA/2)/VN; delta_V_R=-(delta_VA/2)/VN; %delta_V_L normalized, %delta_V_R is normalized
delta_V_L=0/VN; delta_V_R=-delta_VA/VN;

get_Jacobian;
J3=i*omega*blkdiag(tau1*eye(Np),tau2*eye(Np),zeros(Np,Np)); %eqn 4.5, second term on LHS


nbdy=[delta_V_L*hn(1),zeros(1,Np-2),delta_V_R*hn(end)]; %eqn 4.5, second term on RHS
pbdy=[delta_V_L*hp(1),zeros(1,Np-2),delta_V_R*hp(end)]; %eqn 4.5, second term on RHS
Vbdy=t_poi*[delta_V_L,zeros(1,Np-2),delta_V_R]; %eqn 4.5, second term on RHS


bdy=[nbdy,pbdy,Vbdy];
delta_n_p_V=(J1+J2+J3)\transpose(-bdy); %solution of eqn 4.5 
delta_n_p_V=transpose(delta_n_p_V); %column to row
delta_n=delta_n_p_V(1:Np); delta_p=delta_n_p_V(Np+1:2*Np); delta_V=delta_n_p_V(2*Np+1:end);


delta_E=-diff([delta_V_L,delta_V,delta_V_R]*VN)./delta_x;
diff_delta_U=diff([delta_V_L,delta_V,delta_V_R]*(VN/V_thermal));
diff_delta_V=diff([delta_V_L,delta_V,delta_V_R]);


n_all=[delta_n_L,delta_n,delta_n_R]; p_all=[delta_p_L,delta_p,delta_p_R];
delta_Jn=N1*q/(delta_x)*Dn.*(-n_all(1:end-1).*B_minusU_diff + n_all(2:end).*B_U_diff);
delta_Jn=delta_Jn+(N1*q*delta_x/tau1)*diff_delta_V.*hn;

delta_Jp=-N2*q/(delta_x)*Dp.*(-p_all(1:end-1).*B_U_diff + p_all(2:end).*B_minusU_diff);
delta_Jp=delta_Jp-(N2*q*delta_x/tau2)*diff_delta_V.*hp;


delta_JD=-Ks*epsilon0*i*omega*delta_E; %assumed exp(-iwt) dependence

delta_J=delta_Jn+delta_Jp+delta_JD;

delta_J1=mean(delta_J(4:end-4));
Y=delta_J1/(delta_VA);
conductance=real(Y); capacitance=-imag(Y)/omega;




