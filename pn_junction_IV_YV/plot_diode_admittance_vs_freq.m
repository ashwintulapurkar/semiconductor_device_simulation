%%% small signal Y plots as a function of freq, for the last VA
omega_array=10.^linspace(-2,2,40)/tau_p;
for iomega=1:length(omega_array)
    omega=omega_array(iomega);
    get_diode_admittance;
    G_array1(iomega)=conductance; C_array1(iomega)=capacitance;
end
figure(23); loglog(omega_array*tau_p,G_array1,'ro' );  xlabel('\omega\tau_p'); ylabel('G (S/m^2)'); hold on; 
title([str1,', V_A=',num2str(VA),' V']); set(gca,'fontsize', 14); 

figure(24); loglog(omega_array*tau_p,C_array1,'ro' ); xlabel('\omega\tau_p'); ylabel('C (F/m^2)'); hold on; 
title([str1,', V_A=',num2str(VA),' V']); set(gca,'fontsize', 14); 
%add admittance plot as a function of freq (for last VA) from diode equations
get_pn_junction_Y_vs_omega_theory;
figure(23); loglog(omega_array*tau_p,G_diffusion+0*G_RG,'r' ); 
legend('G (num)', 'G (theory)');

figure(24); loglog(omega_array*tau_p,C_diffusion+C_junction,'r',omega_array*tau_p,C_diffusion,'b',omega_array*tau_p,C_junction,'g' ); 
legend('C (num)', 'C (theory)', 'C_{diff} (theory)', 'C_J (theory)');
%%% end small signal Y plots as a function of freq, for the last VA
