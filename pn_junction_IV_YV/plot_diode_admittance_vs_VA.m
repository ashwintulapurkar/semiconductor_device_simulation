%%% small signal Y plots as a function of VA (for a fixed freq)
str1=['N_A=',num2str(NA*1e-6,'%.0e'), '/cc, N_D=',num2str(ND*1e-6,'%.0e'),' /cc'];
str2=[', \omega\tau_p=',num2str(omega*tau_p)];

figure(21); semilogy(VA_array,G_array,'ro' ); xlabel('V_A (V)'); ylabel('G (S/m^2)'); hold on; 
title([str1,str2]);  set(gca,'fontsize', 14); 

figure(22); semilogy(VA_array,C_array,'ro' ); xlabel('V_A (V)'); ylabel('C (F/m^2)'); hold on; 
title([str1,str2]); set(gca,'fontsize', 14); 

%add admittance plot as a function of VA (for a fixed freq) from diode equations
get_pn_junction_Y_vs_VA_theory;
figure(21); semilogy(VA_array,G_diffusion+0*G_RG,'r' ); 
legend('G (num)', 'G diode eqn');  set(gca,'fontsize', 14); 
figure(22); semilogy(VA_array,C_diffusion+C_junction,'r',VA_array,C_diffusion,'b',VA_array,C_junction,'g' ); 
set(gca,'fontsize', 14); legend('C (num)', 'C (theory)', 'C_{diff}(theory)', 'C_J (theory)');

%%% end small signal Y plots as a function of VA (for a fixed freq)
