%%% small signal Y plots for the last VA and one freq. To see that real and imaginary parts of current are conserved.

get_diode_admittance;
str2=['\omega\tau_p=',num2str(omega*tau_p)];
figure(31); plot(x_array_E,real(delta_Jn),'r',x_array_E,real(delta_Jp),'b',x_array_E,real(delta_JD),'g',x_array_E,real(delta_J),'k' );  
set(gca,'fontsize', 14); 
xlabel('x (m)'); ylabel('Re(\deltaJ) (A/m^2)'); xline(-xp); xline(xn); xline(0);
legend('Re(\deltaJ_n)','Re(\deltaJ_p)','Re(\deltaJ_D)','Re(\deltaJ)');
title([str1,', ', str2, ', V_A=',num2str(VA),' V']);

figure(32); plot(x_array_E,imag(delta_Jn),'r',x_array_E,imag(delta_Jp),'b',x_array_E,imag(delta_JD),'g',x_array_E,imag(delta_J),'k' ); 
set(gca,'fontsize', 14); 
xlabel('x (m)'); ylabel('Im(\deltaJ) (A/m^2)'); xline(-xp); xline(xn); xline(0);
legend('Im(\deltaJ_n)','Im(\deltaJ_p)','Im(\deltaJ_D)','Im(\deltaJ)');
title([str1,', ',str2, ', V_A=',num2str(VA),' V']);

return
figure(33); plot(x_array_E,abs(delta_Jn),'r',x_array_E,abs(delta_Jp),'b',x_array_E,abs(delta_JD),'g',x_array_E,abs(delta_J),'k' ); set(gca,'fontsize', 14); 
xlabel('x (m)'); ylabel('abs \delta J (A/m^2)'); xline(-xp); xline(xn); xline(0);
legend('abs(\delta J_n)','abs(\delta J_p)','abs(\delta J_D)','abs(\delta J)');
title(['\omega\tau_p=',num2str(omega*tau_p), ', V_A=',num2str(VA),' V']);

figure(34); plot(x_array_E,real(delta_J),'r',x_array_E,imag(delta_J),'b',x_array_E,abs(delta_J),'g' ); set(gca,'fontsize', 14); 
xlabel('x (m)'); ylabel('abs \delta J (A/m^2)'); xline(-xp); xline(xn); xline(0);
legend('real(\delta J)','imag(\delta J)','abs(\delta J)');
title(['\omega\tau_p=',num2str(omega*tau_p), ', V_A=',num2str(VA),' V']);

return
figure(35); plot(X_array,real([0,delta_n,0])*1e-6,'r',X_array,real([0,delta_p,0])*1e-6,'b' ); set(gca,'fontsize', 14); 
xlabel('x (m)'); ylabel('(/cm^3)'); xline(-xp); xline(xn); xline(0);
legend('real(\delta n)','real(\delta p)');
title(['\omega\tau_p=',num2str(omega*tau_p), ', V_A=',num2str(VA),' V']);

figure(36); plot(X_array,imag([0,delta_n,0])*1e-6,'r',X_array,imag([0,delta_p,0])*1e-6,'b' ); set(gca,'fontsize', 14); 
xlabel('x (m)'); ylabel('(/cm^3)'); xline(-xp); xline(xn); xline(0);
legend('imag(\delta n)','imag(\delta p)');
title(['\omega\tau_p=',num2str(omega*tau_p), ', V_A=',num2str(VA),' V']);

%%% end small signal Y plots for the last VA and one freq
