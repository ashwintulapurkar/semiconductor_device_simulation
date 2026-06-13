% plot various quantities as a function of voltage
str1=['N_A=',num2str(NA*1e-6,'%.0e'), '/cc, N_D=',num2str(ND*1e-6,'%.0e'),' /cc'];

figure(1); plot(VA_array,J_array,'ro','LineWidth',1.,'DisplayName','num'); xlabel('V_A (V)'); ylabel('J (A/m^2)'); hold on;
set(gca,'fontsize', 14);
title(str1);
legend;

figure(2); semilogy(VA_array,abs(J_array),'mo','LineWidth',1.,'DisplayName','num' ); xlabel('V_A (V)'); ylabel('abs J (A/m^2)'); hold on; 
title(str1); set(gca,'fontsize', 14);
legend;


figure(3); plot(VA_array(1:end-1),kT*diff(log(abs(J_array)))./diff(VA_array),'ro','LineWidth',1.,'DisplayName','num' );
xlabel('V_A (V)'); ylabel('(k_BT/q) d ln(I)/dV_A'); hold on;
ylim([0,3]); title(str1); set(gca,'fontsize', 14);
legend;

figure(4); plot(VA_array,Efn_minus_Efp_array,'r' ); xlabel('V_A (V)'); ylabel('(E_{Fn}-E_{Fp})/q (x=0) (V)');
title(str1); set(gca,'fontsize', 14);


figure(5); plot(Efn_minus_Efp_array(1:end-1),kT*diff(log(abs(J_array)))./diff(Efn_minus_Efp_array),'ro' );
xlabel('V=(E_{Fn}-E_{Fp})/q (x=0) (V)'); 
ylabel('(k_BT/q) d ln(I)/dV'); hold on; ylim([0,2]); 
title(str1); set(gca,'fontsize', 14);


% figure(6); semilogy(Efn_minus_Efp_array,abs(J_array),'ro' ); xlabel('(E_{Fn}-E_{Fp})(x=0) (V)'); ylabel('abs J (A/m^2)'); hold on; 
% title(str1); set(gca,'fontsize', 14);

% end plot various quantities as a function of voltage 
