% plot various quantities for the last voltage point
str1=['N_A=',num2str(NA*1e-6,'%.0e'), '/cc, N_D=',num2str(ND*1e-6,'%.0e'),' /cc'];
str2=[str1,', V_A=',num2str(VA),' V'];
figure(11); plot(x_array_E,Jn,'r','LineWidth',1.,'DisplayName','J_n'); hold on;
figure(11); plot(x_array_E,Jp,'b','LineWidth',1.,'DisplayName','J_p');
figure(11); plot(x_array_E,J_total,'g','LineWidth',1.,'DisplayName','J_{total}' );  
xlabel('x (m)'); ylabel('J (A/m^2)');  xline(0,'HandleVisibility','off');  
title(str2); set(gca,'fontsize', 14);
legend;

figure(12); plot(x_array_E,Jn_drift,'r','LineWidth',1.,'DisplayName','J_{n,drift}'); hold on;
figure(12); plot(x_array_E,Jn_diff,'b','DisplayName','J_{n,diff}');
figure(12); plot(x_array_E,Jn,'g','DisplayName','J_{n}' );  xlabel('x (m)'); ylabel('J (A/m^2)');  
xline(0,'HandleVisibility','off'); 
title(str2); set(gca,'fontsize', 14);
legend;


figure(13); plot(x_array_E,Jp_drift,'r','DisplayName','J_{p,drift}'); hold on;
figure(13); plot(x_array_E,Jp_diff,'b','DisplayName','J_{p,diff}');
figure(13); plot(x_array_E,Jp,'g','DisplayName','J_{p}' );  xlabel('x (m)'); ylabel('J (A/m^2)');  xline(0,'HandleVisibility','off');  
hold on; title(str2);  set(gca,'fontsize', 14);
legend;

figure(14); plot(X_array,[V_L,V,V_R]*VN,'r' );  title('V (volts)'); xlabel('x (m)'); ylabel('V (V)'); xline(0); hold on;
title(str2); set(gca,'fontsize', 14);

figure(15); plot(x_array_E,E,'r' );  title('E (V/m)'); xlabel('x (m)'); ylabel('E (V/m)'); xline(0);  hold on;
title(str2); set(gca,'fontsize', 14);

figure(16); semilogy(X_array,N1*1e-6*[n_L,n,n_R],'ro',X_array,N2*1e-6*[p_L,p,p_R],'bo' ); xline(0); legend('n','p'); 
 hold on; xlabel('x (m)'); ylabel('/cc'); title([str1,', V_A=',num2str(VA),' V']); set(gca,'fontsize', 14);

%draw band diagram for last value of VA
Ec_array=Ec-[V_L,V,V_R]*VN; Ev_array=Ev-[V_L,V,V_R]*VN; Ei_array=Ei-[V_L,V,V_R]*VN;
EFn_array=Ei_array+V_thermal*log(N1*[n_L,n,n_R]/ni); EFp_array=Ei_array-V_thermal*log(N2*[p_L,p,p_R]/ni);
figure(17); plot(X_array,Ec_array,'r', X_array,Ev_array,'b',X_array,EFn_array,'g',X_array,EFp_array,'m' ); set(gca,'fontsize', 14);
xlabel('x (m)'); ylabel('ene (eV)'); xline(0);  hold on; title([str1,', V_A=',num2str(VA),' V']);
legend('E_c','E_v','Ef_n','Ef_p');


% end plot various quantities for the last voltage point
