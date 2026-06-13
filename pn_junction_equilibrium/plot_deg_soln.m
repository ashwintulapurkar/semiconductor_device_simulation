%plot deg soln

figure(1); plot(x_array*1e6,phi_array_deg,'r','LineWidth',1.,'DisplayName','deg'); 
xlabel('x (\mu m)'); ylabel('pot (V)'); hold on ;
set(gca,'fontsize', 14); xline(0,'HandleVisibility','off');  

figure(2); plot(x_array_for_E*1e6,E_array_deg,'r','LineWidth',1.,'DisplayName','deg');
xlabel('x (\mu m)'); ylabel('Ele field (V/m)');  hold on ;
set(gca,'fontsize', 14); xline(0,'HandleVisibility','off');   

figure(3); plot(x_array*1e6,rho_array_deg*1e-6,'r','LineWidth',1.,'DisplayName','deg');
xlabel('x (\mu m)'); ylabel('(\rho/q)/cm^3'); hold on ;
set(gca,'fontsize', 14); xline(0,'HandleVisibility','off');

figure(4); plot(x_array*1e6,n_array_deg*1e-6,'r','LineWidth',1.,'DisplayName','deg'); 
xlabel('x (\mu m)'); ylabel('n /cm^3'); hold on ;
set(gca,'fontsize', 14); xline(0,'HandleVisibility','off');  

figure(5); plot(x_array*1e6,p_array_deg*1e-6,'r','LineWidth',1.,'DisplayName','deg'); 
xlabel('x (\mu m)'); ylabel('p /cm^3'); hold on ;
set(gca,'fontsize', 14); xline(0,'HandleVisibility','off');   


% % draw band diagram
Ec=(Ec_L-phi_array_deg).*(1-heaviside(x_array))+ (Ec_R-phi_array_deg).*(heaviside(x_array));
Ev=(Ev_L-phi_array_deg).*(1-heaviside(x_array))+ (Ev_R-phi_array_deg).*(heaviside(x_array));
Ei=(Ei_L-phi_array_deg).*(1-heaviside(x_array))+ (Ei_R-phi_array_deg).*(heaviside(x_array));
EF=EF_L_deg*ones(1,length(x_array));
figure(7) ; plot(x_array*1e6,Ec,'r',x_array*1e6,Ev,'b',x_array*1e6,Ei,'g',x_array*1e6,EF,'m','LineWidth',1.); 
xlabel('x (\mu m)'); ylabel('Ene(eV)');  
set(gca,'fontsize', 14);  hold on ;  
xline(0,'HandleVisibility','off');
title('-nondeg');

