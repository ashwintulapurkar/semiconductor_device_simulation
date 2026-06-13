% add depletion approx plots

[V_depl,E_depl, rho_depl, wL,wR] = fun_get_V_E_rho_from_depletion_approx(x_array_L,x_array_R,Ks_L,Ks_R, ND_L-NA_L, ND_R-NA_R, 0, Vbi_non_deg);

figure(1); plot(x_array*1e6,V_depl,'b--','LineWidth',1.,'DisplayName','depl');
%legend('num', 'depl approx'); 
xline(-wL*1e6,'HandleVisibility','off'); xline(wR*1e6,'HandleVisibility','off');


figure(2); plot(x_array*1e6,E_depl,'b--','LineWidth',1.,'DisplayName','depl');
xline(-wL*1e6,'HandleVisibility','off'); xline(wR*1e6,'HandleVisibility','off');

figure(3); plot(x_array*1e6,rho_depl*1e-6,'b--','LineWidth',1.,'DisplayName','depl');
xline(-wL*1e6,'HandleVisibility','off'); xline(wR*1e6,'HandleVisibility','off'); 


% Ec_depl=(Ec_L-V_depl).*(1-heaviside(x_array))+ (Ec_R-V_depl).*(heaviside(x_array));
% Ev_depl=(Ev_L-V_depl).*(1-heaviside(x_array))+ (Ev_R-V_depl).*(heaviside(x_array));
% Ei_depl=(Ei_L-V_depl).*(1-heaviside(x_array))+ (Ei_R-V_depl).*(heaviside(x_array));

Ec_depl=[Ec_L-V_depl(1:length(x_array_L)), Ec_R-V_depl(length(x_array_L)+1:end)];
Ev_depl=[Ev_L-V_depl(1:length(x_array_L)), Ev_R-V_depl(length(x_array_L)+1:end)];
Ei_depl=[Ei_L-V_depl(1:length(x_array_L)), Ei_R-V_depl(length(x_array_L)+1:end)];


figure(7) ; plot(x_array*1e6,Ec_depl,'r--',x_array*1e6,Ev_depl,'b--',x_array*1e6,Ei_depl,'g--');
xline(-wL*1e6,'HandleVisibility','off'); xline(wR*1e6,'HandleVisibility','off');
title('-nondeg,  -- depl');

