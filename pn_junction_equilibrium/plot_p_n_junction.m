close all;

plot_non_deg_soln;
plot_deg_soln;
plot_depl_soln;

figure(1); legend; title(['N_A=',num2str((NA_L-ND_L)*1e-6),'/cc, N_D=',num2str((ND_R-NA_R)*1e-6), ' /cc']);
figure(2); legend; title(['N_A=',num2str((NA_L-ND_L)*1e-6),'/cc, N_D=',num2str((ND_R-NA_R)*1e-6), ' /cc']);
figure(3); legend; title(['N_A=',num2str((NA_L-ND_L)*1e-6),'/cc, N_D=',num2str((ND_R-NA_R)*1e-6), ' /cc']);
figure(4); legend; title(['N_A=',num2str((NA_L-ND_L)*1e-6),'/cc, N_D=',num2str((ND_R-NA_R)*1e-6), ' /cc']);
figure(5); legend; title(['N_A=',num2str((NA_L-ND_L)*1e-6),'/cc, N_D=',num2str((ND_R-NA_R)*1e-6), ' /cc']);

