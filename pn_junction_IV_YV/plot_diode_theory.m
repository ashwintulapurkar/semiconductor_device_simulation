%add depletion approx plot
get_n_p_J_from_theory;

figure(1); plot(VA_array,J_array1,'b','DisplayName','J (ideal diode)');
figure(1); plot(VA_array,J_array1+J_RG_array,'g','DisplayName','J (including RG)' ); 
legend;

figure(2); semilogy(VA_array,abs(J_array1),'b','DisplayName','J (ideal diode)');
figure(2); semilogy(VA_array,abs(J_array1+J_RG_array),'g','DisplayName','J with RG' ); 
legend;

figure(11); plot(X_array,Jn_depl_approx,'r--','DisplayName','J_n depl');
figure(11); plot(X_array,Jp_depl_approx,'b--','DisplayName','J_p depl');
figure(11); plot(X_array,J_depl_approx,'g--','DisplayName','J_{total} depl' );  
xline(-xp,'HandleVisibility','off'); xline(xn,'HandleVisibility','off');
legend;

figure(12); xline(-xp,'HandleVisibility','off'); xline(xn,'HandleVisibility','off');
figure(13); xline(-xp,'HandleVisibility','off'); xline(xn,'HandleVisibility','off');
figure(14); xline(-xp); xline(xn);
figure(15); xline(-xp); xline(xn);

figure(16); semilogy(X_array,n_depl*1e-6,'r','DisplayName','n depl');
figure(16); semilogy(X_array,p_depl*1e-6,'b','DisplayName','p depl' );
xline(-xp,'HandleVisibility','off'); xline(xn,'HandleVisibility','off');
legend;

figure(17); plot(X_array,Ec_depl,'r--','DisplayName','E_c');
figure(17); plot(X_array,Ev_depl,'b--','DisplayName','E_v');
figure(17); plot(X_array,FN_depl,'g--','DisplayName','EF_n');
figure(17); plot(X_array,FP_depl,'m--','DisplayName','EF_p' ); 
xline(-xp,'HandleVisibility','off'); xline(xn,'HandleVisibility','off');
legend
%end add depletion approx plot


