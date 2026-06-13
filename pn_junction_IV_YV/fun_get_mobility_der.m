function mu_dash = fun_get_mobility_der(mu0, beta, v_sat, E)
E1=abs(E);
mu=mu0./((1+(mu0*E1/v_sat).^beta).^(1/beta));
mu_dash=-(mu*mu0/v_sat).*(mu0*E1/v_sat).^(beta-1)./(1+(mu0*E1/v_sat).^beta);
mu_dash=mu_dash.*sign(E);

%mu_dash=0; %if mobility is independent of E
end
