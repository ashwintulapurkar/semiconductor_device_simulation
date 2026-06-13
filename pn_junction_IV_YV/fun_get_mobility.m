function mu = fun_get_mobility(mu0, beta, v_sat, E)
E=abs(E);
mu=mu0./((1+(mu0*E/v_sat).^beta).^(1/beta)); %eqn 3.33 in Pierret (Semicond. device fundamentals)
%mu=mu0; %to turn off E dependent mobility

end
