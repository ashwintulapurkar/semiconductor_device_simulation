function mu = fun_get_mobility(mu0, beta, v_sat, E)
E=abs(E);
mu=mu0./((1+(mu0*E/v_sat).^beta).^(1/beta));
%mu=mu0;

end
