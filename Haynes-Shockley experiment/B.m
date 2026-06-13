function B1 = B(z)
% Bernoulli function
%if z==0, B1=1;

% if abs(z)<1e-4, B1=1-0.5*z+z.^2/12-z.^4/720;
% else
% B1=z./(exp(z)-1);
% end % if statements in function --> z can not be array
B11=1-0.5*z+z.^2/12-z.^4/720;
B22=z./(exp(z)-1); B22(isnan(B22))=1; B22(isinf(B22))=1; %z=0 gives NaN, small values give inf
B1=B11.*heaviside(1e-4-abs(z))+B22.*heaviside(abs(z)-1e-4);
end

%bernoulli(n,0)*x^n/factorial(n)