function B1 = Bdash(z)

B11=-0.5+z/6-z.^3/180;
B22=(exp(z)-1-z.*exp(z))./(exp(z)-1).^2; B22(isnan(B22))=-.5; B22(isinf(B22))=-.5; 
B1=B11.*heaviside(1e-4-abs(z))+B22.*heaviside(abs(z)-1e-4);

end

%bernoulli(n,0)*x^n/factorial(n)