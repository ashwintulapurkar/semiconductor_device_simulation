function [V,E, rho, wL,wR] = fun_get_V_E_rho_from_depletion_approx(x_array,K1,K2, N_L,N_R, VA, Vbi)
 
epsilon0=8.85e-12;q=1.6e-19; 
%Vbi=+ve -> pn junction,    Vbi=-ve -> np junction

x_array_L=x_array(x_array<0);
x_array_R=x_array(x_array>=0);

Vbi_dash=Vbi-VA; %VA is applied to the left side. Right side is grounded.


wL=sqrt(abs(2*K1*K2*epsilon0*N_R*Vbi_dash/(q*N_L*(K1*abs(N_L)+K2*abs(N_R))))); %depletion reg width on left side
wR=sqrt(abs(2*K1*K2*epsilon0*N_L*Vbi_dash/(q*N_R*(K1*abs(N_L)+K2*abs(N_R))))); %depletion reg width on right side

%x_array_L=flip(x_array_L);
x=[x_array_L,x_array_R]; 

V1=sign(Vbi)*(q*abs(N_L)/(2*K1*epsilon0))*(wL+x).^2;
V2=Vbi_dash-sign(Vbi)*(q*abs(N_R)/(2*K2*epsilon0))*(wR-x).^2; 
V=heaviside(x+wL).*heaviside(-x).*V1 +heaviside(wR-x).*heaviside(x).*V2 +heaviside(x-wR).*Vbi_dash;


E1=-sign(Vbi)*(q*abs(N_L)/(K1*epsilon0))*(wL+x_array_L);
E2=-sign(Vbi)*(q*abs(N_R)/(K2*epsilon0))*(wR-x_array_R);
E=[heaviside(x_array_L+wL).*ceil(heaviside(-x_array_L)).*E1, heaviside(wR-x_array_R).*ceil(heaviside(x_array_R)).*E2];
%E=heaviside(x+wL).*ceil(heaviside(-x)).*E1 +heaviside(wR-x).*floor(heaviside(x)).*E2;

rho=[heaviside(x_array_L+wL).*ceil(heaviside(-x_array_L))*N_L, heaviside(wR-x_array_R).*ceil(heaviside(x_array_R))*N_R];
%rho=heaviside(x+wL).*heaviside(-x).*N_L +heaviside(wR-x).*heaviside(x).*N_R;

%rho=from ionized donors and acceptors
end
