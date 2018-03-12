%   Calculate psi

i=sqrt(-1);
for p=1:8
  disp(['Harmonic number ',num2str(p)])
  omega=p*omega_0;
  tau_alpha=i^(3/2)*R*sqrt(omega*rho/mu);
  Be=tau_alpha*bessel(0,tau_alpha);
  psi(:,p)=(Be-tau_alpha*bessel(0,r_rel*tau_alpha))/(Be-2*bessel(1,tau_alpha));
  end



