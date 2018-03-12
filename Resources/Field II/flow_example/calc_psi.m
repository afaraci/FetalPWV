%   Calculate psi in the Womerslye-Evans model
%
%  Version 1.1, JAJ, June 7, 2011
%  Version 1.2, JAJ, June 7, 2011: Bessel changed to besselj

i=sqrt(-1);
for p=1:8
    disp(['Harmonic number ',num2str(p)])
    omega=p*sim.omega_0;
    tau_alpha=i^(3/2)*sim.R*sqrt(omega*sim.rho/sim.mu);
    Be=tau_alpha*besselj(0,tau_alpha);
    psi(:,p)=(Be-tau_alpha*besselj(0,sim.r_rel*tau_alpha))/(Be-2*besselj(1,tau_alpha));
end