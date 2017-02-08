function y = traffic_flux_qrd(rho,vm,rhom)


%==================================================
% The quadratic form flux function
% Shimao Fan
%==================================================


y = vm*rho.*(1-rho/rhom);