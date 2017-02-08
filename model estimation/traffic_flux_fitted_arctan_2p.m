function f = traffic_flux_fitted_arctan_2p(rho,w,ca,cb,rho0,v0,f0,uf,rhof,rhom)% smooth flux function
%==========================================================================
% w(x) = -c arctan((x-b)/a)
% v(x) = w(x) + beta;
% f(x) = \int(v(s))ds
% shimao fan
% may 28 2013
%==========================================================================
% general work
a = polyval(ca,w); b = polyval(cb,w);  % input data
% more efficient
% a = ca(1)*w+ca(2);
% b = cb(1)*w+cb(2);
f = func_fd_seibold_2p(rho,a,b,rho0,v0,f0,uf,rhof,rhom);