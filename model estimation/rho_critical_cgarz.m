function y = rho_critical_cgarz(a,b,rho0,v0,f0,rhom)

% find critical density of the collapsed model
% shimao fan
% march 24 2014
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

drho = rhom-rho0;

w0 = w(rho0,a,b);
%-------------------------------------------
% transform w to v
I = int_actan(rhom,a,b)-int_actan(rho0,a,b);
k = (v0*drho+f0)./(w0*drho-I); % slop
beta = v0-k.*w0;   % y intersection
%-------------------------------------------
y = a.*tan(beta./k) + b;


% plot(rho, f,'linewidth',3)
% axis([0 rhom 0 10000])
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%////////////////////////////////////////////////////////////
function y = w(rho,a,b)
y = -atan((rho-b)./a);

%--------------------------------------------
function y = int_actan(x,a,b)

y = -a.*((x-b)./a.*atan((x-b)./a)-0.5*log(1+(x-b).^2./(a.^2)));