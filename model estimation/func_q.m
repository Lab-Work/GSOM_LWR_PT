function q = func_q(a,b,rho0,uf,rhof,rhom)

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% function calculate maximum flux in the 2-parameter family of fd curves. 
% a,b --- free parameter
% v0 --- velovity
% rho0 --- free flow threshold density
% u0 --- corresponding flux
% rhom --- maximum density
% Shimao Fan
% May 28 2013
%//////////////////////////////////////////////////////////////////////////
if nargin<1
a = [1,10]; b = [100,110]; rho0 = 67;
uf = 79; rhof = 1266; 
rhom = 800;
end

% calculate parameter
v0 = diff_fd_free(rho0,uf,rhof);
f0 = fd_free(rho0,uf,rhof);
drho = rhom-rho0;
%--------------------
w0 = func_w(rho0,a,b);
I = int_actan(rhom,a,b)-int_actan(rho0,a,b);
% transform w to v = u'(\rho)
k = (v0*drho+f0)./(w0*drho-I);  % slop
beta = v0-k.*w0;                % y intersection

rhoc = a.*tan(beta./k)+b;

q = func_fd_seibold_2p(rhoc,a,b,rho0,v0,f0,uf,rhof,rhom);
% q = 0*a;
% for i = 1:length(a)
%   q(i) = func_fd_seibold(rhoc(i),a(i),b(i),v0,rho0,f0,rhom);
% end

% rho = [0 rhom];
% 
% plot(rho,rho.*q)

%\\\\\\\\\\\\\\\\\\
% free flow branch
%//////////////////
function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);

function y = diff_fd_free(rho,uf,rhof)

y = uf*(1-2*rho/rhof);

function y = func_w(rho,a,b)

y = -atan((rho-b)./a);

function y = int_actan(x,a,b)

y = -a.*((x-b)./a.*atan((x-b)./a)-0.5*log(1+(x-b).^2./(a.^2)));