function f = func_fd_seibold_2p(rho,a,b,rho0,v0,f0,uf,rhof,rhom)

% here, a -- controal the curvature
%       b -- maximum point, or critical density
% [u0,v0,rho0] are initial flux, velocity, and density, which are fd
% parameters.
% Idea of Dr. Benjamin Seibold,

% May 26, 2013
% Shimao Fan, Temple University
%---------------------------------------------------------


%if nargin<1
%a = 10; b = 100; rho0 = 67;
%uf = 79; rhof = 1266; 
%rhom = 800;
%rho = 0:rhom;
%end


% define free branch
% uf = 79; rhof = 1266; 
%v0 = diff_fd_free(rho0,uf,rhof);
%f0 = fd_free(rho0,uf,rhof);
drho = rhom-rho0;

%w = @(rho) -atan((rho-b)./a);
w0 = w(rho0,a,b);
%-------------------------------------------
% transform w to v
I = int_actan(rhom,a,b)-int_actan(rho0,a,b);
k = (v0*drho+f0)./(w0*drho-I); % slop
beta = v0-k.*w0;   % y intersection
%-------------------------------------------
u1 = @(rho) k.*int_actan(rho,a,b)+beta.*rho;
c = f0-u1(rho0);   %value of constant term
f1 = u1(rho)+c;  
%-------------------------------------------
case1 = rho<=rho0;
case2 = rho>rho0;
%-------------------------------------------
f = case1.*fd_free(rho,uf,rhof)+case2.*f1;

%\\\\\\\\\\\\\
% tan(inters/k)
%/////////////
% figure(1)

% plot(rho, f,'linewidth',3)
% axis([0 rhom 0 10000])
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%////////////////////////////////////////////////////////////
function y = w(rho,a,b)
y = -atan((rho-b)./a);
%--------------------------------------------
function y = int_actan(x,a,b)

y = -a.*((x-b)./a.*atan((x-b)./a)-0.5*log(1+(x-b).^2./(a.^2)));
%-------------------------------------------
% free flow branch
function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);
%-------------------------------------------
%function y = diff_fd_free(rho,uf,rhof)

%y = uf*(1-2*rho/rhof);
