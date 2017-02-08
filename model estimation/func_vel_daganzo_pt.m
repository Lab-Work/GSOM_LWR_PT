function vel = func_vel_daganzo_pt(rho,q,vm,rhoc,rhom)

% velocity function of greenshield phase transition model
% rhoc is the critical density
%====================================================
% close all;

% rhom = 800; vm = 100;
% rho = 0:rhom;
% rhoc = rhom/4;
% case1 = rho<rhoc;
% case2 = rho>=rhoc;
% %-----------------------
% n = 15;
% q = linspace(-.5,0.5,n);
% for i = 1:n
%    rhoc = vm.*(1+q(i))./(1+q(i)-vm) 
%----------------------------
 sigma = func_rhoc(rhoc,q,rhom);
 case1 = rho<=sigma;
 case2 = rho>sigma; 
 vel = vf(rho,vm).*case1 + vc(rho,vm,rhoc,rhom,q).*case2;


%\\\\\\\\\\\\
% plotting
%///////////

% plot(rho,rho.*vel,'k-','linewidth',1.5), hold on
% end

% axis([0 rhom 0 1.25*vm])
%=-===========================================
function y = vf(rho,vm)
  y = vm + 0*rho;
  
    function y = vc(rho,vm,rhoc,rhom,q)
        y = vm.*rhoc./(rhom-rhoc).*(rhom./rho -1).*(1+q);
        function y = func_rhoc(rhoc,q,rhom)
            y = rhom.*(1+q).*rhoc./(rhom+rhoc.*q);