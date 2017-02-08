function vel = func_vel_creen_pt(rho,q,rhoc,vm,rhom)
%====================================================
% velocity function of greenshield phase transition model
% rhoc is the critical density, determined by the equilibrum curve
% q == perturbation, which is conserved quantity
% vm, rhom the maximum velocity and density
%====================================================
% close all;

% rhom = 800; vm = 80;
% rho = 0:rhom;
% rhoc = rhom/4;
a = -rhoc*vm./((2*rhoc-3*rhom).^2);
%-----------------------
% n = 15;
% q = linspace(-.5,0.5,n);
% for i = 1:n
%    rhoc = vm.*(1+q(i))./(1+q(i)-vm) 
sigma = func_rhoc(rhoc,a,q,vm,rhom); % calculate new critical density
case1 = rho<sigma;
case2 = rho>=sigma; 
vel = vf(rho,vm).*case1 + vc(rho,vm,rhoc,a,rhom,q).*case2;
% 
% 
% %\\\\\\\\\\\\
% % plotting
% %///////////
% 
% plot(rho,rho.*vel,'k-','linewidth',1.5), hold on
% end

% axis([0 rhom 0 1.25*vm])
%=-===========================================
function y = vf(rho,vm)
  y = vm + 0*rho;
  
    function y = vc(rho,vm,rhoc,a,rhom,q)
        y = (1-rhom./rho).*(a.*(rho-rhoc) + rhoc*vm./(rhoc-rhom)).*(1+q);
        function y = func_rhoc(rhoc,a,q,vm,rhom)
%             y = rhoc-rhoc.*vm.*q./(a.*(rhoc-rhom).*(1+q));
             ca = a;
             cb = -a.*(rhom+rhoc) + rhoc.*vm./(rhoc-rhom)-vm./(1+q);
             cc = a.*rhom.*rhoc -rhoc.*vm.*rhom./(rhoc-rhom);
             
             y = (-cb - sqrt(cb.^2-4*ca.*cc))./(2*ca);