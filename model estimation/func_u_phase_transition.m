function u = func_u_phase_transition(rho,q,u_r,rho_r,rhom)

%--------------------------------------------------------
% this is the velocity function of phase transition model
% shimao fan
% April 10 2013
%--------------------------------------------------------
% quadratic
rhoc = rhoc_q(q,u_r,rho_r);  % find out critical density, given q
% linear form
% rhoc = q/u_f;   
% k = u_r*(1-2*rhoc/rho_r); % Q'(rho_critical)
% k = 0;
case1 = rho<=rhoc;          % free flow regime
case2 = rho>=rhoc;           % cogested flow regime

% free branch
uf = u_r*(1-rho/rho_r);     % free flow velocity branch
% congested branch
ug = (q./(rhoc-rhom)).*(rho-rhom)./rho;
% congested flow velocity branch
% ug = -q./(rhom-rhoc).^2.*(rho+rhoc.^2./rho-2*rhoc) + q./rho;
% ug = -(q+k.*(rhom-rhoc))./(rhom-rhoc).^2.*(rho+rhoc.^2./rho-2*rhoc)...
%     +k.*(rho-rhoc)./rho+ q./rho;

% ug = -(q+k.*(rhom-rhoc))./(rhom-rhoc).^2.*(rho+rhoc.^2./rho-2*rhoc)...
%     +k.*(1-rhoc./rho)+ q./rho;

u = case1.*uf+case2.*ug;    % piecewise continuous velocity function



%\\\\\\\\\\\\\\\\\\\\\\\\\\\
% given q, idenfity rhoc
%///////////////////////////

function y = rhoc_q(q,u_r,rho_r)

y = (rho_r/2)*(1-sqrt(1-4*q/(rho_r*u_r)));
% y = (rho_r/2)*(1+sqrt(1-4*q/(rho_r*u_r)));