function f = traffic_flux_piecewise(rho,rhoc,rhom,vm)


qm = rhoc.*greenshield(rhoc,vm,rhom);
% qm = vm*rhoc;

case1 = rho<=rhoc;          % free flow regime
case2 = rho>rhoc;           % cogested flow regime
% free branch
uf = vm.*rho.*(1-rho./rhom);     % free flow velocity branch
% congested branch
rhom2 = rhoc*5;   % fix the rate rhoc/rhom
ug = (qm./(rhoc-rhom2)).*(rho-rhom2);

f = case1.*uf+case2.*ug;    % piecewise continuous velocity function



%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////

function y = greenshield(rho,um,rhom)

y = um.*(1-rho./rhom);