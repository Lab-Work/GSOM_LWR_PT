function y = traffic_flux_trig(rho,rhoc,qm,rhom)

%==========================================================================
% This is the dayanzo-newell flux function, or the triangle form
% fundamental diagrams
% Shimao Fan, Math department, Temple University
% Aug.08 2012.
%==========================================================================

% find the slop
% left slop
sl = qm/rhoc;   % which is positive
% right slop 
sr = qm/(rhoc-rhom); % is negative.
%---------------------------------------
% there is one discontinuity at 'rhoc'
case1 = rho <= rhoc;
case2 = rho > rhoc;
%---------------------------------------
y1 = sl*(rho-rhoc)+qm;
y2 = sr*(rho-rhoc)+qm;

y = y1.*case1+y2.*case2;