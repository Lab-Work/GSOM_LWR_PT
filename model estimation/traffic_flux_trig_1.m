function y = traffic_flux_trig_1(rho,sl,rhoc,rhom)

%==========================================================================
% This is the dayanzo-newell flux function, or the triangle form
% fundamental diagrams
% Shimao Fan, Math department, Temple University
% Aug.08 2012.
%==========================================================================

% find the slop
% left slop
% sl = qm/rhoc;   % which is positive
qm = sl*rhoc;
% right slop 
sr = qm/(rhoc-rhom); % is negative.
%---------------------------------------
% there is one discontinuity at 'rhoc'
case1 = rho <= rhoc;
case2 = rho > rhoc;
%---------------------------------------
y1 = sl*rho;
% y2 = sr*(rho-rhoc)+qm;
y2 = sr*(rho-rhom);

y = y1.*case1+y2.*case2;