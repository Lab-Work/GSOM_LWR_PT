function f = traffic_flux_fitted(rho,w,c_lambda,c_p,c_alpha,rhom)% smooth flux function

%///////////////////////////////////////////////////////////////////////////
% this function is the smooth traffic flux function, and perform data-fitting process. 
% This function suppose three free parameter as function depends on w,
% i.e., alpha(w), p(w), and lambda(w).
% Shimao Fan
% Last modified on Jan 31 2013.
%///////////////////////////////////////////////////////////////////////////
p = polyval(c_p,w);    % p parameter
lambda = polyval(c_lambda,w);  % \lambda parameter
alpha = polyval(c_alpha,w);     % alpha parameter


% only work for linear case
% lambda = c_lambda(1)*w+c_lambda(2);
% p = c_p(1)*w+c_p(2);
% alpha = c_alpha(1)*w+c_alpha(2);      
      
a = sqrt(1+(p.*lambda).^2);
b = sqrt(1+((1-p).*lambda).^2);
% whos
y = ((rho./rhom)-p).*lambda;

f = alpha.*(a+(b-a).*(rho./rhom)-sqrt(1+y.^2)); 