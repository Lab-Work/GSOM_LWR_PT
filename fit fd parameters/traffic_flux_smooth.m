function f = traffic_flux_smooth(x,lambda,p,alpha,rom)

%==========================================================================
 % define a smooth fundamental diagram function with 3 free parameters.
 % lambda -- is a small number
 % P -- control the critical point, p<.5, on the left, right otherwise
 % alpha -- is the amplitude 
 % rhom -- is the maximum density which only depends on the # of lanes and
 % road conditions.
%==========================================================================

a = sqrt(1+(p.*lambda).^2);
b = sqrt(1+((1-p).*lambda).^2);
y = ((x/rom)-p).*lambda;
f = alpha.*(a+(b-a).*(x/rom)-sqrt(1+y.^2));

