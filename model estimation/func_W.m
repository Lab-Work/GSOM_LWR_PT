 function w = func_W(lambda,p,alpha,rhom)
%======================================================================
% w is the V(rho = 0), i.e., the empty road velocity. Therefore, the 
% quantity w can be represented as the function below
% Shimao Fan, Temple University
% Jan 24 2013 (last modified)
%======================================================================
a = sqrt(1+(p.*lambda).^2);
b = sqrt(1+((1-p).*lambda).^2);

w = alpha.*(b-a+lambda.^2.*p./a)/rhom;
