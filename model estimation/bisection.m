function p = bisection(f,v,tol)
%   Runs the bisection rootfinding algorithm for the function
%   f on the starting interval v=[a,b].
%
%   (C) 2010/09/21 by Benjamin Seibold

% fv = f(v); % values at interval boundaries
% f(v(1))
% f(v(2))
fv(1) = f(v(1));
fv(2) = f(v(2));
% prod(fv)
if prod(fv)>=0, p = v(end); end

k = 1;
nl = 20;      % maximum number of loop
while 1  && k < nl
% while 1
    p = sum(v)/2; fp = f(p); % position and value of midpoint
    if diff(v)<tol&&abs(fp)<tol||fp==0, break, end % stopping criterion
    i = (sign(fp)~=sign(fv(1)))+1; % index to be replaced by midpoint
    v(i) = p; fv(i) = fp; % cut interval in half
    k = k+1;
end
