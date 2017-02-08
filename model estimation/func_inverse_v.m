function w = func_inverse_v(rho,u,w_min,w_max,cl,cp,ca,rhom)

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% For given rho and u, this function try to identify which curve this point
% belongs to. I would like to solve the problem in a flow-density plane,
% that is (rho, u*rho), and find the corresponding 'w' by the Newton's
% Iteration.
% Shimao Fan
% Feb 15 2012
% Temple University
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

tol = 1.0e-02;     % creteria for bisection search
q = u.*rho;   % flux
v = [max(w_min,u),w_max]; % bound for bisection search
f = @(w) q - traffic_flux_fitted(rho,w,cl,cp,ca,rhom);

%\\\\\\\\\\\\\\\\\\\\\
% perform projection
%/////////////////////

case1 = f(w_min)<0;   % data below lower bound
case2 = f(w_max)>0;   % data above upper bound


if case1    % data lay below lower bound
    w = w_min;
else if case2
        w = w_max;
    else
       w = bisection(f,v,tol);
    end
end

    
    % make figure
%     plot(wvec,f(wvec),'k--'), hold on
%     plot(u(i),f(u(i)),'bo'), hold on
%     plot(w(i),f(w(i)),'r*','markersize',12), hold off
%     axis([w_min w_max -20 20])
%     drawnow
% end