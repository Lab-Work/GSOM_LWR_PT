function q = func_inverse_v_gpt(rho,u,q_min,q_max,rhoc,vm,rhom)

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
flux = u.*rho;   % flux
v = [q_min,q_max]; % bound for bisection search
f = @(q) flux - func_flux_creen_pt(rho,q,rhoc,vm,rhom);

%\\\\\\\\\\\\\\\\\\\\\
% perform projection
%/////////////////////

case1 = f(q_min)<0;   % data below lower bound
case2 = f(q_max)>0;   % data above upper bound


if case1    % data lay below lower bound
    q = q_min;
else if case2
        q = q_max;
    else
       q = bisection(f,v,tol);
    end
end

    
    % make figure
%     plot(wvec,f(wvec),'k--'), hold on
%     plot(u(i),f(u(i)),'bo'), hold on
%     plot(w(i),f(w(i)),'r*','markersize',12), hold off
%     axis([w_min w_max -20 20])
%     drawnow
% end