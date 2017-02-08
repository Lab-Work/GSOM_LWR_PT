function q = func_inverse_v_pase_trans(rho,u,q_min,q_max,u_r,rho_r,rhom)

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% For given rho and u, this function try to identify which curve this point
% belongs to. I would like to solve the problem in a flow-density plane,
% that is (rho, u*rho), and find the corresponding 'w' by the Newton's
% Iteration.
% Shimao Fan
% Feb 15 2012
% Temple University
% for phase transition model
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
tol = 1.0e-02;     % creteria for bisection search
v = [q_min,q_max]; % bound for bisection search
f = @(q) u-func_u_phase_transition(rho,q,u_r,rho_r,rhom);

q = bisection(f,v,tol);
