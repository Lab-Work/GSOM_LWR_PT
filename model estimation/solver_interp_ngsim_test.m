function [num_rho,num_vel,erho,evel] = solver_interp_ngsim_test(tfinal,D,V,data_rho,data_vel,x,x0,vm,rhom,ref)

%--------------------------------------------------------------------------
% This is interpolation (linear) process, which is a test case. This is for
% NGSIM trajectory data, note that error are measured over all the study
% area, thus the error is measured in the sense of integration (in space).
% In the Minneapolis validation, the linear interpolation give decent
% results, which is close to the performance of LWR, ARZQ, and better than
% the LWRQ model.
% In this test, we take four piece of input information, that is the:
% density and velocity on right boundary (rho_r,vel_r)
% density and velocity on left boundary  (rho_l,vel_l)
% and then predict the solution along the study area as a linear
% interpolation.
% Shimao Fan.
% Feb 14 2013.
% Temple University
%--------------------------------------------------------------------------
%\\\\\\\\\\\\\\\\\\\\\
% determine dx and dt
%/////////////////////
% tfinal = 14;         % mins
N = length(x);
midd = round(N/2);
h = x(2)-x(1);
k = 1/36000;         % length of time step, does not work start from ref = 8
t0 = 0:k:(tfinal/60);
m1 = length(t0);
% %-------------------------------------------------
dt = 2*k/ref; 

t = 0:dt:(tfinal/60);

r0 = 16;
sg = .00025;      % 25 cm
gh_cell = r0*sg;  % length of ghost cell is 4 meters

BD(:,1) = data_rho(:,1);
BD(:,2) = data_rho(:,end);
BV(:,1) = data_vel(:,1);
BV(:,2) = data_vel(:,end);

D_lbc = spline(t0,BD(1:m1,1),t);
D_rbc = spline(t0,BD(1:m1,2),t);
V_lbc = spline(t0,BV(1:m1,1),t);
V_rbc = spline(t0,BV(1:m1,2),t);

% interpol respect to time for all the data
for i = 1:size(D,2)   % loop respect space grids
   D_data(:,i) = spline(t0,D(1:m1,i),t);
   V_data(:,i) = spline(t0,V(1:m1,i),t);
end

M = ceil(tfinal/(60*dt));
% %------------------------------------------------
% % initial conditions
% U(1,:) = D(1,:);
% U(2,:) = V(1,:);              % initial condition
% %------------------------------------------------

% M = length(t0);
% D_data = D;
% V_data = V;
% D_lbc = BD(:,1);
% D_rbc = BD(:,2);
% V_lbc = BV(:,1);
% V_rbc = BV(:,2);

% 
erho = zeros(1,M);
evel = zeros(1,M);
% x0 = x(1)-h;
% x1 = x(end)+h;

% % initial x
% xi = x(1)-gh_cell;
% xe = x(end)+gh_cell;

xi = x0(1);
xe = x0(end);

for n = 1:M
    
    %///////////////
    % Interpolation
    %///////////////

    % density from linear interpolation
    rho_l = D_lbc(n); rho_r = D_rbc(n);
    % insert into ghost cell [x(0), x(end+1)]
    
    rho = linear(x,xi,xe,rho_l,rho_r);    % x is the space vector
    
    num_rho(n) = rho(midd);
%     velocity predicted as linear interpolation
    vel_l = V_lbc(n); vel_r = V_rbc(n);
    vel = linear(x,xi,xe,vel_l,vel_r);  
    num_vel(n) = vel(midd);
%     vel0 = [vel_l,vel_r];
%     
%     rho = spline(x0,vel0,x);
    
    %\\\\\\\\\\\\\\\\\\
    % Calculate error 
    %//////////////////
%   [e1,e2] = comput_error_nagim_time(D_data(n,:),...
%             V_data(n,:),rho,vel);
     rho_diff = rho-D_data(n,:);
     vel_diff = vel-V_data(n,:);
%     
     e1 = norm(rho_diff,1);
     e2 = norm(vel_diff,1);
         
     erho(n) = e1;
     evel(n) = e2;
     
     %-------------------
     % do rescale
%      rho_diffs = (rho-D_data(n,:))./D_data(n,:);
%      vel_diffs = (vel-V_data(n,:))./V_data(n,:);
% %     
%      es1 = norm(rho_diffs,1);
%      es2 = norm(vel_diffs,1);
%          
%      erhos(n) = es1;
%      evels(n) = es2;
    
end

% \\\\\\\\\
% rescale
% /////////
% erho = erho/ref/rhom;
% evel = evel/ref/vm;

erho = erho/ref;
evel = evel/ref;
%\\\
% erhos = erhos/ref;
% evels = evels/ref;
%//////////////////////////////////////////////////////////////////////////
%//////////////////////////////////////////////////////////////////////////

% linear function with two points
function y = linear(x,x1,x2,y1,y2)

s = (y2-y1)/(x2-x1);
% b = (x2*y1-x1*y2)/(x2-x1);  % equivilent
b = y1-s.*x1;
y = s*x+b;
