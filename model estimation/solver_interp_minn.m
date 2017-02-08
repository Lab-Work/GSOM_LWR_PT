%function [e_rho,e_vel,erho,evel] = solver_interp_minn(tfinal,position,D,Vel,h,w_max)
function [erho,evel] = solver_interp_minn(tfinal,position,D,Vel,h,w_max)

%==========================================================================
% The purpose of this code is to provide a comparision with traffic flow
% models. We will perform a test by applying linear interpolation, there is
% no traffic flow models applied here. Then we would like to see how the
% traffic models improve the prediction accuracy.
% Feb 06 2013
% Shimao Fan. Temple University
%==========================================================================

L = position(end);        % length of road in KM
% tfinal = 1.0;     % final time
N = ceil(L/h);    % number of grid
dx = L/N;         % final grid size
%------------------------------------------
k = .95*h/w_max;   % approximate time step
M = ceil(tfinal/k);         % number of time steps
dt = tfinal/M;              % here M is the number of time steps  
t = 0:dt:tfinal;                % time vector, uniform

%------------------------------------------
nofdata = tfinal*120;   % this is number of data points
inter = round(M/nofdata);   % how to pick the data
%------------------------------------------
% % location of interior sensor stations
% position = ceil(0.616/dx);     % position of sensor 2
% pos = .616;
pos = position(2:end-1);
% position = ceil(pos(2:end-1)/dx);     % position of sensor 2
%-------------------------------
% The boundary conditoins
% DensS1 = B(1,:);      % sensor 1, left bound
% DensS3 = B(3,:);      % sensor 3, right bound
% VelS1 = V(1,:);     % velocity of different time
% VelS3 = V(3,:);     % .. .. ..
% DensS2 = B(2,:);
% VelS2 = V(2,:);
for i = 1:length(position)   % number of sensor
    DensData(i,:) = intpol_bc(D(i,:),t);
    VelData(i,:) = intpol_bc(Vel(i,:),t);
end
% the boundary condition;
D_lbc = DensData(1,:); V_lbc = VelData(1,:);
D_rbc = DensData(end,:); V_rbc = VelData(end,:);

% N = size(D,2);
%--------------------------------
% find boundary condition by 
% the cubic interpolation/time
%--------------------------------
% density
% D_lbc = intpol_bc(DensS1,t);   % inflow
% D_rbc = intpol_bc(DensS3,t);   % outflow
% % velocity
% V_lbc = intpol_bc(VelS1,t);    % inflow
% V_rbc = intpol_bc(VelS3,t);    % outflow
% % sensor station 2
% DensData = intpol_bc(DensS2,t); 
% VelData = intpol_bc(VelS2,t); 



%\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Loop with respcet to time
%///////////////////////////
e_rho = zeros(1,M);
e_vel = zeros(1,M);
rho_num = zeros(1,M);
vel_num = zeros(1,M);
start = 5*inter;

for i = 1:M
% %     % density from linear interpolation
%     rho_l = D_lbc(i); rho_r = D_rbc(i);
%     rho = linear(pos,0,L,rho_l,rho_r);
%     
%     % velocity
%     vel_l = V_lbc(i); vel_r = V_rbc(i);
%     vel = linear(pos,0,L,vel_l,vel_r);
%----------------------------------------
    % density from linear interpolation
    rho_l = D_lbc(i); rho_r = D_rbc(i);
    rho = linear(pos,0,L,rho_l,rho_r);
%     rho_num(i) = rho;

%     velocity
    vel_l = V_lbc(i); vel_r = V_rbc(i);
    vel = linear(pos,0,L,vel_l,vel_r);
%     vel_num(i) = vel;

    data_rho = DensData(2:end-1,i)';
    data_vel = VelData(2:end-1,i)';
    diff_rho = data_rho-rho;
    diff_vel = data_vel-vel;
    
    
    e_rho(i) = norm(diff_rho,1);
    e_vel(i) = norm(diff_vel,1);
end

% error calculation
% start = 10*ceil(M/120);
% b
% start = 1;

% erho = mean(abs(DensData(start:M)-rho_num(start:M)))/rhom;
% evel = mean(abs(VelData(start:M)-vel_num(start:M)))/vm;

% erho = mean(abs(DensS2-rho_num))/rhom;
% evel = mean(abs(VelS2-vel_num))/vm;
% start = 1;
nn = length(position)-2;    % number of sensors validated; the interior ones

erho = mean(e_rho(start:M))/nn;
evel = mean(e_vel(start:M))/nn;


e_rho = e_rho(start:inter:end);
e_vel = e_vel(start:inter:end);
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Subfunctions 
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
function y = intpol_bc(data,t)   % of any order
    % data--- traffic data, traffic density
    % t --- time vector
    % assume picewise constant boundary conditions that came from real data
    tfinal = ceil(t(end)*60); % convert hour to minutes
    t = t*60;           % convert hour to minutes
    dt = 1/2;           % average avery 30 seconds, or half minute
    pl = 0:dt:tfinal;   % time interval we interested in for calculation
    n = length(pl);     % number of data we use
%     [param,res]= polyfit(pl,data,7);  % cl -- coefficient of lambda(w)
%     y = polyval(param,t);
    y = spline(pl,data(1:n),t);
%---------------------------------------
% linear function with two points
function y = linear(x,x1,x2,y1,y2)

s = (y2-y1)/(x2-x1);
b = (x2*y1-x1*y2)/(x2-x1);
y = s*x+b;
