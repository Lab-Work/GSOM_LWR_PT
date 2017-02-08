function P = gen_data_fd_ngsim_free
%==========================================================================
% Perform Weighted Least Squre fitting with NGSIM detector data. This code
% will generate a family of flow-density (FD) curves that with three free
% parameters [lambda,p,alpha].

% Shimao Fan, Temple University
% Nov 29 2011
% Modified on Jan 28 2013
%==========================================================================
close all;
rhom = 800;
% Loading detector data part
load detector_data.dat 
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Deal with the data, find out Vol, and rho.
%////////////////////////////////////////////
B = detector_data;
% Pick the sensor station closest to the NGSIM study area, this is S#7
% i = 7;                         % study sensor station #7
% pick data from sensor #7 and sensor #8
A_indx = find(B(:,1)>6);      % ith detector
% % A_indx = find(B(:,1)==8);      % ith detector
% 
A = B(A_indx,:);
[vol, rho] = flow_density(A);  % find volume and density data
% [vol, rho] = flow_density(B);  % find volume and density data

% remove all NAN
indx = find(~isnan(rho));
D = rho(indx)';                % density data/sensor
Q = vol(indx)';                % flux data/sensor
q_max = max(max(Q));           % ???
% %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Least square fitting process
%//////////////////////////////////

%----------------------------------------
% initial guess [um = 80, rhom = 700].

% x0 = [70,800];   % (lambda, p, alpha);

% daganzo's model
x0 = [80 rhom/4]; % ??? how is initial condition chosen?
%---------------------------------
% for test case (selected samples)
beta = [.001,.01,.1,.4,.5,.8,.9,.99,0.999,0.9999];
%---------------------------------
% another test case (Uniform)
h = 1.0e-03;
% beta = [h,.5,1-h];
beta = linspace(h,1-h,10);
%----------------------------------
% equilibrim cruve, Dangazon's model
beta = 0.5;
% vm = 66.2708;  rhoc = 137.0164;
% x0 = [rhoc vm];
%--------------------------------
x0 = 0;
vm = 66.5019;  rhoc = 135.3308;
% beta = linspace(-1,1,10);
beta = linspace(h,1-h,10);
% beta = 0.5;
%----------------------------------
m = length(beta);
P = zeros(m,1);
% P = zeros(m,2);
% format long
%--------------------------------------------
% determine free branch of the cgarz model
% rhoc = 200;
% beta = 0.5;
% x0 = [80,1260];
% P = zeros(m,2);
% index = find(D<=rhoc);
% D_free = D(index);
% Q_free = Q(index);
%-----------------------------------

for j = 1:m
  %%%%%%%%%%%%%%  
  % 3 parameter
%   %%%%%%%%%%%%%%  
%   [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_invers_q(Q,D,x(1),x(2))...
%    -D),0*D))^2+(1-beta(j))*norm(max(-(func_invers_q(Q,D,x(1),x(2))...
%    -D),0*D))^2,x0);
%----------------------
%%%%%%%%%%%%%%%%%%
   % 2 parameters
%%%%%%%%%%%%%%%%%%   
%   [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_smooth_2(D,x(1),x(2),rhom)...
%    -Q),0*Q))^2+(1-beta(j))*norm(max(-(traffic_flux_smooth_2(D,x(1),x(2),rhom)...
%    -Q),0*Q))^2,x0);NOT THIS

% [x,fval] = fminsearch(@(x)lambda(j)*norm(max((traffic_flux_smooth(D2,x(1),x(2),rhom,alpha)...
%    -Q2),0*Q2))^2+(1-lambda(j))*norm(max(-(traffic_flux_smooth(D2,x(1),x(2),rhom,alpha)...
%    -Q2),0*Q2))^2,x0);
%---------------------------
% 1 parameter quadratic form
% [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(D,x(1),rhom)...
%    -Q),0*Q))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(D,x(1),rhom)...
%    -Q),0*Q))^2,x0);
%-----------------------
% piecewise quadratic form (phase transition model), 1 parameters
% [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig(D,vm,x(1),rhom)...
%    -Q),0*Q))^2+(1-beta(j))*norm(max(-(traffic_flux_trig(D,vm,x(1),rhom)...
%    -Q),0*Q))^2,x0);
%-----------------------
% Daganzo--Newell model
% [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig_1(D,x(1),x(2),rhom)...
%    -Q),0*Q))^2+(1-beta(j))*norm(max(-(traffic_flux_trig_1(D,x(1),x(2),rhom)...
%    -Q),0*Q))^2,x0);

%-----------------------
% greenshield phase transition model
[x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(D,x,rhoc,vm,rhom)...
   -Q),0*Q))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(D,x,rhoc,vm,rhom)...
   -Q),0*Q))^2,x0);
%!!!Greenshield PT!!!

% [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(D,0,x(1),x(2),rhom)...
%    -Q),0*Q))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(D,0,x(1),x(2),rhom)...
%    -Q),0*Q))^2,x0);
%------------------------
% best fitting curve respect to the free data
% a parameter quadratic form
% [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(D_free,x(1),x(2))...
%    -Q_free),0*Q_free))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(D_free,x(1),x(2))...
%    -Q_free),0*Q_free))^2,x0);
% x
%    P(j,:) = x
   
   P(j) = x;

%    j
end
P
% x
%\\\\\\\\\\\\\\\\\\
% Save the results 
%//////////////////
% savefile = sprintf('data_fd_ngsim_test_%03d',rhom);

% savefile = sprintf('data_fd_ngsim_trig_%03d',rhom);

% savefile = sprintf('data_fd_ngsim_2_%03d',rhom);
% savefile = sprintf('data_fd_ngsim_free');
% % 
% save(savefile,'P');   


%\\\\\\\\\\\\\\\\\\\
%  The plots part
% %///////////////////
% grid = 0:100:q_max;    % plot of the sequence of flow-density curves
grid = 0:rhom;
% grid = 
plot(D,Q,'.','color',[0 .5 0],'markersize',8), hold on
% plot(D_s,Q_s,'r.'), hold on
% ---------------------------------------------------------
% %  3 parameters smooth model
% plot(grid,func_flux(grid,P(1,1),P(1,2)),'k-','linewidth',4),hold on
% plot(grid,func_flux(grid,P(2,1),P(2,2)),'m-','linewidth',2.5), hold on
% plot(grid,func_flux(grid,P(3,1),P(3,2)),'r-','linewidth',2.5), hold on
% plot(grid,func_flux(grid,P(4,1),P(4,2)),'b-','linewidth',2.5), hold on
% plot(grid,func_flux(grid,P(5,1),P(5,2)),'k--','linewidth',2.5), hold on
% plot(grid,func_flux(grid,P(6,1),P(6,2)),'m--','linewidth',2.5), hold on
% plot(grid,func_flux(grid,P(7,1),P(7,2)),'r--','linewidth',2.5), hold on
% plot(grid,func_flux(grid,P(8,1),P(8,2)),'b--','linewidth',2.5), hold on
% plot(grid,func_flux(grid,P(9,1),P(9,2)),'k-.','linewidth',2.5), hold on 
% plot(grid,func_flux(grid,P(end,1),P(end,2)),'r-.','linewidth',4), hold on  
% hold off
%-----------------------------------------------
% plot(grid,traffic_flux_smooth_2(grid,P(1,1),P(1,2),rhom),'b-','linewidth',2.5),hold on
% plot(grid,traffic_flux_smooth_2(grid,P(2,1),P(2,2),rhom),'m-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(3,1),P(3,2),rhom),'r-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(4,1),P(4,2),rhom),'k-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(5,1),P(5,2),rhom),'r-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(6,1),P(6,2),rhom),'m-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(7,1),P(7,2),rhom),'b-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(8,1),P(8,2),rhom),'k-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(9,1),P(9,2),rhom),'b--','linewidth',2.5), hold on 
% plot(grid,traffic_flux_smooth_2(grid,P(10,1),P(10,2),rhom),'k--','linewidth',2.5), hold on  
% hold off
%------------------------------------------------
% 1 parameters trig flux functions
% plot(grid,traffic_flux_trig_1(grid,P(end,1),P(end,2),rhom),'b--','linewidth',4),hold on
% plot(grid,traffic_flux_trig_1(grid,P(2,1),P(2,2),rhom),'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig_1(grid,P(3,1),P(3,2),rhom),'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig_1(grid,P(4,1),P(4,2),rhom),'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig_1(grid,P(5,1),P(5,2),rhom),'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig_1(grid,P(6,1),P(6,2),rhom),'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig_1(grid,P(7,1),P(7,2),rhom),'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig_1(grid,P(8,1),P(8,2),rhom),'k-','linewidth',2), hold on 
% plot(grid,traffic_flux_trig_1(grid,P(9,1),P(9,2),rhom),'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig_1(grid,P(1,1),P(1,2),rhom),'k--','linewidth',4),   
% hold off
%-------------------------------------------------
plot(grid,func_flux_creen_pt(grid,P(1),rhoc,vm,rhom),'b--','linewidth',4),hold on
plot(grid,func_flux_creen_pt(grid,P(2),rhoc,vm,rhom),'k-','linewidth',2), hold on
plot(grid,func_flux_creen_pt(grid,P(3),rhoc,vm,rhom),'k-','linewidth',2), hold on
plot(grid,func_flux_creen_pt(grid,P(4),rhoc,vm,rhom),'k-','linewidth',2), hold on
plot(grid,func_flux_creen_pt(grid,P(5),rhoc,vm,rhom),'k-','linewidth',2), hold on
plot(grid,func_flux_creen_pt(grid,P(6),rhoc,vm,rhom),'k-','linewidth',2), hold on
plot(grid,func_flux_creen_pt(grid,P(7),rhoc,vm,rhom),'k-','linewidth',2), hold on
plot(grid,func_flux_creen_pt(grid,P(8),rhoc,vm,rhom),'k-','linewidth',2), hold on 
plot(grid,func_flux_creen_pt(grid,P(9),rhoc,vm,rhom),'k-','linewidth',2), hold on
plot(grid,func_flux_creen_pt(grid,P(10),rhoc,vm,rhom),'k--','linewidth',4),   
hold off

% plot(grid,func_flux_creen_pt(grid,0,P(1,1),P(1,2),rhom),'b--','linewidth',4),hold on

%-----------------------------------------------
axis([0 800 0 q_max])
% title('Traffic data collected from I-80 in Emeryville(D3), CA','fontsize',14)%  title('Fitting to Smooth FDs W.R.T different \gamma','fontsize',14)
title('Phase transition model','fontsize',14)
% legend('\omega = 0.01','\omega = 0.1','\omega = 0.5',...
%     '\omega = 0.9','\omega = 0.99')
%  legend('Data','\gamma = 0.001','\gamma = 0.01','\gamma = 0.1','\gamma = 0.5',...
%      '\gamma = 0.9','\gamma = 0.99','\gamma = 0.999')
%legend('Data','FD')
% text(rhom-20,-300,'\rho_{max}','fontsize',14)
% text(rhom/3,3000,'\Omega_{f}','fontsize',18)
% text(rhom/2,9000,'\Omega_{c}','fontsize',18)
set(gca,'linewidth',2)
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('density \rho','fontsize',14)
ylabel('flow rate Q (veh/h)','fontsize',14)
set(gca,'fontsize',14)
res = 600;
filename_save = sprintf('fig_phase_trans_curves');

set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10 50 res*1.25 res])
set(gca,'position',[.1 .06 .88 .88])
print(gcf,'-dpng',filename_save,'-r290')

% % 
% w = func_W(P(:,1),P(:,2),P(:,3),rhom);
% max(w)
% 
% figure(2)
% plot(beta,w,'k-.')

%\\\\\\\\\\\\\\\\\\
%  sub_functions
%//////////////////

%==========================================================================
%==========================================================================
%==========================================================================
function [Vol, rho] = flow_density(A) % find average velocity and volum

%=========================================
%     Volume over all lanes/30 seconds
%=========================================
Vol = A(:,34)+A(:,35)+A(:,36)+A(:,37)+A(:,38);   % # of cars/30 seconds
% Vol = A(:,29)+A(:,30)+A(:,31)+A(:,32)+A(:,33);   % # of cars/30 seconds

% change unit to # of cars/hour
Vol = 120*Vol;
%======================================
%     Average velocity/30 seconds
%======================================
Vel = A(:,8:5:28);         % store all the velocity in differents lanes
Nv = zeros(size(A(:,1)));  % nonzero velocity
v = zeros(size(A(:,1)));   % velocity vector
for i = 1:length(A(:,1))   % loop over all the rows
    indx = find(Vel(i,:)); % find index of nonzero velociy
    Nv(i) = length(indx);  % # of nonzero velocity
    v(i) = sum(Vel(i,:));
end

V = v./Nv;     % feet/second
% change unit  % km/hour
V = 1.09728*V;
rho = Vol./V;
    

%-------------------------------
function y = func_flux(rho,um,rhom)

y = um*rho.*(1-rho/rhom);




