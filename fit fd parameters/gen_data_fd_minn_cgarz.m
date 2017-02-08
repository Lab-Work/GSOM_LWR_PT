function P = gen_data_fd_minn_cgarz(rho0)

if nargin<1 rho0 = 38.5; end
rhom = 533;
%//////////////////////////////////////////////////////////////////////////
% This code is used to generate FD data from wrighted least squre fitting
% processes. This code can generate fd parameter either with 2 free
% parameters and 3 free parameters in the order [lambda, p, alpha].
% Shimao Fan, Temple University
% Jan 23 2012

% calculate fd for different days for one given sensor station
% k denotes sensor station
% sensor station 1~6 have 4 lanes
% sensor station 7~10 have 3 lanes
% 3 on ramp sensors and 3 off ramp sensor
% for 4 lanes piece of road
%//////////////////////////////////////////////////////////////////////////
% close all;
% clc;
% clear all;
%-----------------------------------------------------
% The maximum traffic density (maximum possible vel can fit into given 
%-----------------------------------------------------
% load data
load Minneapolis_data.mat   % every 30 seconds, V--volume, S---speed
% %---------------------------------------------------
% % deal with data and find out density and flux data D and Q
k = 2;
Q1 = squeeze(V(:,4*(k-1)+1,:));
V1 = squeeze(S(:,4*(k-1)+1,:));
for i = 2:4              % 4 lanes add up
  Q1 = Q1+squeeze(V(:,4*(k-1)+i,:));
  V1 = V1+squeeze(S(:,4*(k-1)+i,:)); 
end
% change units
V1 = 1.609344*V1;        % change into km/hour
Q1 = 120*Q1;             % change into #ofvehicle/hour
V1 = V1/4;               % average velocity of 4 lanes
flux = Q1; vel = V1;
% remove NAN number
density = flux./vel;
flux = flux(~isnan(density));
density = density(~isnan(density));
%------------------------------------------------------ 
% load data_minn_rho.mat   % density data
% load data_minn_flux.mat  % flux data
%==========================================
% The weighted least squre fitting process
%==========================================
%------------------------------------------
h = 1.0e-02;
% h = 1.0e-04;
% h = .01;
% beta = h/2:h:1-h/2;
% beta = .5;
m = 100;
% beta = linspace(h,1-h,m);
%------------------------------------------
% % test beta data
beta = [.001,.01,.03,.05,.1,.3,.5,0.99,0.999,0.9999];
% beta = [.000015,0.9,0.99,0.999,.9999];
%------------------------------------------
m = length(beta);

% beta(51)
%-----------------------------------------
%  3 parameters square root fd (original one)
% x0 = [20,.15,1200];   % (lambda, p, alpha);
%---------------------------------
% 3 parameter arctangent fd
% x0 = [10,100,10];
% P = zeros(m,3);  % initialize the matrix
% lb = [.1, 50 ,1];
% ub = [1000,rhom, 100]; 
%---------------------------------
% 2-p collapsed model, arctangent
% uf = 114; 
% uf = 94;
% rhof = 1026; 
%\\\\\\\\\\\\\\
uf = 95;
rhof = 555;
% rho0 = 38.5;
v0 = diff_fd_free(rho0,uf,rhof);%v0=collapse point velocity
f0 = fd_free(rho0,uf,rhof);%f0=collapse point density

P = zeros(m,2);  % initialize the matrix
x0 = [1 100];
% lb = [1.0e-03, rho0];
% ub = [1000,rhom]; 
%--------------------------------------
% to determine the free threshold rho_f
% x0 = 30;
% P = zeros(m,1);
% x  = 0;     % give a initial value
%--------------------------------------

for j = 1:m
   % 2-p collapsed model
   %----------------------------
   [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
   -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
   -flux),0)))^2,x0);
% pho_f=\pho_max^tilde, rho_0=rho_f£¬uf=v_max freeflow curve needs to be
% determined first
   %----------------------------
   % constrained minimization problem
   
%    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
%    -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
%    -flux),0)))^2,x0,[],[],[],[],lb,ub);
   %----------------------------
   P(j,:) = x;
   j
end

% fd_equi = P;
% rhoc = rho_critical(fd_equi(1),fd_equi(2),rhom)

% data_fd_ngsim_p3 have 3 parameters
% savefile = sprintf('data_fd_minn_%1.0f',rhom);
savefile = sprintf('data_fd_minn_cgarz_%1.0f',rho0);
% % savefile = sprintf('data_fd_minn_arctan_3p');
% 
save(savefile,'P');    % P is the parameters, n cross 3


%\\\\\\\\\\\\\\\\\\
%  The plots part
%//////////////////

figure(1)    % 

grid = 0:rhom;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',5), hold on
% plot(d_s1,f_s1,'r.'), hold on
% plot(d_s3,f_s3,'b*'), hold on

%---------------------------------------------------------
%  3 parameters
% plot(grid,traffic_flux_smooth(grid,P(1,1),P(1,2),P(1,3),rhom),'b-','linewidth',2.5),hold on
% plot(grid,traffic_flux_smooth(grid,P(2,1),P(2,2),P(2,3),rhom),'m-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth(grid,P(3,1),P(3,2),P(3,3),rhom),'r-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth(grid,P(4,1),P(4,2),P(4,3),rhom),'k-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth(grid,P(5,1),P(5,2),P(5,3),rhom),'r-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth(grid,P(6,1),P(6,2),P(6,3),rhom),'m-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth(grid,P(7,1),P(7,2),P(7,3),rhom),'b-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth(grid,P(8,1),P(8,2),P(8,3),rhom),'k-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth(grid,P(9,1),P(9,2),P(9,3),rhom),'b--','linewidth',2.5), hold on 
% plot(grid,traffic_flux_smooth(grid,P(10,1),P(10,2),P(10,3),rhom),'k--','linewidth',2.5), hold on  
% hold off
%--------------------------------------------------------
% % 2 parameters
% plot(grid,traffic_flux_smooth_2(grid,P(1,1),P(1,2),rhom),'b-','linewidth',2.5),hold on
% plot(grid,traffic_flux_smooth_2(grid,P(2,1),P(2,2),rhom),'r-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(3,1),P(3,2),rhom),'k-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(4,1),P(4,2),rhom),'r-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(5,1),P(5,2),rhom),'m-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(6,1),P(6,2),rhom),'b-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(7,1),P(7,2),rhom),'k-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(8,1),P(8,2),rhom),'b--','linewidth',2.5), hold on 
% plot(grid,traffic_flux_smooth_2(grid,P(9,1),P(9,2),rhom),'m-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(end,1),P(end,2),rhom),'k--','linewidth',2.5), hold on  
% hold off

% % 1 parameters
% plot(grid,traffic_flux_smooth_2(grid,P(1),p,rhom),'b-','linewidth',2.5),hold on
% plot(grid,traffic_flux_smooth_2(grid,P(2),p,rhom),'r-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(3),p,rhom),'k-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(4),p,rhom),'r-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(5),p,rhom),'m-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(6),p,rhom),'b-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(7),p,rhom),'k-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(8),p,rhom),'b--','linewidth',2.5), hold on 
% plot(grid,traffic_flux_smooth_2(grid,P(9),p,rhom),'m-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(end),p,rhom),'k--','linewidth',2.5), hold on  
% hold off
%--------------------------------------------------------
plot(grid,func_fd_seibold_2p(grid,P(1,1),P(1,2),rho0,v0,f0,uf,rhof,rhom),'r--','linewidth',4),hold on
plot(grid,func_fd_seibold_2p(grid,P(2,1),P(2,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p(grid,P(3,1),P(3,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p(grid,P(4,1),P(4,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p(grid,P(5,1),P(5,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p(grid,P(6,1),P(6,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p(grid,P(7,1),P(7,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p(grid,P(8,1),P(8,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p(grid,P(9,1),P(9,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on 
plot(grid,func_fd_seibold_2p(grid,P(end,1),P(end,2),rho0,v0,f0,uf,rhof,rhom),'b-.','linewidth',4),
hold off
%-------------------------------------------------------
% 3- arctangent function
% plot(grid,func_fd_seibold_3p(grid,P(1,1),P(1,2),P(1,3),rhom),'k--','linewidth',4),hold on
% plot(grid,func_fd_seibold_3p(grid,P(2,1),P(2,2),P(2,3),rhom),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_3p(grid,P(3,1),P(3,2),P(3,3),rhom),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_3p(grid,P(4,1),P(4,2),P(4,3),rhom),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_3p(grid,P(5,1),P(5,2),P(5,3),rhom),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_3p(grid,P(6,1),P(6,2),P(6,3),rhom),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_3p(grid,P(7,1),P(7,2),P(7,3),rhom),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_3p(grid,P(8,1),P(8,2),P(8,3),rhom),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_3p(grid,P(9,1),P(9,2),P(9,3),rhom),'k-','linewidth',4), hold on 
% plot(grid,func_fd_seibold_3p(grid,P(end,1),P(end,2),P(end,3),rhom),'k-.','linewidth',4),
% hold off
%--------------------------------------------------------
% plot(grid,func_det_rhof(grid,P(1),uf,rhof,rhom),'r--','linewidth',4),hold on
% plot(grid,func_det_rhof(grid,P(2),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(grid,func_det_rhof(grid,P(3),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(grid,func_det_rhof(grid,P(4),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(grid,func_det_rhof(grid,P(5),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(grid,func_det_rhof(grid,P(6),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(grid,func_det_rhof(grid,P(7),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(grid,func_det_rhof(grid,P(8),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(grid,func_det_rhof(grid,P(9),uf,rhof,rhom),'k-','linewidth',4), hold on 
% plot(grid,func_det_rhof(grid,P(end),uf,rhof,rhom),'b-.','linewidth',4),
% hold off
%--------------------------------------------------------
axis([0 rhom 0 11000])
% title('The G-AR model','fontsize',14)
title('Phase transition model','fontsize',14)%  title('Fitting to Smooth FDs W.R.T different \gamma','fontsize',14)
% legend('\omega = 0.01','\omega = 0.1','\omega = 0.5',...
%     '\omega = 0.9','\omega = 0.99')
%  legend('Data','\gamma = 0.001','\gamma = 0.01','\gamma = 0.1','\gamma = 0.5',...
%      '\gamma = 0.9','\gamma = 0.99','\gamma = 0.999')
%legend('Data','FD')
text(rhom-10,-200,'\rho_{max}','fontsize',14)
text(rhom/3,2000,'\Omega_{f}','fontsize',18)
text(rhom/2,8000,'\Omega_{c}','fontsize',18)

set(gca,'linewidth',2)
set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('density \rho','fontsize',14)
ylabel('flow rate Q (veh/h)','fontsize',14)
set(gca,'fontsize',14)
res = 600;
set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10 50 res*1.25 res])
set(gca,'position',[.1 .06 .88 .88])


%\\\\\\\\\\\\
%   figure(2)
% %////////////
% % plot(beta,P,'ro')
% plot(density,flux./density,'.','color',[0,0.5,0]), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(end),rhom)./grid,'k--','linewidth',4),hold on
% plot(grid,traffic_flux_trig(grid,vm,P(2),rhom)./grid,'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(3),rhom)./grid,'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(4),rhom)./grid,'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(5),rhom)./grid,'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(6),rhom)./grid,'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(7),rhom)./grid,'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(8),rhom)./grid,'k-','linewidth',2), hold on 
% plot(grid,traffic_flux_trig(grid,vm,P(9),rhom)./grid,'k-','linewidth',2), hold on
% plot(grid,traffic_flux_trig(grid,vm,P(1),rhom)./grid,'k-','linewidth',4),   
% hold off

%\\\\\\\\\\\\
% figure(3)
% %\\\\\\\\\\\\
% w = func_W(P,p,rhom,rhom);
% 
% plot(beta,P,'r--.')

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%//////////////////////////////////////////////////////////////////////////

%\\\\\\\\\\\\\\\\\\\\\\\\\\
%subfunctions
%//////////////////////////

% free flow branch
function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);

function y = diff_fd_free(rho,uf,rhof)

y = uf*(1-2*rho/rhof);