
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
%Ye Sun
%RELAX RHOM
%//////////////////////////////////////////////////////////////////////////
% close all;
% clc;
% clear all;
%---------------------------------------------------
% The maximum traffic density (maximum possible vel can fit into given 
%---------------------------------------------------
% load data
close all;
% Loading detector data part
load Minneapolis_data.mat   % every 30 seconds, V--volume, S---speed
% %---------------------------------------------------
% % deal with data and find out density and flux data D and Q
k = 2;
fd=1;
ld=79;
len=6;
Q1 = squeeze(V([fd:len:ld],4*(k-1)+1,:));
V1 = squeeze(S([fd:len:ld],4*(k-1)+1,:));
for i = 2:4              % 4 lanes add up
  Q1 = Q1+squeeze(V([fd:len:ld],4*(k-1)+i,:));
  V1 = V1+squeeze(S([fd:len:ld],4*(k-1)+i,:)); 
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
index0 = find(density(:,1)<=40 | flux(:,1)>=3000);

flux = flux(index0);
density = density(index0);
% index1= find(density(:,1)<=60 | flux(:,1)>=4000 | density(:,1)>=260)
% 
% flux = flux(index1);
% density = density(index1);

% f1=min(flux);
% flux=flux-f1;
% r1=min(density);
% density=density-r1;
size(density)

% load data_minn_rho.mat   % density data
% load data_minn_flux.mat  % flux data
%==========================================
% The weighted least squre fitting process
%==========================================

h = 1.0e-03;
%h = 1.0e-03;
beta = linspace(h,1-h,100);
m=length(beta);
eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);
%------------------------------------------


count=0;

load CGARZ_eq_6para_updated.mat;
rho0 = para_eq(3,1);
% uf = 114; 
uf = para_eq(4,1);
rhof = para_eq(5,1);
rho_m=para_eq(6,1);

load CGARZ_100beta_2para_updated.mat;




grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.5,0.5,0.5],'markersize',6), hold on



for i=1:m
    if mod(i-1,10)==0 || i==m       
%        plot(grid,func_fd_seibold_2p_updated(grid,P2(i,1),P2(i,2),rho0,uf,rhof,rho_m),'k-','linewidth',0.05),hold on
    end
end

plot(grid,func_fd_seibold_2p_updated(grid,P2(21,1),P2(21,2),rho0,uf,rhof,rho_m),'k-','linewidth',2.5),hold on
plot(grid,func_fd_seibold_2p_updated(grid,P2(80,1),P2(80,2),rho0,uf,rhof,rho_m),'k-','linewidth',2.5),hold on
plot(grid,func_fd_seibold_2p_updated(grid,para_eq(1,1),para_eq(2,1),rho0,uf,rhof,rho_m),'r-','linewidth',2.5),hold on
plot(rho0,5,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
hold off

