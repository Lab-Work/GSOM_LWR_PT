function P = gen_data_fd_pt_ngsim(rhom)
if nargin<1 rhom = 800;j = 2; end

%==========================================================================
% Perform Weighted Least Squre fitting with NGSIM detector data. This code
% will generate a family of flow-density (FD) curves that with three free
% parameters [lambda,p,alpha].

% Shimao Fan, Temple University
% Nov 29 2011
% Modified on Jan 28 2013
%-------------------------------------------------------
% generate fd curves for Generic Phase Transition Model.
% Shimao, May 23 2013
%--------------------------------------------------------
% clc;
close all;

% Loading detector data part
load detector_data.dat 
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Deal with the data, find out Vol, and rho.
%////////////////////////////////////////////
B = detector_data;
% Pick the sensor station closest to the NGSIM study area, this is S#7
% pick data from sensor #7 and sensor #8
A_indx = find(B(:,1)>6);      % ith detector
A = B(A_indx,:);
[vol, rho,v] = flow_density(A);  % find volume and density data

% remove all NAN
indx = find(~isnan(rho));
D = rho(indx)';                % density data/sensor
Q = vol(indx)';                % flux data/sensor
V = v(indx)';                  % flux data/sensor
%----------------------------------------
% %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% Least square fitting process
%//////////////////////////////////

%----------------------------------------
% piecewise FD model
% fd_test(rho,sl,rhol,rhoc,rhom)
% sl -- slop of the free branch
% rhol -- intersection of linear and quadratic on left

% sl = 75;
% % eps = 0.01;
% eps = 1/rhom;
% % rhol = rhom/10;   % fix this value
% rhol = 58;
% x0 = rhom/5;         % critical density, initial guess
% h = 1.0e-03;
% beta = linspace(h,1-h,100);
%------------------------------------
% Seibold's 2 parameter model
% define free branch
% uf = 79; rhof = 1266; 
% uf = 73.6; rhof = 1266; 
% uf = 68;   rhof = 1057;
% uf = 71;  rhof = 1266;
uf = 78.7280;  rhof = 634.9755;
rho0 = 74.6000;
% rho0 = 67;  % I dont expect this parameter affect the result
% Rho0 = linspace(50,100,10);
% rho0 = Rho0(j);
% rho0 = 56;   % rho threshold
% rho0 = 63;
v0 = diff_fd_free(rho0,uf,rhof);
f0 = fd_free(rho0,uf,rhof);
% v0 = 79; rho0 = 58; f0 = v0*rho0;
% q_max = f0+v0*(rhom-rho0);
% D = [D,rhom];
% Q = [Q,q_max];
x0 = [1,100];
%\\\\\\\\\\\\\\\
% select data
%///////////////
% indx1 = find(D>rho0);
% D1 = D(indx1);
% Q1 = Q(indx1);

% beta = [.0001,.001,.01,.03,.005,.1,.3,.5,0.99,.999]; % sample test case
h = 1.0e-03;
m = 100;
beta = linspace(h,1-h,m);

% beta = [0.01,.01,.02,.03,.1,.3,.5,0.6,.9,0.99]; % sample test case

% lower bound and upper bound
lb = [1.0e-03, rho0];
% lb = [1.0, rho0];
ub = [1000,rhom]; 
m = length(beta);
P = zeros(m,2);

%\\\\\\\\\\\\\\\\\\\
% 3 parameter model
%///////////////////
% x0 = [10,100,10];
% P = zeros(m,3);

%\\\\\\\\\\\\\\\\\\\\\\\\\
% determine rho thresold
%\\\\\\\\\\\\\\\\\\\\\\\\
% x0 = 40;
% P = zeros(m,1);
% x = 0;
% 
% beta0 = 10.^(linspace(-3,-4,m));
% beta = 1-beta0;
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% m = 0;
for j = 1:m
    
    %=====================
    % piecewise quadratic
    %=====================
%    [x,fval] = fminsearch(@(x)beta(j)*norm(max((fd_test(D,sl,rhol,x,rhom)...
%    -Q),0*Q))^2+(1-beta(j))*norm(max(-(fd_test(D,sl,rhol,x,rhom)...
%    -Q),0*Q))^2,x0);

%    [x,fval] = fminsearch(@(x)beta(j)*norm(max((fd_piecewise_quadratic(D,sl,rhol,x,rhom,eps)...
%    -Q),0*Q))^2+(1-beta(j))*norm(max(-(fd_piecewise_quadratic(D,sl,rhol,x,rhom,eps)...
%    -Q),0*Q))^2,x0);

%%%%%%%%%%%%%%%%%%%%%
% seibold's 2-p model\\\fminunc
%%%%%%%%%%%%%%%%%%%%%
% non-constrained minimization problem
%    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_fd_seibold(D1,x(1),x(2),v0,rho0,f0,rhom)...
%    -Q1),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold(D1,x(1),x(2),v0,rho0,f0,rhom)...
%    -Q1),0)))^2,x0);
%--------------------------------------
%    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_fd_seibold_2p(D,x(1),x(2),rho0,uf,rhof,rhom)...
%    -Q),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(D,x(1),x(2),rho0,uf,rhof,rhom)...
%    -Q),0)))^2,x0);
% % constrained minimization problem
   [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(D,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
   -Q),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(D,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
   -Q),0)))^2,x0,[],[],[],[],lb,ub);

%%%%%%%%%%%%%%%%%%%%%
% seibold's 3-p model\\\fminunc
%%%%%%%%%%%%%%%%%%%%%
%    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_fd_seibold_3p(D,x(1),x(2),x(3),rhom)...
%    -Q),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_3p(D,x(1),x(2),x(3),rhom)...
%    -Q),0)))^2,x0);
%\\\\\\\\\\\\\\\\\\\\\\\\\\\
% determine rho threshold
%///////////////////////////
%    indx1 = find(D>x);
%    D1 = D(indx1);
%    Q1 = Q(indx1);
%    %-------------------
%    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_det_rhof(D1,x,uf,rhof,rhom)...
%    -Q1),0)))^2+(1-beta(j))*(norm(max(-(func_det_rhof(D1,x,uf,rhof,rhom)...
%    -Q1),0)))^2,x0);

   P(j,:) = x;
 % P(j,2) = x;   
   j
%    fval
end

% P(:,1) = 20;
% x
%\\\\\\\\\\\\\\\\\\
% Save the results 
%//////////////////
% savefile = sprintf('data_fd_ngsim_rhom_%03d',rhom);

% savefile = sprintf('data_fd_ngsim_trig_%03d',rhom);

% savefile = sprintf('data_fd_ngsim_pt_piecewise');
% savefile = sprintf('data_fd_ngsim_pt_seibold_2p');

% savefile = sprintf('data_fd_ngsim_arctan_2p_test');

savefile = sprintf('data_fd_cgarz_ngsim_%1.0f',rhom);

% % % 
save(savefile,'P');   

% load data_fd_ngsim_arctan_3p.mat

figure(1)
%\\\\\\\\\\\\\\\\\\\
%  The plots part
%///////////////////
rho = 0:rhom;
rho1 = 0:rho0;
rho2 = rho0:.1:rhom;    % plot of the sequence of flow-density curves
bound = linear_bound(rho2,v0,rho0,f0);
plot(rho2,bound,'b-','linewidth',6), hold on
plot(D,Q,'.','color',[0 .4 0],'markersize',8), hold on
plot(rho1,fd_free(rho1,uf,rhof),'r-','linewidth',4), hold on
% plot(rho1,fd_free(rho1-10,uf,rhof),'b--','linewidth',4), hold on
% ---------------------------------------------------------
%  piecewise C1 model
% plot(grid,fd_test(grid,sl,rhol,P(1),rhom),'b-','linewidth',4),hold on
% plot(grid,fd_test(grid,sl,rhol,P(2),rhom),'m-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(3),rhom),'r-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(4),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(5),rhom),'r-.','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(6),rhom),'m-.','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(7),rhom),'b-.','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(8),rhom),'k-.','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(9),rhom),'b--','linewidth',4), hold on 
% plot(grid,fd_test(grid,sl,rhol,P(10),rhom),'k--','linewidth',4), hold on  
% hold off

% plot(grid,fd_test(grid,sl,rhol,P(1),rhom),'k-','linewidth',4),hold on
% plot(grid,fd_test(grid,sl,rhol,P(2),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(3),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(4),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(5),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(6),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(7),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(8),rhom),'k-','linewidth',4), hold on
% plot(grid,fd_test(grid,sl,rhol,P(9),rhom),'k-','linewidth',4), hold on 
% plot(grid,fd_test(grid,sl,rhol,P(end),rhom),'k-','linewidth',4), hold on  
% hold off
%\\\\\\\\\\\\\\\\\\\\\\

% plot(rho2,func_fd_seibold(rho2,P(1,1),P(1,2),v0,rho0,f0,rhom),'r--','linewidth',4),hold on
% plot(rho2,func_fd_seibold(rho2,P(2,1),P(2,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on
% plot(rho2,func_fd_seibold(rho2,P(3,1),P(3,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on
% plot(rho2,func_fd_seibold(rho2,P(4,1),P(4,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on
% plot(rho2,func_fd_seibold(rho2,P(5,1),P(5,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on
% plot(rho2,func_fd_seibold(rho2,P(6,1),P(6,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on
% plot(rho2,func_fd_seibold(rho2,P(7,1),P(7,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on
% plot(rho2,func_fd_seibold(rho2,P(8,1),P(8,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on
% plot(rho2,func_fd_seibold(rho2,P(9,1),P(9,2),v0,rho0,f0,rhom),'k-','linewidth',4), hold on 
% plot(rho2,func_fd_seibold(rho2,P(end,1),P(end,2),v0,rho0,f0,rhom),'b-.','linewidth',4),
% hold off
%\\\\\\\\\\\\\\\\\\\\\
plot(rho,func_fd_seibold_2p(rho,P(1,1),P(1,2),rho0,v0,f0,uf,rhof,rhom),'r--','linewidth',4),hold on
plot(rho,func_fd_seibold_2p(rho,P(2,1),P(2,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(rho,func_fd_seibold_2p(rho,P(3,1),P(3,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(rho,func_fd_seibold_2p(rho,P(4,1),P(4,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(rho,func_fd_seibold_2p(rho,P(5,1),P(5,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(rho,func_fd_seibold_2p(rho,P(6,1),P(6,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(rho,func_fd_seibold_2p(rho,P(7,1),P(7,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(rho,func_fd_seibold_2p(rho,P(8,1),P(8,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on
plot(rho,func_fd_seibold_2p(rho,P(9,1),P(9,2),rho0,v0,f0,uf,rhof,rhom),'k-','linewidth',4), hold on 
plot(rho,func_fd_seibold_2p(rho,P(end,1),P(end,2),rho0,v0,f0,uf,rhof,rhom),'b-.','linewidth',4),
hold off
%\\\\\\\\\\\\\\\\\\\\\\\
% plot(rho,func_fd_seibold_3p(rho,P(1,1),P(1,2),P(1,3),rhom),'k--','linewidth',4),hold on
% plot(rho,func_fd_seibold_3p(rho,P(2,1),P(2,2),P(2,3),rhom),'k-','linewidth',4), hold on
% plot(rho,func_fd_seibold_3p(rho,P(3,1),P(3,2),P(3,3),rhom),'k-','linewidth',4), hold on
% plot(rho,func_fd_seibold_3p(rho,P(4,1),P(4,2),P(4,3),rhom),'k-','linewidth',4), hold on
% plot(rho,func_fd_seibold_3p(rho,P(5,1),P(5,2),P(5,3),rhom),'k-','linewidth',4), hold on
% plot(rho,func_fd_seibold_3p(rho,P(6,1),P(6,2),P(6,3),rhom),'k-','linewidth',4), hold on
% plot(rho,func_fd_seibold_3p(rho,P(7,1),P(7,2),P(7,3),rhom),'k-','linewidth',4), hold on
% plot(rho,func_fd_seibold_3p(rho,P(8,1),P(8,2),P(8,3),rhom),'k-','linewidth',4), hold on
% plot(rho,func_fd_seibold_3p(rho,P(9,1),P(9,2),P(9,3),rhom),'k-','linewidth',4), hold on 
% plot(rho,func_fd_seibold_3p(rho,P(end,1),P(end,2),P(end,3),rhom),'k-.','linewidth',4),
% hold off

% plot(rho,func_det_rhof(rho,P(1),uf,rhof,rhom),'r--','linewidth',4),hold on
% plot(rho,func_det_rhof(rho,P(2),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(rho,func_det_rhof(rho,P(3),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(rho,func_det_rhof(rho,P(4),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(rho,func_det_rhof(rho,P(5),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(rho,func_det_rhof(rho,P(6),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(rho,func_det_rhof(rho,P(7),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(rho,func_det_rhof(rho,P(8),uf,rhof,rhom),'k-','linewidth',4), hold on
% plot(rho,func_det_rhof(rho,P(9),uf,rhof,rhom),'k-','linewidth',4), hold on 
% plot(rho,func_det_rhof(rho,P(end),uf,rhof,rhom),'b-.','linewidth',4),
% hold off

%-----------------------------------------------
axis([0 rhom 0 14000])
% title('Traffic data collected from I-80 in Emeryville(D3), CA','fontsize',14)%  title('Fitting to Smooth FDs W.R.T different \gamma','fontsize',14)
% title('Phase transition model','fontsize',14)
% legend('\omega = 0.01','\omega = 0.1','\omega = 0.5',...
%     '\omega = 0.9','\omega = 0.99')
%  legend('Data','\gamma = 0.001','\gamma = 0.01','\gamma = 0.1','\gamma = 0.5',...
%      '\gamma = 0.9','\gamma = 0.99','\gamma = 0.999')
%legend('Data','FD')
text(rhom-20,-300,'\rho_{max}','fontsize',14)
% text(rhom/3,3000,'\Omega_{f}','fontsize',18)
% text(rhom/2,9000,'\Omega_{c}','fontsize',18)
set(gca,'linewidth',2)
% set(gca,'xtick',[])
% set(gca,'ytick',[])
xlabel('density \rho','fontsize',14)
ylabel('flow rate Q (veh/h)','fontsize',14)
set(gca,'fontsize',14)
res = 600;
filename_save = sprintf('fig_fd_rhom_%03d',rhom);

set(gcf,'paperpositionmode','auto')
set(gcf,'position',[10 50 res*1.15 res])
set(gca,'position',[.08 .06 .90 .88])
print(gcf,'-dpng',filename_save,'-r290')


%-----------------------------
% figure(2)
% plot(beta,P(:,1),'r--.','linewidth',3)
% 
% 
% figure(3)
% plot(beta,P(:,2),'r--.','linewidth',3)
% figure(2)
% q = func_q(P(:,1),P(:,2),rho0,uf,rhof,rhom);
% plot(q,P(:,1),'k--.','linewidth',4), hold on
% plot(q,P(:,2),'r-','linewidth',4), hold on
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% w = func_w_arctan(P(:,1),P(:,2),P(:,3),rhom);
% 
% plot(beta(1:end),w,'r-','linewidth',3), hold on
% plot(beta(1:end),P(:,1),'k--','linewidth',4), hold on
% plot(beta(1:end),P(:,2),'k-','linewidth',4), hold on
% plot(beta(1:end),P(:,3),'k--','linewidth',4), hold off
% 
% figure(3)
% plot(w,P(:,1),'k--.','linewidth',4), hold on
% plot(w,P(:,2),'r-','linewidth',4), hold on
% plot(w,P(:,3),'k--.','linewidth',4), hold off

%\\\\\\\\\\\\\\\\\\
%  sub_functions
%//////////////////

%==========================================================================
%==========================================================================
%==========================================================================
function [Vol, rho, V] = flow_density(A) % find average velocity and volum

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
    
% free flow branch
function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);

function y = diff_fd_free(rho,uf,rhof)

y = uf*(1-2*rho/rhof);

function y = linear_bound(rho,v0,rho0,f0)

y = f0+v0*(rho-rho0);




