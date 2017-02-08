function P1 = gen_data_fd_ngsim_sensors_relax_ic_cgarz_newFFpara1


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
clear;
load CGARZ_varyingIC-realistic_rhom0-ngsim-at.mat
beta = [.001,.01,.03,.05,.1,.3,.5,0.9,0.99];


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
density = rho(indx)';                % density data/sensor
flux = vol(indx)';                % flux data/sensor
q_max = max(max(flux));           % 
 

% load data_minn_rho.mat   % density data
% load data_minn_flux.mat  % flux data
%==========================================
% The weighted least squre fitting process
%==========================================

h = 1.0e-02;
% h = 1.0e-03;
% h = .01;
% beta = h/2:h:1-h/2;
% beta = .5;
%m = 100;
% beta = linspace(h,1-h,2);
% beta = [h,.5,1-h];
% % test beta data
beta = [.001,.01,.03,.05,.1,.3,.5,0.9,0.99];
% beta = [.000015,0.9,0.99,0.999,.9999];
%------------------------------------------
%YE GARZ-------------------------------------------------------
%  3 parameters square root fd (original one)
% lambda0=20;
% p0=0.5;
% alpha0=1000;
% rhom0=10;
% x0 = [lambda0,p0,alpha0,rhom0];   % (lambda, p, alpha,rhom);
% x1 = [lambda0,p0,alpha0];   % (lambda, p, alpha);
% beta = linspace(h,1-h,100);
%YE GARZ-------------------------------------------------------

%---------------------------------
% 3 parameter arctangent fd
% x0 = [10,100,10];
% P = zeros(m,3);  % initialize the matrix
% lb = [.1, 50 ,1];
% ub = [1000,rhom, 100]; 
%---------------------------------
% 2-p collapsed model, arctangent
% beta = linspace(h,1-h,100);
% 
% % uf = 114; rhof = 1026; 
% % uf = 104.4861; rhof = 566.9842;
% % rho0 = 48.5;
% % rho0 = 33;
%YE CGARZ-------------------------------------------------------

rho0 = 30.000000005234; uf = 76.7063878200036;  rhof = 622.7951046329340;
v0 = diff_fd_free(rho0,uf,rhof);
f0 = fd_free(rho0,uf,rhof);
rho_m=777.2992521252818;

lb = [1.0e-03,rho0];
ub = [1000,2000]; 
lb1 = [1.0e-03, rho0];
ub1 = [1000,2000]; 

%YE CGARZ-------------------------------------------------------
%--------------------------------------
% newell-danganzo model
% x0 = [100 rhom/4];
% beta = 0.5;
%--------------------------------------
% free flow branch of cgarz model
% vm = 93.6380;   rhoc = 81.6558;
% x0 = [vm rhom];
% beta = 0.5;
% index = find(density<=rhoc);
% D_free = density(index);
% Q_free = flux(index);
%---------------------------------------
% determine standard curve in gpt
% x0 = [rhom/4 95];
% beta = 0.5;
%-------
% a family of perturbed curves
% rhoc = 80.6972;   vm = 93.7165;
% beta = [h 1-h];
% x0=0;
%-----------------------------------
% a quadratic flux curve
%x0 = 90;
%beta = 0.5;
%---------------------------------
%m = length(beta);
%P = zeros(m,1);  % initialize the matrix


%---------------------------------

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%YE GARZ-------------------------------------------------------
% options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
% %options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
%  lbeq = [0,0,0,0];
% %options = optimset('MaxFunEvals', 1*10^10);
%  [x,fval] = fmincon(@(x)beta(7)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
%     -flux),0*flux))^2+(1-beta(7))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
%     -flux),0*flux))^2,x0,[],[],[],[],lbeq,[],[],options);
% %  [x,fval] = fminsearch(@(x)beta(7)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
% %     -flux),0*flux))^2+(1-beta(7))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
% %     -flux),0*flux))^2,x0,options);
% 
%   P(7,1) = x(1);
%   P(7,2) = x(2);
%   P(7,3) = x(3);
%   
%   rho_m = x(4);
%   P(7,1)
%   P(7,2)
%   P(7,3)
%   rho_m
%YE GARZ-------------------------------------------------------

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\


% P = zeros(m,3);  % initialize the matrix
%YE CGARZ--------------------------------count=0;

%YE CGARZ-------------------------------------------------------
initials = [1,10,50,100,500,1000];
initialm = [50,100,500,1000,2000];
% initials = [100];
% initialm = [100];
% initialr = [500,600];
N1 = length(initials);
N2=length(initialm);
N=N1*N2;
m = length(beta);
P2 = zeros(N,3*9+2);  % initialize the matrix-----------------------

%YE GARZ-----------------------------------------------------------------
% initiall = [1,50,100,500,1000,4000];
% initialp = [0.1,0.5,0.9];
% initiala = [1,50,100,500,1000,4000];
% initialr = [400,600,800];
% % initiall = [50];
% % initialp = [0.5];
% % initiala = [1000];
% % initialr = [600,400];
% lbeq = [0,0,0,0];
% lb1 = [0,0,0];
% N1=length(initiall);
% N2=length(initialp);
% N3=length(initiala);
% N4=length(initialr);
% N=N1*N2*N3*N4;
% m = length(beta);
% P = zeros(N,3*(m)+5+m);  % initialize the matrix
% P_current = zeros(1,3*(m)+5+m)
% % P = zeros(m,3);  % initialize the matrix
%YE GARZ-----------------------------------------------------------------

%YE CGARZ----------------------------------------------------------------
for iter1 = 1:N1
    for iter2 = 1:N2
        
            count=count+1;
            
            
            x01 = [initials(iter1) initialm(iter2)];


            
            for j = 1:m
               
                if j~=7 && j~=1 && j~=9
                    %        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
                    %        options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
                    %options = optimset('MaxFunEvals', 1*10^10);
                    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rho_m)...
                        -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rho_m)...
                        -flux),0*flux)))^2,x01,[],[],[],[],lb1,ub1);%CGARZ!!!
                    
                    %THIS ONE!!! GARZ!!!
                    %-----------------------
                    % triangular flux function
                    %-----------------------
                    % % triangle, 2 parameters
                    % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig(density,x(1),x(2),rhom)...
                    %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_trig(density,x(1),x(2),rhom)...
                    %    -flux),0*flux))^2,x0);
                    % P(j,:) = x;
                    %---------------------------
                    % Daganzo--Newell model
                    % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig_1(density,x(1),x(2),rhom)...
                    %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_trig_1(density,x(1),x(2),rhom)...
                    %    -flux),0*flux))^2,x0);
                    %------------------------------------------------
                    % 3-p arctangent function
                    %    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
                    %    -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
                    %    -flux),0)))^2,x0);
                    
                    %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
                    %    -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
                    %    -flux),0)))^2,x0,[],[],[],[],lb,ub);
                    %---------------------------
                    % 2-p collapsed model
                    %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
                    %    -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
                    %    -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);%CGARZ!!!
                    % ----------------------------
                    % % constrained minimization problem
                    %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,uf,rhof,rhom)...
                    %    -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,uf,rhof,rhom)...
                    %    -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);
                    %---------------------------
                    
                    % % determine the threshold cretical density, from free to conjected
                    %    indx1 = find(density>x);
                    %    D1 = density(indx1);
                    %    Q1 = flux(indx1);
                    %    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_det_rhof(D1,x,uf,rhof,rhom)...
                    %    -Q1),0)))^2+(1-beta(j))*(norm(max(-(func_det_rhof(D1,x,uf,rhof,rhom)...
                    %    -Q1),0)))^2,x0);
                    %----------------------------
                    %------------------------
                    % best fitting curve respect to the free data/cgarz model
                    % a parameter quadratic form
                    % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(D_free,x(1),x(2))...
                    %    -Q_free),0*Q_free))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(D_free,x(1),x(2))...
                    %    -Q_free),0*Q_free))^2,x0);
                    %-----------------------
                    % quadratic form flux function
                    %[x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(density,x,rhom)...
                    %   -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(density,x,rhom)...
                    %   -flux),0*flux))^2,x0);
                    %--------------------------------
                    % greenshield phase transition model
                    % first determine the standard curve
                    % [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(density,0,x(1),x(2),rhom)...
                    %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(density,0,x(1),x(2),rhom)...
                    %    -flux),0*flux))^2,x0);
                    %--------------------------------
                    % determine a family of curves, parametrized by beta---q(link to q)
                    % [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(density,x,rhoc,vm,rhom)...
                    %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(density,x,rhoc,vm,rhom)...
                    %    -flux),0*flux))^2,x0);
                    %---------------------------------
                    % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_piecewise(D,x,rhom0,vm)...
                    %    -Q),0*Q))^2+(1-beta(j))*norm(max(-(traffic_flux_piecewise(D,x,rhom0,vm)...
                    %    -Q),0*Q))^2,x0);
                    
                    %   P(j,:) = x
                    
                    P2(count,(j-1)*3+1) = x(1);
                    P2(count,(j-1)*3+2) = x(2);
                    P2(count,(j-1)*3+3)=fval;
                    P2(count,28)=initials(iter1);
                    P2(count,29)=initialm(iter2);
                    count
                    j
                end
          
        end
    end
    
    
    
end
%YE CGARZ-------------------------------------------------------------



% fd_equi = P;
% rhoc = rho_critical(fd_equi(1),fd_equi(2),rhom)

% data_fd_ngsim_p3 have 3 parameters
% savefile = sprintf('data_fd_pt_minn_%1.0f',rhom);
%strsv = ['rhom0_',num2str(rhom0),'sigma0_',num2str(sigma0),'mu0_',num2str(mu0)]
savefile = sprintf('CGARZ_varyingIC-realistic_rhom0_fval-temp-updated');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'P2');    % P is the parameters, n cross 3

% load data_fd_minn_arctan.mat


%\\\\\\\\\\\\\\\\\\
%  The plots part
%//////////////////

figure(1)    % 

grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
% plot(d_s1,f_s1,'r.'), hold on
% plot(d_s3,f_s3,'b*'), hold on

%---------------------------------------------------------
%  4 parameters
%  plot(grid,traffic_flux_smooth(grid,P(1,1),P(1,2),P(1,3),rho_m),'b-','linewidth',2.5),hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,4),P(1,5),P(1,6),rho_m),'m-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,7),P(1,8),P(1,9),rho_m),'r-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,10),P(1,11),P(1,12),rho_m),'k-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,13),P(1,14),P(1,15),rho_m),'r-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,16),P(1,17),P(1,18),rho_m),'m-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,19),P(1,20),P(1,21),rho_m),'b-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,22),P(1,23),P(1,24),rho_m),'k-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(1,25),P(1,26),P(1,27),rho_m),'b--','linewidth',2.5), hold on 
% % plot(grid,traffic_flux_smooth(grid,P(10,1),P(10,2),P(10,3),rhom),'k--','linewidth',2.5), hold on  
% % hold off
% 
% %  3 parameters
% %  plot(grid,traffic_flux_smooth(grid,P(1,1),P(1,2),P(1,3),rhom),'b-','linewidth',2.5),hold on
% %  plot(grid,traffic_flux_smooth(grid,P(2,1),P(2,2),P(2,3),rhom),'m-','linewidth',2.5), hold on
% %  plot(grid,traffic_flux_smooth(grid,P(3,1),P(3,2),P(3,3),rhom),'r-','linewidth',2.5), hold on
% %  plot(grid,traffic_flux_smooth(grid,P(4,1),P(4,2),P(4,3),rhom),'k-','linewidth',2.5), hold on
% %  plot(grid,traffic_flux_smooth(grid,P(5,1),P(5,2),P(5,3),rhom),'r-.','linewidth',2.5), hold on
% %  plot(grid,traffic_flux_smooth(grid,P(6,1),P(6,2),P(6,3),rhom),'m-.','linewidth',2.5), hold on
% %  plot(grid,traffic_flux_smooth(grid,P(7,1),P(7,2),P(7,3),rhom),'b-.','linewidth',2.5), hold on
% %  plot(grid,traffic_flux_smooth(grid,P(8,1),P(8,2),P(8,3),rhom),'k-.','linewidth',2.5), hold on
% %  plot(grid,traffic_flux_smooth(grid,P(9,1),P(9,2),P(9,3),rhom),'b--','linewidth',2.5), hold on 
% % % plot(grid,traffic_flux_smooth(grid,P(10,1),P(10,2),P(10,3),rhom),'k--','linewidth',2.5), hold on  
% %  hold off
% %--------------------------------------------------------
% % % 2 parameters
% % plot(grid,traffic_flux_smooth_2(grid,P(1,1),P(1,2),rho_m),'b-','linewidth',2.5),hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(2,1),P(2,2),rho_m),'r-','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(3,1),P(3,2),rho_m),'k-','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(4,1),P(4,2),rho_m),'r-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(5,1),P(5,2),rho_m),'m-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(6,1),P(6,2),rho_m),'b-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(7,1),P(7,2),rho_m),'k-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(8,1),P(8,2),rho_m),'b--','linewidth',2.5), hold on 
% % plot(grid,traffic_flux_smooth_2(grid,P(9,1),P(9,2),rho_m),'m-','linewidth',2.5), hold on
% % % plot(grid,traffic_flux_smooth_2(grid,P(end,1),P(end,2),rhom),'k--','linewidth',2.5), hold on  
% %hold off
% 
% % % 1 parameters
% % plot(grid,traffic_flux_smooth_2(grid,P(1),p,rhom),'b-','linewidth',2.5),hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(2),p,rhom),'r-','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(3),p,rhom),'k-','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(4),p,rhom),'r-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(5),p,rhom),'m-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(6),p,rhom),'b-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(7),p,rhom),'k-.','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(8),p,rhom),'b--','linewidth',2.5), hold on 
% % plot(grid,traffic_flux_smooth_2(grid,P(9),p,rhom),'m-','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(end),p,rhom),'k--','linewidth',2.5), hold on  
% % hold off
%--------------------------------------------------------
%plot(grid,func_fd_seibold_2p(grid,P1(1,1),P1(1,2),rho0,v0,f0,uf,rhof,P1(1,22)),'r--','linewidth',4),hold on
plot(grid,func_fd_seibold_2p_updated(grid,P2(2,4),P2(2,5),rho0,uf,rhof,rho_m),'k-','linewidth',2), hold on
plot(grid,func_fd_seibold_2p_updated(grid,P2(2,7),P2(2,8),rho0,uf,rhof,rho_m),'k-','linewidth',2), hold on
plot(grid,func_fd_seibold_2p_updated(grid,P2(2,10),P2(2,11),rho0,uf,rhof,rho_m),'k-','linewidth',2), hold on
plot(grid,func_fd_seibold_2p_updated(grid,P2(2,13),P2(2,14),rho0,uf,rhof,rho_m),'k-','linewidth',2), hold on
plot(grid,func_fd_seibold_2p_updated(grid,P2(2,16),P2(2,17),rho0,uf,rhof,rho_m),'k-','linewidth',2), hold on
%plot(grid,func_fd_seibold_2p(grid,P1(1,13),P1(1,14),rho0,v0,f0,uf,rhof,P1(1,22)),'m-','linewidth',4), hold on
plot(grid,func_fd_seibold_2p_updated(grid,P2(2,22),P2(2,23),rho0,uf,rhof,rho_m),'k-','linewidth',2), hold on
%plot(grid,func_fd_seibold_2p(grid,P1(1,17),P1(1,18),rho0,v0,f0,uf,rhof,P1(1,22)),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(end,1),P(end,2),rho0,v0,f0,uf,rhof,rhom),'b-.','linewidth',4),
hold off
% %-------------------------------------------------------
% % 3- arctangent function
% % plot(grid,func_fd_seibold_3p(grid,P(1,1),P(1,2),P(1,3),rhom),'k--','linewidth',4),hold on
% % plot(grid,func_fd_seibold_3p(grid,P(2,1),P(2,2),P(2,3),rhom),'k-','linewidth',4), hold on
% % plot(grid,func_fd_seibold_3p(grid,P(3,1),P(3,2),P(3,3),rhom),'k-','linewidth',4), hold on
% % plot(grid,func_fd_seibold_3p(grid,P(4,1),P(4,2),P(4,3),rhom),'k-','linewidth',4), hold on
% % plot(grid,func_fd_seibold_3p(grid,P(5,1),P(5,2),P(5,3),rhom),'k-','linewidth',4), hold on
% % plot(grid,func_fd_seibold_3p(grid,P(6,1),P(6,2),P(6,3),rhom),'k-','linewidth',4), hold on
% % plot(grid,func_fd_seibold_3p(grid,P(7,1),P(7,2),P(7,3),rhom),'k-','linewidth',4), hold on
% % plot(grid,func_fd_seibold_3p(grid,P(8,1),P(8,2),P(8,3),rhom),'k-','linewidth',4), hold on
% % plot(grid,func_fd_seibold_3p(grid,P(9,1),P(9,2),P(9,3),rhom),'k-','linewidth',4), hold on 
% % plot(grid,func_fd_seibold_3p(grid,P(end,1),P(end,2),P(end,3),rhom),'k-.','linewidth',4),
% % hold off
% %--------------------------------------------------------
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(1),rhom),'r--','linewidth',4),hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(2),rhom),'k-','linewidth',4), hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(3),rhom),'k-','linewidth',4), hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(4),rhom),'k-','linewidth',4), hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(5),rhom),'k-','linewidth',4), hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(6),rhom),'k-','linewidth',4), hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(7),rhom),'k-','linewidth',4), hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(8),rhom),'k-','linewidth',4), hold on
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(9),rhom),'k-','linewidth',4), hold on 
% % plot(grid,traffic_flux_trig(grid,vm,rhof,P(end),rhom),'b-.','linewidth',4),
% % hold off
% %--------------------------------------------------------
% % plot(grid,traffic_flux_trig_1(grid,P(1,1),P(1,2),rhom),'b--','linewidth',4),hold on
% %--------------------------------------------------------
% % plot(grid,func_flux_creen_pt(grid,0,P(1,1),P(1,2),rhom),'b--','linewidth',4),hold on
% %--------------------------------------------------------
% %plot(grid,func_flux_creen_pt(grid,P(1),rhoc,vm,rhom),'b--','linewidth',4),hold on
% %plot(grid,func_flux_creen_pt(grid,P(2),rhoc,vm,rhom),'k-','linewidth',4),hold on
% 
% 
% axis([0 rho_m+25 0 11000])
% % title('The G-AR model','fontsize',14)
% title('GARZ','fontsize',14)%  title('Fitting to Smooth FDs W.R.T different \gamma','fontsize',14)
% % legend('\omega = 0.01','\omega = 0.1','\omega = 0.5',...
% %     '\omega = 0.9','\omega = 0.99')
% %  legend('Data','\gamma = 0.001','\gamma = 0.01','\gamma = 0.1','\gamma = 0.5',...
% %      '\gamma = 0.9','\gamma = 0.99','\gamma = 0.999')
% %legend('Data','FD')
% %GARZ-----------------------------------------------
% % strl = ['\lambda_0 = ',num2str(lambda0)];
% % strp = ['p_0 = ',num2str(p0)]
% % stra = ['\alpha_0 = ',num2str(alpha0)]
% % strr = ['\rho_{m0} = ',num2str(rhom0)]
% % text(rhom-10,-200,'\rho_{max}','fontsize',14)
% % text(rhom/3,2000,'\Omega_{f}','fontsize',18)
% % text(rhom/2,8000,'\Omega_{c}','fontsize',18)
% % text(rhom*0.6,8000,strl,'fontsize',14)
% % text(rhom*0.6,9000,strp,'fontsize',14)
% % text(rhom*0.6,10000,stra,'fontsize',14)
% % text(rhom*0.6,7000,strr,'fontsize',14)
% % text(rho_m-10,-200,num2str(rho_m),'fontsize',14)
% %GARZ-----------------------------------------------
% 
% %CGARZ-----------------------------------------------
% % strs = ['\sigma_0 = ',num2str(sigma0)];
% % strm = ['\mu_0 = ',num2str(mu0)]
% % strr = ['\rho_{m0} = ',num2str(rhom0)]
% % text(rhom/3,2000,'\Omega_{f}','fontsize',18)
% % text(rhom/2,8000,'\Omega_{c}','fontsize',18)
% % text(rho_m*0.85,8000,strs,'fontsize',14)
% % text(rho_m*0.85,9000,strm,'fontsize',14)
% % text(rho_m*0.85,7000,strr,'fontsize',14)
% text(P(1,22)-25,-200,num2str(P(1,22)),'fontsize',14)
% %CGARZ-----------------------------------------------
% 
% set(gca,'linewidth',2)
% set(gca,'xtick',[])
% % set(gca,'ytick',[])
% xlabel('density \rho','fontsize',14)
% ylabel('flow rate Q (veh/h)','fontsize',14)
% set(gca,'fontsize',14)
% res = 600;
% set(gcf,'paperpositionmode','auto')
% set(gcf,'position',[10 50 res*1.25 res])
% set(gca,'position',[.1 .06 .88 .88])


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

% free flow branch
function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);

function y = diff_fd_free(rho,uf,rhof)

y = uf*(1-2*rho/rhof);

function y = linear_bound(rho,v0,rho0,f0)

y = f0+v0*(rho-rho0);

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

