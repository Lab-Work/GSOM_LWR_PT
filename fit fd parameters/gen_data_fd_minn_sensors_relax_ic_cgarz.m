function P = gen_data_fd_minn_sensors_relax_ic_cgarz(rhom)

if nargin<1 rhom = 533; end
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

rho0 = 36; uf = 95;  rhof = 555;
v0 = diff_fd_free(rho0,uf,rhof);
f0 = fd_free(rho0,uf,rhof);

lb = [1.0e-03, rho0,rho0];
ub = [1000,2000,2000]; 
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

count=0;

%YE CGARZ-------------------------------------------------------
initials = [1,10,50,100,500,1000,4000];
initialm = [1,10,50,100,500,1000,4000];
initialr = [1,10,50,100,500,1000,4000];
N1=length(initials);
N2=length(initialm);
N3=length(initialr);
N=N1*N2*N3;
m = length(beta);
P = zeros(N,6);  % initialize the matrix
% P = zeros(m,3);  % initialize the matrix
%YE CGARZ-------------------------------------------------------

%YE GARZ-----------------------------------------------------------------
% initiall = [1,10,50,100,500,1000,4000];
% initialp = [0.1,0.5,0.9];
% initiala = [1,10,50,100,500,1000,4000];
% initialr = [1,50,100,500,1000,4000];
% lbeq = [0,0,0,0];
% lb1 = [0,0,0];
% N1=length(initiall);
% N2=length(initialp);
% N3=length(initiala);
% N4=length(initialr);
% N=N1*N2*N3*N4;
% m = length(beta);
% P = zeros(N,3*(m)+5);  % initialize the matrix
% P_current = zeros(1,3*(m)+5)
% % P = zeros(m,3);  % initialize the matrix
%YE GARZ-----------------------------------------------------------------

%YE CGARZ----------------------------------------------------------------
for iter1 = 1:N1
    for iter2 = 1:N2
        for iter3 = 1:N3
            count=count+1;
            
            x0 = [initials(iter1) initialm(iter2) initialr(iter3)];
            %options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5,);
            %options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
            %options = optimset('MaxFunEvals', 1*10^10);
            [x,fval] = fmincon(@(x) beta(7)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,x(3))...
                -flux),0*flux)))^2+(1-beta(7))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,x(3))...
                -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);%CGARZ!!!
            
            P(count,4) = initials(iter1);
            P(count,5) = initialm(iter2);
            P(count,6) = initialm(iter3);
            P(count,1) = x(1);
            P(count,2) = x(2);
            P(count,3) = x(3);
            rho_m = x(3);
            count
            rho_m

            
%             for j = 1:m
%                 lbeq_1 = [0,0];
%                 if j~=7
%                     %        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
%                     %        options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
%                     %options = optimset('MaxFunEvals', 1*10^10);
%                     [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rho_m)...
%                         -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rho_m)...
%                         -flux),0*flux)))^2,x01,[],[],[],[],lb1,ub1);%CGARZ!!!
%                     
%                     %THIS ONE!!! GARZ!!!
%                     %-----------------------
%                     % triangular flux function
%                     %-----------------------
%                     % % triangle, 2 parameters
%                     % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig(density,x(1),x(2),rhom)...
%                     %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_trig(density,x(1),x(2),rhom)...
%                     %    -flux),0*flux))^2,x0);
%                     % P(j,:) = x;
%                     %---------------------------
%                     % Daganzo--Newell model
%                     % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig_1(density,x(1),x(2),rhom)...
%                     %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_trig_1(density,x(1),x(2),rhom)...
%                     %    -flux),0*flux))^2,x0);
%                     %------------------------------------------------
%                     % 3-p arctangent function
%                     %    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                     %    -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                     %    -flux),0)))^2,x0);
%                     
%                     %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                     %    -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                     %    -flux),0)))^2,x0,[],[],[],[],lb,ub);
%                     %---------------------------
%                     % 2-p collapsed model
%                     %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
%                     %    -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
%                     %    -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);%CGARZ!!!
%                     % ----------------------------
%                     % % constrained minimization problem
%                     %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,uf,rhof,rhom)...
%                     %    -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,uf,rhof,rhom)...
%                     %    -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);
%                     %---------------------------
%                     
%                     % % determine the threshold cretical density, from free to conjected
%                     %    indx1 = find(density>x);
%                     %    D1 = density(indx1);
%                     %    Q1 = flux(indx1);
%                     %    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_det_rhof(D1,x,uf,rhof,rhom)...
%                     %    -Q1),0)))^2+(1-beta(j))*(norm(max(-(func_det_rhof(D1,x,uf,rhof,rhom)...
%                     %    -Q1),0)))^2,x0);
%                     %----------------------------
%                     %------------------------
%                     % best fitting curve respect to the free data/cgarz model
%                     % a parameter quadratic form
%                     % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(D_free,x(1),x(2))...
%                     %    -Q_free),0*Q_free))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(D_free,x(1),x(2))...
%                     %    -Q_free),0*Q_free))^2,x0);
%                     %-----------------------
%                     % quadratic form flux function
%                     %[x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(density,x,rhom)...
%                     %   -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(density,x,rhom)...
%                     %   -flux),0*flux))^2,x0);
%                     %--------------------------------
%                     % greenshield phase transition model
%                     % first determine the standard curve
%                     % [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(density,0,x(1),x(2),rhom)...
%                     %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(density,0,x(1),x(2),rhom)...
%                     %    -flux),0*flux))^2,x0);
%                     %--------------------------------
%                     % determine a family of curves, parametrized by beta---q(link to q)
%                     % [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(density,x,rhoc,vm,rhom)...
%                     %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(density,x,rhoc,vm,rhom)...
%                     %    -flux),0*flux))^2,x0);
%                     %---------------------------------
%                     % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_piecewise(D,x,rhom0,vm)...
%                     %    -Q),0*Q))^2+(1-beta(j))*norm(max(-(traffic_flux_piecewise(D,x,rhom0,vm)...
%                     %    -Q),0*Q))^2,x0);
%                     
%                     %   P(j,:) = x
%                     
%                     P(count,(j-1)*2+1) = x(1);
%                     P(count,(j-1)*2+2) = x(2);
%                     count
%                     j
%                 end
%             end
        end
    end
    
    
    
end
%YE CGARZ-------------------------------------------------------------


% for iter1 = 1:N1
%     for iter2 = 1:N2
%         for iter3 = 1:N3
%             for iter4 =1:N4
%                 count=count+1;
%                 
%                 x0 = [initiall(iter1) initialp(iter2) initiala(iter3) initialr(iter4)];
%                 x01 = [initiall(iter1) initialp(iter2) initiala(iter3)];
%                 %options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5,);
%                 %options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
%                 %options = optimset('MaxFunEvals', 1*10^10);
%                 [x,fval] = fmincon(@(x)beta(7)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
%                     -flux),0*flux))^2+(1-beta(7))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
%                     -flux),0*flux))^2,x0,[],[],[],[],lbeq,[],[]);
%                 
%                 P(count,7*3-2) = x(1);
%                 P(count,7*3-1) = x(2);
%                 P(count,7*3) = x(3);
%                 P(count,m*3+5) = x(4);
%                 P(count,m*3+1) = initiall(iter1);
%                 P(count,m*3+2) = initialp(iter2);
%                 P(count,m*3+3) = initiala(iter3);
%                 P(count,m*3+4) = initialr(iter4);
%                 rho_m = x(4);
%                 rho_m
%                 
%                 
%                 for j = 1:m
%                     if j~=7
%                         %        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
%                         %        options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
%                         %options = optimset('MaxFunEvals', 1*10^10);
%                         [x,fval] = fmincon(@(x)beta(j)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),rho_m)...
%                             -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),rho_m)...
%                             -flux),0*flux))^2,x01,[],[],[],[],lb1,[],[]);
%                         
%                         %THIS ONE!!! GARZ!!!
%                         %-----------------------
%                         % triangular flux function
%                         %-----------------------
%                         % % triangle, 2 parameters
%                         % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig(density,x(1),x(2),rhom)...
%                         %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_trig(density,x(1),x(2),rhom)...
%                         %    -flux),0*flux))^2,x0);
%                         % P(j,:) = x;
%                         %---------------------------
%                         % Daganzo--Newell model
%                         % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_trig_1(density,x(1),x(2),rhom)...
%                         %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_trig_1(density,x(1),x(2),rhom)...
%                         %    -flux),0*flux))^2,x0);
%                         %------------------------------------------------
%                         % 3-p arctangent function
%                         %    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                         %    -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                         %    -flux),0)))^2,x0);
%                         
%                         %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                         %    -flux),0)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_3p(density,x(1),x(2),x(3),rhom)...
%                         %    -flux),0)))^2,x0,[],[],[],[],lb,ub);
%                         %---------------------------
%                         % 2-p collapsed model
%                         %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
%                         %    -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,v0,f0,uf,rhof,rhom)...
%                         %    -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);%CGARZ!!!
%                         % ----------------------------
%                         % % constrained minimization problem
%                         %    [x,fval] = fmincon(@(x) beta(j)*(norm(max((func_fd_seibold_2p(density,x(1),x(2),rho0,uf,rhof,rhom)...
%                         %    -flux),0*flux)))^2+(1-beta(j))*(norm(max(-(func_fd_seibold_2p(density,x(1),x(2),rho0,uf,rhof,rhom)...
%                         %    -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);
%                         %---------------------------
%                         
%                         % % determine the threshold cretical density, from free to conjected
%                         %    indx1 = find(density>x);
%                         %    D1 = density(indx1);
%                         %    Q1 = flux(indx1);
%                         %    [x,fval] = fminsearch(@(x) beta(j)*(norm(max((func_det_rhof(D1,x,uf,rhof,rhom)...
%                         %    -Q1),0)))^2+(1-beta(j))*(norm(max(-(func_det_rhof(D1,x,uf,rhof,rhom)...
%                         %    -Q1),0)))^2,x0);
%                         %----------------------------
%                         %------------------------
%                         % best fitting curve respect to the free data/cgarz model
%                         % a parameter quadratic form
%                         % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(D_free,x(1),x(2))...
%                         %    -Q_free),0*Q_free))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(D_free,x(1),x(2))...
%                         %    -Q_free),0*Q_free))^2,x0);
%                         %-----------------------
%                         % quadratic form flux function
%                         %[x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_qrd(density,x,rhom)...
%                         %   -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_qrd(density,x,rhom)...
%                         %   -flux),0*flux))^2,x0);
%                         %--------------------------------
%                         % greenshield phase transition model
%                         % first determine the standard curve
%                         % [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(density,0,x(1),x(2),rhom)...
%                         %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(density,0,x(1),x(2),rhom)...
%                         %    -flux),0*flux))^2,x0);
%                         %--------------------------------
%                         % determine a family of curves, parametrized by beta---q(link to q)
%                         % [x,fval] = fminsearch(@(x)beta(j)*norm(max((func_flux_creen_pt(density,x,rhoc,vm,rhom)...
%                         %    -flux),0*flux))^2+(1-beta(j))*norm(max(-(func_flux_creen_pt(density,x,rhoc,vm,rhom)...
%                         %    -flux),0*flux))^2,x0);
%                         %---------------------------------
%                         % [x,fval] = fminsearch(@(x)beta(j)*norm(max((traffic_flux_piecewise(D,x,rhom0,vm)...
%                         %    -Q),0*Q))^2+(1-beta(j))*norm(max(-(traffic_flux_piecewise(D,x,rhom0,vm)...
%                         %    -Q),0*Q))^2,x0);
%                         
%                         %   P(j,:) = x
%                         
%                         P(count,(j-1)*3+1) = x(1);
%                         P(count,(j-1)*3+2) = x(2);
%                         P(count,(j-1)*3+3) = x(3);
%                         P_current=P(count,:)
%                         strsv = ['GARZ_varyingIC',num2str(count)]
%                         savefile = sprintf(strsv);
%                         save(savefile,'P_current'); 
%                         count
%                         j
%                     end
%                 end
%             end
%         end
%     end
%     
%     
%     
% end

% fd_equi = P;
% rhoc = rho_critical(fd_equi(1),fd_equi(2),rhom)

% data_fd_ngsim_p3 have 3 parameters
% savefile = sprintf('data_fd_pt_minn_%1.0f',rhom);
%strsv = ['rhom0_',num2str(rhom0),'sigma0_',num2str(sigma0),'mu0_',num2str(mu0)]
savefile = sprintf('GARZ_varyingIC_eq');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'P');    % P is the parameters, n cross 3

% load data_fd_minn_arctan.mat


%\\\\\\\\\\\\\\\\\\
%  The plots part
%//////////////////


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
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(1),rhom),'r--','linewidth',4),hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(2),rhom),'k-','linewidth',4), hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(3),rhom),'k-','linewidth',4), hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(4),rhom),'k-','linewidth',4), hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(5),rhom),'k-','linewidth',4), hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(6),rhom),'k-','linewidth',4), hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(7),rhom),'k-','linewidth',4), hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(8),rhom),'k-','linewidth',4), hold on
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(9),rhom),'k-','linewidth',4), hold on 
% plot(grid,traffic_flux_trig(grid,vm,rhof,P(end),rhom),'b-.','linewidth',4),
% hold off
%--------------------------------------------------------
% plot(grid,traffic_flux_trig_1(grid,P(1,1),P(1,2),rhom),'b--','linewidth',4),hold on
%--------------------------------------------------------
% plot(grid,func_flux_creen_pt(grid,0,P(1,1),P(1,2),rhom),'b--','linewidth',4),hold on
%--------------------------------------------------------
%plot(grid,func_flux_creen_pt(grid,P(1),rhoc,vm,rhom),'b--','linewidth',4),hold on
%plot(grid,func_flux_creen_pt(grid,P(2),rhoc,vm,rhom),'k-','linewidth',4),hold on





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