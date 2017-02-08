function [e_rho,e_vel]=comput_minn_validation_cluster(day)
% if nargin<2 models = 2; resol = 1; proj = 1; end
if nargin<1 day = 14; end
clc;
resol = 7;
proj = 1;%or 2
   tic
%==========================================================================
% We do Validation respect different model for more days from 4pm to 5pm. 
% Scalar model, AR model, and our generalized AR model. We tested totally
% 76 days of year 2003 of Minneapolis data.
% Shimao Fan, March 13 2012
% Temple University.
% July 31, 2012.
% j specity different days
%****************************
% switch models 
%   case 1 --- garz model
%   case 2 --- arz model
%   case 3 --- arzq model
%   case 4 --- lwr model
%   case 5 --- lwrq model
%   case 6 --- interpolation, no model applied here.
%****************************
% this code can perform two test, detector counts the number of sensor
% stations of the test case.
% if detector == 3, it is the original 3-sensor tests
% if detector == 6, it is the new 6-sensor tests.
%****************************
% parameter == 1    % model with 3 free parameters
% parameter == 2    % model with 2 free parameters
%****************************
% finaly cleaned on Jan 24 2013
%==========================================================================
% clc;
close all;
%----------------------------------------------------
% same fixed parameters
%rhom = 533;         % maximum density
%----------------------------------------------------
% Load data in the three test sensors
% density data
day
load data_dens_minn_2.mat             % sensor 2 data from 4pm tp 5pm/D
D = Dens;
load data_dens_minn_1.mat             % sensor 1, inflow density/LBD
LBD = Dens;
load data_dens_minn_3.mat             % sensor 3, outflow density/RBD
RBD = Dens;
% velocity data
load data_vel_minn_1.mat              % sensor 1, inflow density/LBV
LBV = Vel;
load data_vel_minn_3.mat              % sensor 3, outflow density/RBV
RBV = Vel;
load data_vel_minn_2.mat              % sensor 2 data from 4pm tp 5pm/Vel

days = size(Vel,1);                   % number of days
%-----------------------------------------
rho_thresh = 20;   % 200 meter apart
indexx = find(D>rho_thresh);
%---------------------------------------
D_remov = D(indexx);
V_remov = Vel(indexx);
data_rho_remov = sort(reshape(D_remov,1,[]),'descend');
data_vel_remov = sort(reshape(V_remov,1,[]),'descend');
likehood_up = ceil(.001*length(data_vel_remov));   % 99%
likehood_low = ceil(.999*length(data_vel_remov));
delta_rho = data_rho_remov(likehood_up);
ess_vel(1) = data_vel_remov(likehood_up);
ess_vel(2) = data_vel_remov(likehood_low);
delta_vel = ess_vel(1)-ess_vel(2);
%REMOVE DENSITY VERY SMALL!!!
%-----------------------------------------------------
% data selection for congested days
% density_level = 80;
% n = ceil(size(D,2)/2);    % select 4pm-5pm
% max_rho = mean(D(1:days,1:n)');  % return max density for each day
% % low density data at sensor station 2
% % indx = find(max_rho<density_level);
% % high density data at sensor 2
% indx = find(max_rho>=density_level);
% D = D(indx,:);
% LBD = LBD(indx,:);
% RBD = RBD(indx,:);
% LBV = LBV(indx,:);
% RBV = RBV(indx,:);
% Vel = Vel(indx,:);
%------------------------------------------------------
% % make selection for special patterns, jump from free to congested
% load data_minn_pattern.mat   % in matrix shimao
% indx = shimao(:,1);
% ts = 30;   % 20 data points, 10 minutes
% dc = shimao(:,2)+6;  % the critial points where data jumps
% D = D(indx,:);
% LBD = LBD(indx,:);
% RBD = RBD(indx,:);
% LBV = LBV(indx,:);
% RBV = RBV(indx,:);
% Vel = Vel(indx,:);
% days = size(D,1);  % so far, we have the correct days
%------------------------------------------
% then select the correct time interval
% for j = 1:days   % number of days in analysis/20 minutes interval
%     D1(j,:) = D(j,dc(j)-ts:dc(j)+ts);
%     LBD1(j,:) = LBD(j,dc(j)-ts:dc(j)+ts);
%     RBD1(j,:) = RBD(j,dc(j)-ts:dc(j)+ts);
% 
%     RBV1(j,:) = RBV(j,dc(j)-ts:dc(j)+ts);
%     LBV1(j,:) = LBV(j,dc(j)-ts:dc(j)+ts);
%     Vel1(j,:) = Vel(j,dc(j)-ts:dc(j)+ts);
% end
% at this point, we have selected data in approperiate patterns in D1,... 
%--------------------------------------------------------
% load fd data from fitting process
% garz with arctangent form flux function

% load data_fd_minn_arctan.mat
% fd_arctan = P;
%-------------------------------
% load data_fd_minn_arctan_2p.mat
% load data_fd_minn_cgarz_38.mat
%load data_fd_minn_cgarz_34.mat
load CGARZ_100beta_2para_updated.mat

fd_collap = P2;
% load fd data for phase transition model
load data_fd_pt_minn_533.mat   % fitted at the sensor 2 data.
fd_pt = P';   % one parameter about rhoc, critical density
%-------------------------------
% load data_fd_minn_test.mat
%load data_fd_minn_533.mat 
load GARZ_para-0406.mat % fitted at the sensor 2 data.

% load data_fd_rhom_minn_trb.mat
%**************************************************************
%**** calculate fitting pramaters for the Minneapolis data  ***
%**************************************************************
% calculate maximum velocity
% first locate the equilibrium curve.
midd = ceil(size(para,1)/2);         % the equilibium curve
equi_lda = para(midd,1);
equi_p = para(midd,2);
equi_alpha = para(midd,3);
fd_equi = [equi_lda,equi_p,equi_alpha];

load GARZ_rhom.mat;
rhom=rhom_garz(1,1);
vm = u_max(fd_equi(1),fd_equi(2),fd_equi(3),rhom); % 93.3846 % maximum velocity
% rhoc = rho_critical(equi_lda,equi_p,rhom);
% q_max = traffic_flux_smooth(rhoc,equi_lda,equi_p,equi_alpha,rhom);
% vm_q = 4*q_max/rhom;
% vm_q = vm;
%-----------------------------------------------
% find out empty road velocity w   GARZ
lambda = para(:,1);
p = para(:,2);
alpha = para(:,3);
w = func_W(lambda,p,alpha,rhom);  % values of w(\omega)
w_min = min(w);
w_max = max(w);

%\\\\\\\\\\\\\\\\\\\\\\\
% fitting fd parameters
%///////////////////////
order = 2;

cl = polyfit(w(1:end),lambda(1:end),order);
cp = polyfit(w(1:end),p(1:end),order);
ca = polyfit(w(1:end),alpha(1:end),order);

% cl = polyfit(w(2:end),lambda(2:end),order);
% cp = polyfit(w(2:end),p(2:end),order);
% ca = polyfit(w(2:end),alpha(2:end),order);
%------------------------------------------------
% % 3-p arctangent function
% w_arctan = func_w_arctan(fd_arctan(:,1),fd_arctan(:,2),fd_arctan(:,3),rhom);
% 
% % c_sigma= polyfit(w_arctan,fd_arctan(:,1),order-1);
% % c_mu = polyfit(w_arctan,fd_arctan(:,2),order-1);
% % c_alpha = polyfit(w_arctan,fd_arctan(:,3),order-1);
% 
% c_sigma= polyfit(w_arctan(2:end-1),fd_arctan(2:end-1,1),order);
% c_mu = polyfit(w_arctan(2:end-1),fd_arctan(2:end-1,2),order);
% c_alpha = polyfit(w_arctan(2:end-1),fd_arctan(2:end-1,3),order);
% 
% w_min_arc = min(w_arctan);
% w_max_arc = max(w_arctan);
%---------------------------------------------
% fitting with respect 2-parameter arctangent flux
% rho0 = 38.5; 
% rho0 = 38;

load CGARZ_eq_6para_updated.mat;
rho0 = para_eq(3,1);
% uf = 114; 
uf = para_eq(4,1);
rhof = para_eq(5,1);
rhom_cgarz=para_eq(6,1);
% rhof = 1026;
f0 = fd_free(rho0,uf,rhof);
v0 = diff_fd_free(rho0,uf,rhof);
q = func_q(fd_collap(:,1),fd_collap(:,2),rho0,uf,rhof,rhom_cgarz);

% ca_2p= polyfit(q,fd_collap(:,1),order-1);
% cb_2p = polyfit(q,fd_collap(:,2),order-1);



ca_2p= polyfit(q(1:end-1),fd_collap(1:end-1,1),order-1);
cb_2p = polyfit(q(1:end-1),fd_collap(1:end-1,2),order-1);

midd_arc = ceil(size(fd_collap,1)/2);

q_eq_arc = q(midd_arc);   % equilibrium cuve
q_min_arc = min(q);
q_max_arc = max(q);
%-----------------------------------------------
% phase transition model
% beta = linspace(.01,1-.01,100);
% ct = polyfit(beta,fd_pt,order-1);

%-----------------------------------------------
dist = [0.616 0.598 0.72 0.6975 0.6975];
pos = [0,dist(1),sum(dist(1:2))];
%-----------------------------------------------
vec_h = .032./2.^(0:6);   % in km.
h = vec_h(1);

%*********************************
%****  Numerical solvers  ********
%*********************************
e_rho = zeros(1,8);
e_vel = zeros(1,8);
% e = zeros(1,days);
%------------------------------------------------
cost = 0;    % compute the time used for the program

tfinal = 1;

i = day;
% tfinal = 1/2;   % 1/3 hour, 20 minutes

% for i = 1:days
%    tic
   
%    B(1,:) = LBD1(i,:);
%    B(2,:) = D1(i,:);
%    B(3,:) = RBD1(i,:);
%    V(1,:) = LBV1(i,:);
%    V(2,:) = Vel1(i,:);
%    V(3,:) = RBV1(i,:);
   %----------------------
   B(1,:) = LBD(i,:);
   B(2,:) = D(i,:);
   B(3,:) = RBD(i,:);
   V(1,:) = LBV(i,:);
   V(2,:) = Vel(i,:);
   V(3,:) = RBV(i,:);
   %--------------------------------------------  
%    switch models
%        case 1    % garz model
%            if parameter == 1
            % cgarz
%            modeltest=1
            [erho,evel] = solver_garz_minn_collapsed_ctm(tfinal,pos,B,V,ca_2p,...
            cb_2p,h,rho0,v0,f0,uf,rhof,q_min_arc,q_eq_arc,q_max_arc,w_max,rhom_cgarz,proj);
        
%             [erho,evel] = solver_garz_minn_collapsed(tfinal,pos,B,V,ca_2p,...
%             cb_2p,h,rho0,v0,f0,uf,rhof,q_min_arc,q_eq_arc,q_max_arc,w_max,rhom,proj)
            erho

            e_rho(1) = erho; % error from density
            e_vel(1) = evel; % error from velocity
%             % garz
%           modeltest=2
           [erho,evel] = solver_garz_minn_ctm(tfinal,pos,B,V,cl,...
             cp,ca,h,w_min,w_max,rhom);
            e_rho(2) = erho; % error from density
            e_vel(2) = evel; % error from velocity
%             
%             [erho,evel] = solver_garz_minn(tfinal,pos,B,V,cl,cp,ca,h,w_min,...
%                w_max,rhom);
%             e_rho(2) = erho; % error from density
%             e_vel(2) = evel; % error from velocity
            
            % arz 
%            [erho,evel] = solver_arz_minn(tfinal,pos,fd_equi,h,B,V,w_max,rhom);          
%             e_rho(3) = erho; % error from density
%             e_vel(3) = evel; % error from velocity
           [erho,evel] = solver_arz_minn_ctm(tfinal,pos,fd_equi,h,B,V,vm,w_max,rhom);          
            e_rho(3) = erho; % error from density
            e_vel(3) = evel; % error from velocity  

% %            [erho,evel] = solver_arzq_minn(tfinal,pos,vm,w_max,rhom,h,B,V);      
% %           [erho,evel] = solver_arzq_minn_ctm(tfinal,pos,94.1755,w_max,545.8877,h,B,V);         
%             [erho,evel] = solver_arzq_minn_ctm(tfinal,pos,115.7694,w_max,279.1618,h,B,V); 
%             e_rho(4) = erho; % error from density
%             e_vel(4) = evel; % error from velocity  
%             % lwr
           [erho,evel] = solver_lwr_minn_godunov(tfinal,pos,B,V,fd_equi,...
               h,w_max,rhom);
            e_rho(5) = erho; % error from density
            e_vel(5) = evel; % error from velocity%  
            % lwrq
%            [erho,evel] = solver_lwrq_minn_godunov(tfinal,pos,B,V,115.7694,w_max,h,279.1618);
%             e_rho(6) = erho; % error from density
%             e_vel(6) = evel; % error from velocity
%             % interp
             [erho,evel] = solver_interp_minn(tfinal,pos,B,V,h,w_max);  
             e_rho(7) = erho; % error from density
             e_vel(7) = evel; % error from velocity
%            phase-transition
           load GPT_q_0406;
           fd_pt1=func_rhoc(80.8879,para,94.1755,545.8811);



           u_r = 94.1755; rho_r = 2*rhof; 
% %        rho_r = rhof;
           qmax = u_r*fd_pt1.*(1-fd_pt1/rho_r);
           q_min = min(qmax)
           q_max = max(qmax)
           w_min=min(para);
           w_max=max(para);
%            [erho,evel] = solver_phase_transition_minn_ctm(tfinal,pos,B,V,u_r,...
%                rho_r,h,w_max,q_min,q_max,545.8811);
           [erho,evel] = solver_gpt_green_minn_ctm(tfinal,pos,B,V,...
    h,q_min,q_max,80.8879,94.1755,545.8811);



%            u_r = 95; rho_r = 533; 
% %            rho_r = rhof;
% 
% %            [erho,evel]=solver_phase_trans_minn(tfinal,pos,B,V,ct,...
% %                h,w_max,uf,rhom);
%            %---------------------------------------
% %            [erho,evel]=solver_phase_trans_minn_q(tfinal,pos,B,V,u_r,rho_r,...
% %                h,w_max,q_min,q_eq,q_max,rhom);
%            %---------------------------------------   
% %            [erho,evel] = solver_phase_transition_minn(tfinal,pos,B,V,u_r,...
% %                rho_r,h,w_max,q_min,q_max,rhom)
%            %---------------------------------------
%  %           [erho,evel] = solver_phase_transition_minn(tfinal,pos,B,V,93,...
%   %             rho_r,h,w_max,3967.0033,9981.3,300);
%                       [erho,evel] = solver_phase_transition_minn_ctm(tfinal,pos,B,V,u_r,...
%                 rho_r,h,w_max,967.0033,9981.3,533);
            e_rho(8) = erho % error from density
            e_vel(8) = evel % error from velocity 
            erho
            evel
           
%           q_eq = qmax(2);
%            [erho,evel]=solver_phase_trans_minn(tfinal,pos,B,V,ct,...
%                h,w_max,uf,rhom);
           %---------------------------------------
%            [erho,evel]=solver_phase_trans_minn_q(tfinal,pos,B,V,u_r,rho_r,...
%                h,w_max,q_min,q_eq,q_max,rhom);
           %---------------------------------------   
%            [erho,evel] = solver_phase_transition_minn(tfinal,pos,B,V,u_r,...
%                rho_r,h,w_max,q_min,q_max,rhom)
           %---------------------------------------   


 %           [erho,evel] = solver_gpt_green_minn_ctm(tfinal,pos,B,V,h,q_min,q_max,80.8879,94.1755,545.8811);
  

            e_rho(8) = erho; % error from density
            e_vel(8) = evel; % error from velocity 
            e_rho
            e_vel
%    end
      % same numerical result
%    rho_num(i,:) = rho;
%    vel_num(i,:) = vel;
%-------------------------------------- 
%    e_rho(i) = erho; % error from density
%    e_vel(i) = evel; % error from velocity
   %-------------------------------------
%    e_rho(i) = erho/delta_rho; % error from density
%    e_vel(i) = evel/delta_vel; % error from velocity
%    e(i) = e_rho(i)/delta_rho+e_vel(i)/delta_vel % the total error 
%    i
%    cost = cost+toc;
% 
% end

% error = mean(e);
% Report run time of solution
% fprintf('Time spent on computation: %0.3f seconds\n',cost)
%--------------------------------------
% % file names
% j = resol;

%    switch models
%        case 1    % garz model
           savefile1 = sprintf('data_error_minn_rho_%02d',day);   % error density
           savefile2 = sprintf('data_error_minn_vel_%02d',day);% vel error
%            savefile3 = sprintf('data_error_minn_garz_pattern_%01d',j);% total error
%            savefile4 = sprintf('data_num_minn_garz_rho_%02d',j);% total error
%            savefile5 = sprintf('data_num_minn_garz_vel_%02d',j);% total error
%            
%            
%        case 2    % arz model
%            savefile1 = sprintf('data_error_minn_arz_rho');   % error density
%            savefile2 = sprintf('data_error_minn_arz_vel');   % vel error
% %            savefile3 = sprintf('data_error_minn_arz_pattern_%01d',j);% total error
% %            savefile4 = sprintf('data_num_minn_arz_rho_%02d',j);% total error
% %            savefile5 = sprintf('data_num_minn_arz_vel_%02d',j);% total error           
%        case 3    % arzq model
%            savefile1 = sprintf('data_error_minn_arzq_rho');   % error density
%            savefile2 = sprintf('data_error_minn_arzq_vel');   % vel error
% %            savefile3 = sprintf('data_error_minn_arzq_pattern_%01d',j);% total error 
% %            savefile4 = sprintf('data_num_minn_arzq_rho_%02d',j);% total error
% %            savefile5 = sprintf('data_num_minn_arzq_vel_%02d',j);% total error           
%        case 4    % lwr model  
%            savefile1 = sprintf('data_error_minn_lwr_rho');   % error density
%            savefile2 = sprintf('data_error_minn_lwr_vel');   % vel error
% %            savefile3 = sprintf('data_error_minn_lwr_pattern_%01d',j);% total error
% %            savefile4 = sprintf('data_num_minn_lwr_rho_%02d',j);% total error
% %            savefile5 = sprintf('data_num_minn_lwr_vel_%02d',j);% total error           
%        case 5    % lwrq model
%            savefile1 = sprintf('data_error_minn_lwrq_rho');   % error density
%            savefile2 = sprintf('data_error_minn_lwrq_vel');   % vel error
% %            savefile3 = sprintf('data_error_minn_lwrq_pattern_%01d',j);% total error
% %            savefile4 = sprintf('data_num_minn_lwrq_rho_%02d',j);% total error
% %            savefile5 = sprintf('data_num_minn_lwrq_vel_%02d',j);% total error           
%        case 8    % test case, no traffic model
%            savefile1 = sprintf('data_error_minn_interp_rho');   % error density
%            savefile2 = sprintf('data_error_minn_interp_vel');   % vel error
% %            savefile3 =  sprintf('data_error_minn_interp_pattern_%01d',j);   % vel error
% %            savefile4 = sprintf('data_num_minn_interp_rho');% total error 
% %            savefile5 = sprintf('data_num_minn_interp_vel');% total error 
% 
%        case 9    % garz with arctangent flux
%            savefile1 = sprintf('data_error_minn_garz_rho_arctan');   % error density
%            savefile2 = sprintf('data_error_minn_garz_vel_arctan');   % vel error
% %            savefile3 = sprintf('data_error_minn_garz_arctan_%01d',j);% total error
%        case 7    % collapsed garz
%            savefile1 = sprintf('data_error_minn_cgarz_rho');   % error density
%            savefile2 = sprintf('data_error_minn_cgarz_vel');   % vel error
% %            savefile3 = sprintf('data_error_minn_cgarz_%01d_%01d',j,proj);% total error
%            
%        case 6    % collapsed garz
%            savefile1 = sprintf('data_error_minn_gpt_rho');   % error density
%            savefile2 = sprintf('data_error_minn_gpt_vel');   % vel error
% %            savefile3 = sprintf('data_error_minn_gpt_%01d',j);% total error
%  
%    end
   
%--------------------------------------   
% save error result
%  e_rho
%  e_vel
 
% save(savefile1,'e_rho');       % save numerical result
% save(savefile2,'e_vel');  
% % save(savefile3,'e');       % save numerical result
% % save(savefile4,'rho_num');  
% % save(savefile5,'vel_num');  
%----------------------------------------------

 cost = cost+toc;
fprintf('Time spent on computation: %0.3f seconds\n',cost)

%/////////////////
% sub functions
%==========================================================================
%==========================================================================
function y = u_max(lambda,p,alpha,rhom)
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = alpha.*(b-a+(lambda.^2.*p)./a)/rhom;
      %compute maximum velocity, first take derivative of eq3.3 in thesis,
      %than set rho=0
%----------------------------------------------      
% free flow branch
function y = fd_free(rho,uf,rhof)  % quadratic form
   y = uf*rho.*(1-rho/rhof);
%----------------------------------------------
function y = diff_fd_free(rho,uf,rhof)
   y = uf*(1-2*rho/rhof);      
   
   function y = func_rhoc(rhoc,q,vm,rhom)
       a = -rhoc.*vm./((2*rhoc-3*rhom).^2);
       ca = a;
       cb = -a.*(rhom+rhoc) + rhoc.*vm./(rhoc-rhom)-vm./(1+q);
       cc = a.*rhom.*rhoc -rhoc.*vm.*rhom./(rhoc-rhom);
       y = (-cb - sqrt(cb.^2-4*ca.*cc))./(2*ca);


     
