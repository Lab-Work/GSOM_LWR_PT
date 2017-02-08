function error = comput_minn_validation(models,j,proj)
if nargin<2 models = 10; j = 1; proj = 1; end
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
rhom = 533;         % maximum density
%----------------------------------------------------
% Load data in the three test sensors
% density data
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
%-----------------------------------------------------
% load fd data from fitting process
% GARZ of the original 3 parameter model

% load data_fd_minn_533.mat   % beta = linspace(1e-04,1-1e-04,1000).
% load data_fd_minn_500.mat   % beta = linspace(1e-04,1-1e-04,1000).
% load data_fd_minn_arctan.mat
% fd_arctan = P;
load data_fd_minn_arctan_2p.mat
fd_collap = P;
load data_fd_minn_test.mat
% load data_fd_minn_533.mat
% load data_fd_rhom_minn_trb.mat
%**************************************************************
%**** calculate fitting pramaters for the Minneapolis data  ***
%**************************************************************
% calculate maximum velocity
% first locate the equilibrium curve.
midd = ceil(size(P,1)/2);         % the equilibium curve
equi_lda = P(midd,1);
equi_p = P(midd,2);
equi_alpha = P(midd,3);
fd_equi = [equi_lda,equi_p,equi_alpha];
% fd_equi = gen_data_fd_minn_trb(rhom);
% fd_equi = FD_EQ(8,:);
vm = u_max(fd_equi(1),fd_equi(2),fd_equi(3),rhom); % 93.3846 % maximum velocity
% rhoc = rho_critical(equi_lda,equi_p,rhom);
% q_max = traffic_flux_smooth(rhoc,equi_lda,equi_p,equi_alpha,rhom);
% vm_q = 4*q_max/rhom;
% vm_q = vm;
%-----------------------------------------------
% find out empty road velocity w
lambda = P(:,1);
p = P(:,2);
alpha = P(:,3);
w = func_W(lambda,p,alpha,rhom);  % values of w(\omega)
w_min = min(w);
w_max = max(w);

%\\\\\\\\\\\\\\\\\\\\\\\
% fitting fd parameters
%///////////////////////
order = 2;

% cl = polyfit(w(2:end),lambda(2:end),order);
% cp = polyfit(w(2:end),p(2:end),order);
% ca = polyfit(w(2:end),alpha(2:end),order);

% lambda0 = [lambda(1),lambda(midd),lambda(end)];
% p0 = [p(1),p(midd),p(end)];
% alpha0 = [alpha(1),alpha(midd),alpha(end)];
% w0 = [w(1),w(midd),w(end)];
 
% lambda0 = [lambda(2),lambda(midd),lambda(end)];
% p0 = [p(2),p(midd),p(end)];
% alpha0 = [alpha(2),alpha(midd),alpha(end)];
% w0 = [w(2),w(midd),w(end)];
%  
% % 
% cl = polyfit(w0,lambda0,order);
% cp = polyfit(w0,p0,order);
% ca = polyfit(w0,alpha0,order);

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
rho0 = 38.5; 
uf = 114; rhof = 1026;
f0 = fd_free(rho0,uf,rhof);
v0 = diff_fd_free(rho0,uf,rhof);
q = func_q(fd_collap(:,1),fd_collap(:,2),rho0,uf,rhof,rhom);

% ca_2p= polyfit(q,fd_collap(:,1),order);
% cb_2p = polyfit(q,fd_collap(:,2),order);

ca_2p= polyfit(q(1:end-1),fd_collap(1:end-1,1),order);
cb_2p = polyfit(q(1:end-1),fd_collap(1:end-1,2),order);

midd_arc = ceil(size(fd_collap,1)/2)+1;

q_eq_arc = q(midd_arc);   % equilibrium cuve
q_min_arc = min(q);
q_max_arc = max(q);
%---------------------------------------------
dist = [0.616 0.598 0.72 0.6975 0.6975];
pos = [0,dist(1),sum(dist(1:2))];
%-----------------------------------------------
vec_h = .032./2.^(0:6);   % in km.
h = vec_h(j);
%*********************************
%****  Numerical solvers  ********
%*********************************
e_rho = zeros(1,days);
e_vel = zeros(1,days);
e = zeros(1,days);
%------------------------------------------------
cost = 0;    % compute the time used for the program
% tfinal = 1/3;   % 20 minutes
tfinal = 1;
% rho_num
for i = 1:days
   tic
   B(1,:) = LBD(i,:);
   B(2,:) = D(i,:);
   B(3,:) = RBD(i,:);
   V(1,:) = LBV(i,:);
   V(2,:) = Vel(i,:);
   V(3,:) = RBV(i,:);
   %--------------------------------------------  
   switch models
       case 1  % cgarz
            [erho,evel] = solver_garz_minn_collapsed_ctm(tfinal,pos,B,V,ca_2p,...
            cb_2p,h,rho0,v0,f0,uf,rhof,q_min_arc,q_eq_arc,q_max_arc,w_max,rhom,proj);
        
%             [erho,evel] = solver_garz_minn_collapsed(tfinal,pos,B,V,ca_2p,...
%             cb_2p,h,rho0,v0,f0,uf,rhof,q_min_arc,q_eq_arc,q_max_arc,w_max,rhom,proj)
       case 2
%             % garz
           [erho,evel] = solver_garz_minn_ctm(tfinal,pos,B,V,cl,...
             cp,ca,h,w_min,w_max,rhom);
       case 3  % arz model
           [erho,evel] = solver_arz_minn_ctm(tfinal,pos,fd_equi,h,B,V,vm,w_max,rhom);          
       case 4  % arzq model
%            [erho,evel] = sol  ver_arzq_minn(tfinal,pos,vm,w_max,rhom,h,B,V);      
           [erho,evel] = solver_arzq_minn_ctm(tfinal,pos,vm,w_max,rhom,h,B,V);          

       case 5     % lwr
           [erho,evel] = solver_lwr_minn_godunov(tfinal,pos,B,V,fd_equi,...
               h,w_max,rhom);
       case 6
            % lwrq
           [erho,evel] = solver_lwrq_minn_godunov(tfinal,pos,B,V,vm,w_max,h,rhom);
       case 7
            % interp
            [erho,evel] = solver_interp_minn(tfinal,pos,B,V,h,w_max);  
       case 8 % phase-transition
           u_r = 95; rho_r = 1.0e+06*rhof; 
%            rho_r = rhof;
           qmax = uf*fd_pt.*(1-fd_pt/rho_r);
           q_min = min(qmax);
           q_max = max(qmax);
           q_eq = qmax(2);

           [erho,evel] = solver_phase_transition_minn_ctm(tfinal,pos,B,V,u_r,...
               rho_r,h,w_max,q_min,q_max,rhom);
       case 10
%            [erho,evel] = solver_arz_composit_minn(tfinal,pos,B,V,...
%              h,w_min,w_max,rhom);
%            vm = w_max;
%            vm = 95;
           vm = 106;
           [erho,evel] = f_arz(tfinal,pos,B,V,...
             h,vm,rhom);
%            [erho,evel] = solver_arz_composit_minn(tfinal,pos,B,V,...
%              h,vm,vm,rhom);
           
   end
      % same numerical result
%    rho_num(i,:) = rho;
%    vel_num(i,:) = vel;
%--------------------------------------

% save error data
%    
   e_rho(i) = erho; % error from density
   e_vel(i) = evel; % error from velocity
%--------------------------------------   
   e_rho(i) = erho/delta_rho; % error from density
   e_vel(i) = evel/delta_vel; % error from velocity
   e(i) = e_rho(i)+e_vel(i) % the total error 
%    i
   cost = cost+toc

end

error = mean(e);
% Report run time of solution
fprintf('Time spent on computation: %0.3f seconds\n',cost)
%--------------------------------------
% % file names
   switch models
       case 1    % garz model
           savefile1 = sprintf('data_error_minn_cgarz_rho');   % error density
           savefile2 = sprintf('data_error_minn_cgarz_vel');% vel error
           savefile3 = sprintf('data_error_minn_garz_%02d',j);% total error
           savefile4 = sprintf('data_num_minn_garz_rho_%02d',j);% total error
           savefile5 = sprintf('data_num_minn_garz_vel_%02d',j);% total error
           
           
       case 2    % arz model
           savefile1 = sprintf('data_error_minn_garz_rho');   % error density
           savefile2 = sprintf('data_error_minn_garz_vel');   % vel error
           savefile3 = sprintf('data_error_minn_arz_%02d',j);% total error
           savefile4 = sprintf('data_num_minn_arz_rho_%02d',j);% total error
           savefile5 = sprintf('data_num_minn_arz_vel_%02d',j);% total error           
       case 3    % arzq model
           savefile1 = sprintf('data_error_minn_arz_rho');   % error density
           savefile2 = sprintf('data_error_minn_arz_vel');   % vel error
           savefile3 = sprintf('data_error_minn_arzq_%01d',j);% total error 
           savefile4 = sprintf('data_num_minn_arzq_rho_%02d',j);% total error
           savefile5 = sprintf('data_num_minn_arzq_vel_%02d',j);% total error           
       case 4    % lwr model  
           savefile1 = sprintf('data_error_minn_arzq_rho');   % error density
           savefile2 = sprintf('data_error_minn_arzq_vel');   % vel error
           savefile3 = sprintf('data_error_minn_lwr_%02d',j);% total error
           savefile4 = sprintf('data_num_minn_lwr_rho_%02d',j);% total error
           savefile5 = sprintf('data_num_minn_lwr_vel_%02d',j);% total error           
       case 5    % lwrq model
           savefile1 = sprintf('data_error_minn_lwr_rho');   % error density
           savefile2 = sprintf('data_error_minn_lwr_vel');   % vel error
           savefile3 = sprintf('data_error_minn_lwrq_%02d',j);% total error
           savefile4 = sprintf('data_num_minn_lwrq_rho_%02d',j);% total error
           savefile5 = sprintf('data_num_minn_lwrq_vel_%02d',j);% total error           
       case 6    % test case, no traffic model
           savefile1 = sprintf('data_error_minn_lwrq_rho');   % error density
           savefile2 = sprintf('data_error_minn_lwrq_vel');   % vel error
           savefile3 =  sprintf('data_error_minn_interp_%02d',j);   % vel error
           savefile4 = sprintf('data_num_minn_interp_rho');% total error 
           savefile5 = sprintf('data_num_minn_interp_vel');% total error 

       case 7    % garz with arctangent flux
           savefile1 = sprintf('data_error_minn_interp_rho');   % error density
           savefile2 = sprintf('data_error_minn_interp_vel');   % vel error
           savefile3 = sprintf('data_error_minn_garz_arctan_%01d',j);% total error
       case 8    % collapsed garz
           savefile1 = sprintf('data_error_minn_gpt_rho');   % error density
           savefile2 = sprintf('data_error_minn_gpt_vel');   % vel error
           savefile3 = sprintf('data_error_minn_garz_collapsed_%02d',proj);% total error

 
   end
   
%--------------------------------------   
% save error result
 
% save(savefile1,'e_rho');       % save numerical result
% save(savefile2,'e_vel');  
% save(savefile3,'e');       % save numerical result
% save(savefile4,'rho_num');  
% save(savefile5,'vel_num');  
%----------------------------------------------

%/////////////////
% sub functions
%==========================================================================
%==========================================================================
function y = u_max(lambda,p,alpha,rhom)
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = alpha.*(b-a+(lambda.^2.*p)./a)/rhom;
%----------------------------------------------      
% free flow branch
function y = fd_free(rho,uf,rhof)  % quadratic form
   y = uf*rho.*(1-rho/rhof);
%----------------------------------------------
function y = diff_fd_free(rho,uf,rhof)
   y = uf*(1-2*rho/rhof);      
     
