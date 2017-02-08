function [er_mean, ev_mean] = comput_ngsim_validation_cgarz_paper(dataset,models,k)
if nargin<1 dataset = 1; models = 1;k = 3; end
%========================================================================
% This code perform validation with NGSIM traffic data. Then perform
% comparison of GARZ model with classical model ARZ, ARZQ, LWR, and LWRQ.
% This validation processes are respect to 3 NGSIM data set, that is
% 4:4:15, 5:00-5:1 5, and 5:15--5:30
%**************************************
% switch dataset
%    case 1 -- 4:00-4:15
%    case 2 -- 5:00-5:15
%    case 3 -- 5:15-5:30
% end
%**************************************
% models 
%    case 1 -- garz
%    case 2 -- arz
%    case 3 -- arzq
%    case 4 -- lwr
%    case 5 -- lwrq
%    case 6 -- test case, linear interpolation
%    k --- refine the grids
%--------------------------------------
% Shimao Fan, Feb 10, 2012
% Math Department, Temple University.
% Last modified on Jan 24 2013
%========================================================================
% clc;
close all;
%-----------------------------------------
% load traffic data 
switch dataset
    case 1
        tfinal = 13.5;

%       load data_rho_ngsim_smooth_1_3.mat     % D / dt = 12 seconds
%       load data_vel_ngsim_smooth_1_3.mat    % V 
%       
      load data_density_ngsim_smooth_1_3.mat     % D / 15 m
      load data_velocity_ngsim_smooth_1_3.mat    % V 

%       load data_density_ngsim_ref_1.mat     % D / 15 m
%       load data_velocity_ngsim_ref_1.mat    % V 

    case 2
        tfinal = 13;

%       load data_rho_ngsim_smooth_2_3.mat     % D / dt = 12 seconds
%       load data_vel_ngsim_smooth_2_3.mat    % V 
%       
      load data_density_ngsim_smooth_2_3.mat     % D / 15 m
      load data_velocity_ngsim_smooth_2_3.mat    % V 

%       load data_density_ngsim_ref_2.mat     % D / 15 m
%       load data_velocity_ngsim_ref_2.mat    % V 
    
    case 3   
        tfinal = 12.5;

%       load data_rho_ngsim_smooth_3_3.mat     % D / dt = 12 seconds
%       load data_vel_ngsim_smooth_3_3.mat    % V 
%       
      load data_density_ngsim_smooth_3_3.mat     % D / 15 m
      load data_velocity_ngsim_smooth_3_3.mat    % V 

%       load data_density_ngsim_ref_3.mat     % D / 15 m
%       load data_velocity_ngsim_ref_3.mat    % V
     
end

D0 = D;
V0 = V;
%--------------------------------------------
% calculate the rescale parameter in the error measurement
% absolute average
rho_avg = 459.6;
vel_avg = 66.5;
%------------------------------------------------
load GPT_q_0411-ngsim.mat
fd_pt = para;
load data_fd_ngsim_arctan_3p.mat
fd_arctan = P;

% load data_fd_ngsim_arctan_2p.mat
load CGARZ_100beta_2para-ngsim.mat;
%load CGARZ_100beta_2para-ngsim_updated-0516.mat
fd_collap = P2;

% load data_fd_ngsim_rhom_800.mat      % matrix P
load GARZ_para-0411-ngsim_updated.mat      % matrix P
load GARZ_rhom-0411-ngsim_updated.mat

rhom = para0(4);
%-----------------------------------------
sg = .00025;               % finiest grids size/25 cm 
refine = 2.^(0:5);         % refine the girds, if k = 6
r0 = refine(end-1);          % 16, 16*sg = .004, 4 meters
ref = refine(k);           % grid size  ref = 2^5
g_size = refine(7-k);      % pick all data refine(1), 7-k follows
a = .022; b = .5;        % study area in km
xi = a:sg:b;      % truncate the study area into a smaller region
% trancate the study area tranc = 1
init = 14/.25;   % cut 14 meters both left and right
x = xi(init+1:g_size:end-init);    % spacial grids  xi = .036:sg:.486;

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
N0 = 57;   % need to change this number
%N0 = length(x)
%----------------------------------------------
BD(:,1) = D(:,init);         % begining 49 th feet position
BD(:,2) = D(:,end-init+1);      % end of study area(1651) position
BV(:,1) = V(:,init);
BV(:,2) = V(:,end-init+1);

% the simulation traffic data [50 1650]
D = D(:,init+1:g_size:end-init);  % data for calculation of numerical results
V = V(:,init+1:g_size:end-init);
%----------------------------------------------
% w vectors
lambda = para(:,1);
p = para(:,2);
alpha = para(:,3);
w = func_W(lambda,p,alpha,rhom);
%%%%%%%%%%%%%%%%
w_max = max(w);
w_min = min(w);
% vm = w_max;
%%%%%%%%%%%%%%%%
% calculate the maximum velocity
n = size(para,1);
midd = ceil(n/2);   % pick the middle points
lambda_equi = lambda(midd);
p_equi = p(midd);
alpha_equi = alpha(midd);
fd_equi = [lambda_equi,p_equi,alpha_equi];
vm = u_max(lambda_equi,p_equi,alpha_equi,rhom); % the maximum velocity(equi)

% calculate the critical density

% vm_q = 4*q_max/rhom;
% % vm_q = 59;
% vm_q = vm;
% %------------------------------------------
% polynomial fitting process
order = 10;
% 
% cl= polyfit(w,lambda,order);
% cp = polyfit(w,p,order);
% ca = polyfit(w,alpha,order);



% 
cl= polyfit(w(1:end),lambda(1:end),order)
cp = polyfit(w(1:end),p(1:end),order)
ca = polyfit(w(1:end),alpha(1:end),order)

grid = (min(w(1:end))):(max(w(1:end-1))); 

subplot(2,3,1);
plot(w(1:end),lambda(1:end),'.','color',[0.2,0.2,0.2],'markersize',6), hold on
plot(grid, polyval(cl,grid) ,'b-','linewidth',1),hold off

subplot(2,3,2);
plot(w(1:end),p(1:end),'.','color',[0.2,0.2,0.2],'markersize',6), hold on
plot(grid, polyval(cp,grid) ,'b-','linewidth',1),hold off

subplot(2,3,3);
plot(w(1:end),alpha(1:end),'.','color',[0.2,0.2,0.2],'markersize',6), hold on
plot(grid, polyval(ca,grid) ,'b-','linewidth',1),hold off


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% replace w by q (maximum flux)
%////////////////////////////////
rho_c = rho_critical(lambda,p,rhom);   % critical density
q = traffic_flux_smooth(rho_c,lambda,p,alpha,rhom); % maximum flux

clq= polyfit(q(2:end),lambda(2:end),order);
cpq = polyfit(q(2:end),p(2:end),order);
caq = polyfit(q(2:end),alpha(2:end),order);


%------------------------------
% fitting with respect 3-parameter arctangent flux
%w_arctan = func_w_arctan(fd_arctan(:,1),fd_arctan(:,2),fd_arctan(:,3),rhom);

% c_sigma= polyfit(w_arctan,fd_arctan(:,1),order);
% c_mu = polyfit(w_arctan,fd_arctan(:,2),order);
% c_alpha = polyfit(w_arctan,fd_arctan(:,3),order);

% c_sigma= polyfit(w_arctan(2:end),fd_arctan(2:end,1),order);
% c_mu = polyfit(w_arctan(2:end),fd_arctan(2:end,2),order);
% c_alpha = polyfit(w_arctan(2:end),fd_arctan(2:end,3),order);
% 
% w_min_arc = min(w_arctan);
% w_max_arc = max(w_arctan);

%-----------------------------------------
% fitting with respect 2-parameter arctangent flux
% rho0 = 67; 
% uf = 79; rhof = 1266;
% uf = 71; rhof = 1266;
% uf = 68;   rhof = 1057;
% uf = 73.6; rhof = 1266;
%-----------------------
% new parameters

load CGARZ_eq_6para-ngsim_updated.mat;
%load CGARZ_eq_6para-ngsim_updated-0516.mat

uf = para_eq(4,1);  rhof = para_eq(5,1);
rho0 = para_eq(3,1);
f0 = fd_free(rho0,uf,rhof);
v0 = diff_fd_free(rho0,uf,rhof);
rhom_cgarz=para_eq(6,1);
q = func_q(fd_collap(:,1),fd_collap(:,2),rho0,uf,rhof,rhom_cgarz);

order1=10;

ca_2p= polyfit(q,fd_collap(:,1),order1);
cb_2p = polyfit(q,fd_collap(:,2),order1);

grid1 = (min(q)):(max(q)); 

subplot(2,3,4);
plot(q,fd_collap(:,1),'.','color',[0.2,0.2,0.2],'markersize',6), hold on
plot(grid1, polyval(ca_2p,grid1) ,'b-','linewidth',1),hold off

subplot(2,3,5);
plot(q,fd_collap(:,2),'.','color',[0.2,0.2,0.2],'markersize',6), hold on
plot(grid1, polyval(cb_2p,grid1) ,'b-','linewidth',1),hold off

% ca_2p= polyfit(q(1:end-1),fd_collap(1:end-1,1),order);
% cb_2p = polyfit(q(1:end-1),fd_collap(1:end-1,2),order);

% ca_2p= polyfit(q(2:end),fd_collap(2:end,1),order);
% cb_2p = polyfit(q(2:end),fd_collap(2:end,2),order);


midd_arc = ceil(size(fd_collap,1)/2)+1;

q_eq_arc = q(midd_arc);   % equilibrium cuve
q_min_arc = min(q);
q_max_arc = max(q);

vm_q=79.95;
rhom_q=442.6201;
%-----------------------------------------
%*********************************
%****  Numerical solvers  ********
%*********************************

%\\\\\\\\\\\\\\\\\\\\
% main algorithm \\\\
%////////////////////

    switch models
        case 1 % garz
%            [num_rho,num_vel,erho,evel] = solver_garz_ngsim(tfinal,D,V,BD,BV,...
%                 cl,cp,ca,x,ref,w_min,w_max,rhom);
% num_rho            
           [num_rho,num_vel,erho,evel] = solver_garz_ngsim_ctm(tfinal,D,V,BD,BV,...
               cl,cp,ca,x,ref,w_min,w_max,rhom);
           
       

%            [erho,evel] = solver_garz_ngsim(D,V,BD,BV,...
%                 clq,cpq,caq,x,ref,q_min,q_max,rhom);
            
           savefile = sprintf('data_error_ngsim_garz_avg_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_garz_rho_%01d',dataset);   % vel error
           savefile2 = sprintf('data_error_ngsim_garz_vel_%01d',dataset); 
%            savefile4 = sprintf('data_error_ngsim_garz_avg3_%01d',dataset);   % vel error
%            savefile5 = sprintf('data_error_ngsim_garz_avg4_%01d',dataset);   % vel error          
        case 2 % arz
           [num_rho,num_vel,erho,evel] = solver_arz_ngsim_ctm(tfinal,fd_equi,D,...
                          V,BD,BV,x,ref,vm,rhom);
           savefile = sprintf('data_error_ngsim_arz_avg_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_arz_rho_%01d',dataset);   % vel error
           savefile2 = sprintf('data_error_ngsim_arz_vel_%01d',dataset); 
%            savefile4 = sprintf('data_error_ngsim_arz_avg3_%01d',dataset); 
%            savefile5 = sprintf('data_error_ngsim_arz_avg4_%01d',dataset); 
          
        case 3 % arzq
         
           [num_rho,num_vel,erho,evel] = solver_arzq_ngsim_ctm(tfinal,D,V,BD,BV,x,ref,vm_q,rhom_q);
           savefile = sprintf('data_error_ngsim_arzq_avg_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_arzq_rho_%01d',dataset);   % vel error
           savefile2 = sprintf('data_error_ngsim_arzq_vel_%01d',dataset); 
%            savefile4 = sprintf('data_error_ngsim_arzq_avg3_%01d',dataset); 
%            savefile5 = sprintf('data_error_ngsim_arzq_avg4_%01d',dataset); 
           
        case 4 % lwr
%            [erho,evel,erhos,evels] = solver_lwr_ngsim(fd_equi,D,V,BD,x,vm,rhom,ref);
%            [num_rho,num_vel,erho,evel] = solver_lwr_ngsim_godunov(fd_equi,D,V,BD,x,vm,rhom,ref);
           [num_rho,num_vel,erho,evel] = solver_lwr_ngsim_godunov(tfinal,fd_equi,D,V,BD,x,vm,rhom,ref);

           savefile = sprintf('data_error_ngsim_lwr_avg_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_lwr_rho_%01d',dataset);   % vel error
           savefile2 = sprintf('data_error_ngsim_lwr_vel_%01d',dataset);       % total error
%            savefile4 = sprintf('data_error_ngsim_lwr_avg3_%01d',dataset);       % total error
%            savefile5 = sprintf('data_error_ngsim_lwr_avg4_%01d',dataset);       % total error
        
        case 5 % lwrq
            
%            [erho,evel,erhos,evels] = solver_lwrq_ngsim(D,V,BD,x,vm,vm_q,rhom,ref);
           [num_rho,num_vel,erho,evel] = solver_lwrq_ngsim_godunov(tfinal,D,V,BD,x,vm_q,rhom_q,ref);
           savefile = sprintf('data_error_ngsim_lwrq_avg_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_lwrq_rho_%01d',dataset);   % vel error
           savefile2 = sprintf('data_error_ngsim_lwrq_vel_%01d',dataset);       % total error
%            savefile2 = sprintf('data_error_ngsim_lwrq_avg1_%01d',dataset);   % vel error
%            savefile3 = sprintf('data_error_ngsim_lwrq_avg2_%01d',dataset);       % total error  
%            savefile4 = sprintf('data_error_ngsim_lwrq_avg3_%01d',dataset);       % total error  
%            savefile5 = sprintf('data_error_ngsim_lwrq_avg4_%01d',dataset);       % total error             
        case 7  % phase transition model
            load GPT_q_0516-ngsim_newptfd.mat;
            load GPT_q_0516-ngsim-paraeq.mat;
            
            vm =para0(1,1); rho_tilde=para0(2,1); rhom_gpt=para0(3,1);
            %            q_min = -0.7096;
            %            q_max = 0.5602;
            q_min = -.4072; q_max = .4349;
            
            [erho,evel] = solver_phase_trans_ngsim_ctm(tfinal,D,V,BD,BV,...
                x,ref,vm,rho_tilde,min(para),max(para),rhom_gpt);
            
%            vm = 80.7027; rhoc = 132.0909; rhom_gpt=885.4343;
% %            q_min = -0.7096;
% %            q_max = 0.5602;
%            q_min = -.4072; q_max = .4349;
%            
%            [erho,evel] = solver_phase_trans_ngsim_ctm(tfinal,D,V,BD,BV,...
%     x,ref,vm,545.8773,5467.2,10865,rhom_gpt);
           

%            [num_rho,num_vel,erho,evel] = solver_gpt_green_minn_ctm(tfinal,D,V,BD,BV,...
%                  x,ref,q_min,q_max,rhoc,vm,rhom_gpt);
             
%            savefile = sprintf('data_error_ngsim_phase_trans_avg_%01d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_gpt_rho_%01d',dataset);   % vel error
           savefile2 = sprintf('data_error_ngsim_gpt_vel_%01d',dataset); 
%            savefile4 = sprintf('data_error_ngsim_phase_trans_avg3_%01d',dataset); 
%            savefile5 = sprintf('data_error_ngsim_phase_trans_avg4_%01d',dataset);            
           
%         case 8  % phase transition model
% %            [erho,evel] = solver_phase_trans_ngsim(D,V,BD,BV,ct,...
% %                     x,ref,sl,vm,rhom);
% %            u_r = 73; rho_r = 1266;
% %            u_r = 79; rho_r = 700;
%            u_r = vm; rho_r = 1.0e+06*1266;
%            
% %            qmax = uf*fd_pt.*(1-fd_pt/rho_r);
%            q_min = min(fd_pt);
%            q_max = max(fd_pt);%            
% %            q_min = 4000; q_max = 12000;
%            [num_rho,num_vel,erho,evel] = solver_phase_trans_ngsim_ctm(tfinal,D,V,BD,BV,...
%             x,ref,u_r,rho_r,q_min,q_max,vm,rhom);
% %            savefile = sprintf('data_error_ngsim_phase_trans_avg_%01d',dataset);   % error density
%            savefile1 = sprintf('data_error_ngsim_gpt_rho_%01d',dataset);   % vel error
%            savefile2 = sprintf('data_error_ngsim_gpt_vel_%01d',dataset); 
% %            savefile4 = sprintf('data_error_ngsim_phase_trans_avg3_%01d',dataset); 
% %            savefile5 = sprintf('data_error_ngsim_phase_trans_avg4_%01d',dataset); 
           
        case 9 % garz with arctangent flux
           [erho,evel] = solver_garz_ngsim_arctan(D,V,BD,BV,c_sigma,c_mu,c_alpha,...
              x,ref,w_min_arc,w_max_arc,rhom);
%             [erho,evel] = solver_garz_ngsim_2(D,V,BD,BV,...
%                 cl,cp,x,ref,w_min,w_max,vm,rhom);
           savefile = sprintf('data_error_ngsim_garz_arctan_avg_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_garz_arctan_rho_%01d',dataset);   % error density
           savefile2 = sprintf('data_error_ngsim_garz_arctan_vel_%01d',dataset);   % error density

%            savefile2 = sprintf('data_error_ngsim_garz_avg1_%01d',dataset);   % vel error
%            savefile3 = sprintf('data_error_ngsim_garz_avg2_%01d',dataset); 
        case 6 % garz with arctangent flux
           [num_rho,num_vel,erho,evel] = solver_garz_ngsim_collapsed_ctm(tfinal,D,V,BD,BV,ca_2p,cb_2p,...
            rho0,v0,f0,uf,rhof,x,ref,q_min_arc,q_eq_arc,q_max_arc,rhom_cgarz);
        
           savefile = sprintf('data_error_ngsim_cgarz_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_cgarz_rho_%01d',dataset);   % error density
           savefile2 = sprintf('data_error_ngsim_cgarz_vel_%01d',dataset);   % error density

        case 8 % linear interpolation
%            [erho,evel] = solver_interp_ngsim(tfinal,D,V,BD,BV,x,a,b,vm,rhom,ref);
           [num_rho,num_vel,erho,evel] = solver_interp_ngsim_test(tfinal,D,V,D0,V0,x,xi,vm,rhom,ref);
           savefile = sprintf('data_error_ngsim_interp_avg_%03d',dataset);   % error density
           savefile1 = sprintf('data_error_ngsim_interp_rho_%01d',dataset);   % vel error
           savefile2 = sprintf('data_error_ngsim_interp_vel_%01d',dataset); 
%            savefile4 = sprintf('data_error_ngsim_interp_avg3_%01d',dataset); 
%            savefile5 = sprintf('data_error_ngsim_interp_avg4_%01d',dataset); 
        case 10  % test case , w = rho_1/rho, the fraction of cars
          [erho,evel] = solver_arz_composit_ngsim(tfinal,D,V,BD,BV,...
             x,ref,w_min,w_max,rhom);
    end
    

    
    erho = erho/(5*N0);
    evel = evel/N0;
    er_mean=mean(erho)
    ev_mean=mean(evel)
    
%     e_avg0 = erho/rhom+evel/vm;
%     e_mean(1) = mean(e_avg0); % the original error
%     
%     erhos = erhos/N0;
%     evels = evels/N0;
%     e_avg1 = evels+erhos;
%     e_mean(2) = mean(e_avg1); % rescale with data
% %------------------------------------------
%     e_avg2 = erho/rho_avg1+evel/vel_avg1;
%     e_mean(3) = mean(e_avg2);   % absolute average
% %------------------------------------------
%     e_avg3 = erho/rho_avg2+evel/vel_avg2;
%     e_mean(4) = mean(e_avg3);
    % rescale with likelyhood maximum density and velocity
%     e_avg = erho/rho_avg+evel/vel_avg;
%     e_mean = mean(e_avg)
    
%\\\\\\\\\\\\\\\\\\\\
% Save results
% %////////////////////

% save(savefile,'e_avg');
save(savefile1,'erho')
save(savefile2,'evel')
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% compute maximum velocity
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

function y = u_max(lambda,p,alpha,rhom)
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = alpha.*(b-a+(lambda.^2.*p)./a)/rhom;
%----------------------------------------------      
% determine the critical density, for any given smooth curve
function y = rho_critical(lambda,p,rhom)
      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = rhom*((b-a)./(lambda.*(sqrt(lambda.^2-(b-a).^2)))+p); 
%----------------------------------------------      
% free flow branch
function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);
%----------------------------------------------
function y = diff_fd_free(rho,uf,rhof)

y = uf*(1-2*rho/rhof);      

