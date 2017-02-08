function P = gen_data_fd_ngsim_gpt_newfd(rhom)

close all;
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
density = rho(indx)';                % density data/sensor
flux = vol(indx)';                % flux data/sensor
q_max = max(max(flux));           % ???
 

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


%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

count=0;

%YE GPT-----------------------------------------------------------------
initialur = [10,50,100];
initialrr = [1000,2000];
initialrm = [400,600,800];
initialq = [5000,7000,9000];
% initialur = [50];
% initialrr = [2000];
% initialrm = [400,600];
% initialq = [8000];

lb=[0,0,0,0];
A=[0 -1 1 0];
b=0;
A1=[-1];
b1=1
N1=length(initialur);
N2=length(initialrr);
N3=length(initialrm);
N4=length(initialq);
N=N1*N2*N3*N4;
m = length(beta);
P = zeros(N,9);  % initialize the matrix
% P = zeros(m,3);  % initialize the matrix
%YE GPT-----------------------------------------------------------------

for iter1 = 1:N1
    for iter2 = 1:N2
        for iter3 = 1:N3
            for iter4 =1:N4
                
                count=count+1;
                
                x0 = [initialur(iter1) initialrr(iter2) initialrm(iter3) initialq(iter4)];
              
                options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
                %options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
                %options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
                
                [x,fval] = fmincon(@(x)0.5*norm(max((func_flux_phase_transition(density,x(4),x(1),x(2),x(3))...
                    -flux),0*flux))^2+(0.5)*norm(max(-(func_flux_phase_transition(density,x(4),x(1),x(2),x(3))...
                    -flux),0*flux))^2,x0,A,b,[],[],lb,[],[],options);

                rho_m=x(3);
                rho_m
                
                P(count,1) = x(1);
                P(count,2) = x(2);
                P(count,3) = x(3);
                P(count,4) = x(4);
                P(count,5) = initialur(iter1);
                P(count,6) = initialrr(iter2);
                P(count,7) = initialrm(iter3);
                P(count,8) = initialq(iter4);
                P(count,9) = fval;
                
            end
        end
    end
    
    
    
end

P_sort=sort(P,9);

% fd_equi = P;
% rhoc = rho_critical(fd_equi(1),fd_equi(2),rhom)

% data_fd_ngsim_p3 have 3 parameters
% savefile = sprintf('data_fd_pt_minn_%1.0f',rhom);
%strsv = ['rhom0_',num2str(rhom0),'sigma0_',num2str(sigma0),'mu0_',num2str(mu0)]
savefile = sprintf('GPT_varyingIC_real_fval_ngsim_newptfd-fixrhotilde');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'P');    % P is the parameters, n cross 3

% load data_fd_minn_arctan.mat


%\\\\\\\\\\\\\\\\\\
%  The plots part
%//////////////////

figure(1)    % 

grid = 0:P_sort(1,3);    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
% plot(d_s1,f_s1,'r.'), hold on
% plot(d_s3,f_s3,'b*'), hold on

%---------------------------------------------------------
%--------------------------------------------------------
 plot(grid,func_flux_phase_transition(grid,P_sort(1,4),P_sort(1,1),P_sort(1,2),P_sort(1,3)),'m--','linewidth',4),hold on
%--------------------------------------------------------
hold off

%  3 parameters
%  plot(grid,traffic_flux_smooth(grid,P(1,1),P(1,2),P(1,3),rhom),'b-','linewidth',2.5),hold on
%  plot(grid,traffic_flux_smooth(grid,P(2,1),P(2,2),P(2,3),rhom),'m-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(3,1),P(3,2),P(3,3),rhom),'r-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(4,1),P(4,2),P(4,3),rhom),'k-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(5,1),P(5,2),P(5,3),rhom),'r-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(6,1),P(6,2),P(6,3),rhom),'m-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(7,1),P(7,2),P(7,3),rhom),'b-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(8,1),P(8,2),P(8,3),rhom),'k-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P(9,1),P(9,2),P(9,3),rhom),'b--','linewidth',2.5), hold on 
% % plot(grid,traffic_flux_smooth(grid,P(10,1),P(10,2),P(10,3),rhom),'k--','linewidth',2.5), hold on  
%  hold off
%--------------------------------------------------------
% % 2 parameters
% plot(grid,traffic_flux_smooth_2(grid,P(1,1),P(1,2),rho_m),'b-','linewidth',2.5),hold on
% plot(grid,traffic_flux_smooth_2(grid,P(2,1),P(2,2),rho_m),'r-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(3,1),P(3,2),rho_m),'k-','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(4,1),P(4,2),rho_m),'r-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(5,1),P(5,2),rho_m),'m-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(6,1),P(6,2),rho_m),'b-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(7,1),P(7,2),rho_m),'k-.','linewidth',2.5), hold on
% plot(grid,traffic_flux_smooth_2(grid,P(8,1),P(8,2),rho_m),'b--','linewidth',2.5), hold on 
% plot(grid,traffic_flux_smooth_2(grid,P(9,1),P(9,2),rho_m),'m-','linewidth',2.5), hold on
% % plot(grid,traffic_flux_smooth_2(grid,P(end,1),P(end,2),rhom),'k--','linewidth',2.5), hold on  
%hold off

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
% plot(grid,func_fd_seibold_2p(grid,P(1,1),P(1,2),rho0,v0,f0,uf,rhof,P(1,22)),'r--','linewidth',4),hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,3),P(2,4),rho0,v0,f0,uf,rhof,P(1,22)),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,5),P(1,6),rho0,v0,f0,uf,rhof,P(1,22)),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,7),P(1,8),rho0,v0,f0,uf,rhof,P(1,22)),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,9),P(1,10),rho0,v0,f0,uf,rhof,P(1,22)),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,11),P(1,12),rho0,v0,f0,uf,rhof,P(1,22)),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,13),P(1,14),rho0,v0,f0,uf,rhof,P(1,22)),'m-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,15),P(1,16),rho0,v0,f0,uf,rhof,P(1,22)),'k-','linewidth',4), hold on
% plot(grid,func_fd_seibold_2p(grid,P(1,17),P(1,18),rho0,v0,f0,uf,rhof,P(1,22)),'k-','linewidth',4), hold on
%plot(grid,func_fd_seibold_2p(grid,P(end,1),P(end,2),rho0,v0,f0,uf,rhof,rhom),'b-.','linewidth',4),

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

% legend('\omega = 0.01','\omega = 0.1','\omega = 0.5',...
%     '\omega = 0.9','\omega = 0.99')
%  legend('Data','\gamma = 0.001','\gamma = 0.01','\gamma = 0.1','\gamma = 0.5',...
%      '\gamma = 0.9','\gamma = 0.99','\gamma = 0.999')
%legend('Data','FD')
%GARZ-----------------------------------------------
% strl = ['\lambda_0 = ',num2str(lambda0)];
% strp = ['p_0 = ',num2str(p0)]
% stra = ['\alpha_0 = ',num2str(alpha0)]
% strr = ['\rho_{m0} = ',num2str(rhom0)]
% text(rhom-10,-200,'\rho_{max}','fontsize',14)
% text(rhom/3,2000,'\Omega_{f}','fontsize',18)
% text(rhom/2,8000,'\Omega_{c}','fontsize',18)
% text(rhom*0.6,8000,strl,'fontsize',14)
% text(rhom*0.6,9000,strp,'fontsize',14)
% text(rhom*0.6,10000,stra,'fontsize',14)
% text(rhom*0.6,7000,strr,'fontsize',14)
% text(rho_m-10,-200,num2str(rho_m),'fontsize',14)
%GARZ-----------------------------------------------

%CGARZ-----------------------------------------------
% strs = ['\sigma_0 = ',num2str(sigma0)];
% strm = ['\mu_0 = ',num2str(mu0)]
% strr = ['\rho_{m0} = ',num2str(rhom0)]
% text(rhom/3,2000,'\Omega_{f}','fontsize',18)
% text(rhom/2,8000,'\Omega_{c}','fontsize',18)
% text(rho_m*0.85,8000,strs,'fontsize',14)
% text(rho_m*0.85,9000,strm,'fontsize',14)
% text(rho_m*0.85,7000,strr,'fontsize',14)
%CGARZ-----------------------------------------------




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
    