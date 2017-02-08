function plot_garz_freq_ngsim_100beta_temp
clear;
load GARZ_varyingIC-realistic_rhom0_fval_ngsim.mat
h = 1.0e-03;
beta = linspace(h,1-h,100);
m=length(beta);


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
figure(1)


rho_m=809.2603;


grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on

para0=zeros(4,1);
para0(4,1)=rho_m;
para0(1,1)=24.1319;
para0(2,1)=0.1611;
para0(3,1)=1450.9;



para=zeros(m,3);


lb1 = [0,0,0];
x01= [para0(1,1) para0(2,1) para0(3,1)];
eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);



for j = 1:m
    j
    
    
        %        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
        %        options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
        %options = optimset('MaxFunEvals', 1*10^10);
        [x,fval] = fmincon(@(x)beta(j)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),rho_m)...
            -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),rho_m)...
            -flux),0*flux))^2,x01,[],[],[],[],lb1,[],[]);
        para(j,1)=x(1);
        para(j,2)=x(2);
        para(j,3)=x(3);
        
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
    
        
    if j~=eq_id
        plot(grid,traffic_flux_smooth(grid,para(j,1),para(j,2),para(j,3),rho_m),'b-','linewidth',0.05),hold on
    else
        plot(grid,traffic_flux_smooth(grid,para(j,1),para(j,2),para(j,3),rho_m),'m-','linewidth',0.05),hold on
    end
        
    
end

savefile = sprintf('GARZ_para-0411-ngsim_updated');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'para');    % P is the parameters, n cross 3


savefile = sprintf('GARZ_rhom-0411-ngsim_updated');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'para0');    % P is the parameters, n cross 3









% grid = 0:rho_m;    % plot of the sequence of flow-density curves
% plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
% 
% for i =1:9
%     P_sort=sortrows(P,32+i);
%     plot(grid,traffic_flux_smooth(grid,P_sort(1,(i-1)*3+1),P_sort(1,(i-1)*3+2),P_sort(1,(i-1)*3+3),rho_m),'m-','linewidth',2.5),hold on
% end
% hold off
% 
% % for iter =1:10;
% %     scatter3(iter,iter,iter),hold on;
% % end
% % hold off
% 
% %         xlabel('initial condition index','fontsize',14);
% %         ylabel('\rho_{max}','fontsize',14);
% %         strt = ['Parameter-GARZ-rhom ',' \beta = ',num2str(beta(i)),' count=',num2str(count),'/',num2str(882)];
% %         strts = ['Parameter-GARZ-rhom-log',num2str(i),'.fig'];
% %         title(strt,'fontsize',14);
% 
% 
% 
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
