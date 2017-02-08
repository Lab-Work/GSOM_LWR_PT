clear;
load GARZ_varyingIC-realistic_rhom0_fval-0406-combine.mat
h = 1.0e-03;
%h = 1.0e-03;
beta = linspace(h,1-h,100);
m=length(beta);
eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);

rhom_garz=zeros(1,1);

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
size(density);
figure(1) 

P_sort=sortrows(P,m*3+5);
P_sort_rm=zeros(length(P(:,1)),2);
count=1;
for i =1:(length(P(:,1))-1)
    if P_sort(i,m*3+5)>P_sort(i+1,m*3+5)-0.5 && P_sort(i,m*3+5)<P_sort(i+1,m*3+5)+0.5
        count=count+1;
        P_sort_rm(i,1)=P_sort(i,m*3+5);
        P_sort_rm(i,2)=count;
    else
        P_sort_rm(i,1)=P_sort(i,m*3+5);
        P_sort_rm(i,2)=count;
        count=1;
    end
        
end
P_sort_rm_end=sortrows(P_sort_rm,2);
%rho_m=P_sort_rm_end(length(P(:,1)),1);
rho_m=491.5110;
num=P_sort_rm_end(length(P(:,1)),2);

index0 = find(P_sort_rm(:,2)==num)+1;
index=max(index0);
P_rm_0=P_sort( [index-num+1:index] , : );

para=zeros(m,3);
rhom_garz(1,1)=rho_m;


P_rm=P_rm_0;
for j=1:3
    P_sort=sortrows(P_rm,(eq_id-1)*3+j);
    P_sort_1=zeros(length(P_rm(:,1)),2);
    count=1;
    for i =1:(length(P_rm(:,1))-1)
        if P_sort(i,(eq_id-1)*3+j)>P_sort(i+1,(eq_id-1)*3+j)-0.1 && P_sort(i,(eq_id-1)*3+j)<P_sort(i+1,(eq_id-1)*3+j)+0.1
            count=count+1;
            P_sort_1(i,1)=P_sort(i,(eq_id-1)*3+j);
            P_sort_1(i,2)=count;
        else
            P_sort_1(i,1)=P_sort(i,(eq_id-1)*3+j);
            P_sort_1(i,2)=count;
            count=1;
        end
        
    end
    P_sort_1_end=sortrows(P_sort_1,2);
    para(eq_id,j)=P_sort_1_end(length(P_rm(:,1)),1);
    num=P_sort_1_end(length(P_rm(:,1)),2);
    
    index0 = find(P_sort_1(:,2)==num)+1;
    index=max(index0);
    P_rm=P_sort([index-num+1:index],:);
    
end

para(eq_id,1)=28.2751;
para(eq_id,2)=0.1669;
para(eq_id,3)=1033.6;

load data_dens_minn_1;
load data_vel_minn_1;
flux_1=Dens.*Vel;


    


grid = 0:rho_m;    % plot of the sequence of flow-density curves



lb1 = [0,0,0];
x01= [para(eq_id,1) para(eq_id,2) para(eq_id,3)];
eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);

load GARZ_para-0406.mat;


for j = 1:m
    j
    if j~=eq_id
    
        %        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
        %        options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
        %options = optimset('MaxFunEvals', 1*10^10);
%         [x,fval] = fmincon(@(x)beta(j)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),rho_m)...
%             -flux),0*flux))^2+(1-beta(j))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),rho_m)...
%             -flux),0*flux))^2,x01,[],[],[],[],lb1,[],[]);
%         para(j,1)=x(1);
%         para(j,2)=x(2);
%         para(j,3)=x(3);
        
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
    end
        
    if j~=eq_id
        plot(grid,traffic_flux_smooth(grid,para(j,1),para(j,2),para(j,3),rho_m),'m-','linewidth',0.01),hold on
    else
        plot(grid,traffic_flux_smooth(grid,para(j,1),para(j,2),para(j,3),rho_m),'m-','linewidth',0.01),hold on
    end
        
    
end
id_plot=26;
plot(Dens(id_plot,:)',flux_1(id_plot,:)','.','color',[0.2,0.2,0.2],'markersize',6), hold on

savefile = sprintf('GARZ_para-0406');
% savefile = sprintf('data_fd_minn_garz');
% 
%save(savefile,'para');    % P is the parameters, n cross 3


savefile = sprintf('GARZ_rhom');
% savefile = sprintf('data_fd_minn_garz');
% 
%save(savefile,'rhom_garz');    % P is the parameters, n cross 3
% for iter=1:9
%     P_rm=P_rm_0;
%     for j=1:3
%         P_sort=sortrows(P_rm,(iter-1)*3+j);
%         P_sort_1=zeros(length(P_rm(:,1)),2);
%         count=1;
%         for i =1:(length(P_rm(:,1))-1)
%             if P_sort(i,(iter-1)*3+j)>P_sort(i+1,(iter-1)*3+j)-0.1 && P_sort(i,(iter-1)*3+j)<P_sort(i+1,(iter-1)*3+j)+0.1
%                 count=count+1;
%                 P_sort_1(i,1)=P_sort(i,(iter-1)*3+j);
%                 P_sort_1(i,2)=count;
%             else
%                 P_sort_1(i,1)=P_sort(i,(iter-1)*3+j);
%                 P_sort_1(i,2)=count;
%                 count=1;
%             end
%             
%         end
%         P_sort_1_end=sortrows(P_sort_1,2);
%         para(iter,j)=P_sort_1_end(length(P_rm(:,1)),1);
%         num=P_sort_1_end(length(P_rm(:,1)),2);
%         
%         index0 = find(P_sort_1(:,2)==num)+1;
%         index=max(index0);
%         P_rm=P_sort([index-num+1:index],:);
%         
%     end
%     if iter~=7
%         plot(grid,traffic_flux_smooth(grid,para(iter,1),para(iter,2),para(iter,3),rho_m),'b-','linewidth',2.5),hold on
%     else
%         plot(grid,traffic_flux_smooth(grid,para(iter,1),para(iter,2),para(iter,3),rho_m),'m-','linewidth',2.5),hold on
%     end
%     
% end
% hold off
% % 
% % test=1;
% % for i=i:length(grid)-1
% %     for j=1:8
% %         if traffic_flux_smooth(i,para(j,1),para(j,2),para(j,3),rho_m)>=traffic_flux_smooth(i,para(j+1,1),para(j+1,2),para(j+1,3),rho_m)
% %             test=test*1;
% %         else
% %             i
% %             j
% %             test=test*0;
% %         end
% %     end
% % end
% % test
% % rho_m
% % para
% 
% 
% 
% 
% 
% 
% 
% 
% % grid = 0:rho_m;    % plot of the sequence of flow-density curves
% % plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
% % 
% % for i =1:9
% %     P_sort=sortrows(P,32+i);
% %     plot(grid,traffic_flux_smooth(grid,P_sort(1,(i-1)*3+1),P_sort(1,(i-1)*3+2),P_sort(1,(i-1)*3+3),rho_m),'m-','linewidth',2.5),hold on
% % end
% % hold off
% % 
% % % for iter =1:10;
% % %     scatter3(iter,iter,iter),hold on;
% % % end
% % % hold off
% % 
% % %         xlabel('initial condition index','fontsize',14);
% % %         ylabel('\rho_{max}','fontsize',14);
% % %         strt = ['Parameter-GARZ-rhom ',' \beta = ',num2str(beta(i)),' count=',num2str(count),'/',num2str(882)];
% % %         strts = ['Parameter-GARZ-rhom-log',num2str(i),'.fig'];
% % %         title(strt,'fontsize',14);
% % 
% % 
% % 
% % 
