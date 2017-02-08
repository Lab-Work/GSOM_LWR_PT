function plot_gpt_freq_ngsim_100beta_newptfd
clear;
load GPT_varyingIC_real_fval_ngsim.mat
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

P_sort=sortrows(P,12);
P_sort_rm=zeros(length(P(:,1)),2);
count=1;
for i =1:(length(P(:,1))-1)
    if P_sort(i,12)>P_sort(i+1,12)-0.1 && P_sort(i,12)<P_sort(i+1,12)+0.1
        count=count+1;
        P_sort_rm(i,1)=P_sort(i,12);
        P_sort_rm(i,2)=count;
    else
        P_sort_rm(i,1)=P_sort(i,12);
        P_sort_rm(i,2)=count;
        count=1;
    end
    
end
P_sort_rm_end=sortrows(P_sort_rm,2);
rho_m=P_sort_rm_end(length(P(:,1)),1);
num=P_sort_rm_end(length(P(:,1)),2);

index0 = find(P_sort_rm(:,2)==num)+1;
index=max(index0);
P_rm_0=P_sort( [index-num+1:index] , : );

para0=zeros(1,3);
for j=1:2
    P_sort=sortrows(P_rm_0,9+j);
    P_sort_1=zeros(length(P_rm_0(:,1)),2);
    count=1;
    for i =1:(length(P_rm_0(:,1))-1)
        if P_sort(i,9+j)>P_sort(i+1,9+j)-0.1 && P_sort(i,9+j)<P_sort(i+1,9+j)+0.1
            count=count+1;
            P_sort_1(i,1)=P_sort(i,9+j);
            P_sort_1(i,2)=count;
        else
            P_sort_1(i,1)=P_sort(i,9+j);
            P_sort_1(i,2)=count;
            count=1;
        end
        
    end
    P_sort_1_end=sortrows(P_sort_1,2);
    para0(1,j)=P_sort_1_end(length(P_rm_0(:,1)),1);
    num=P_sort_1_end(length(P_rm_0(:,1)),2);
    
    index0 = find(P_sort_1(:,2)==num)+1;
    index=max(index0);
    P_rm_0=P_sort([index-num+1:index],:);
    
end

rho_c=para0(1,1);
v_m=para0(1,2);
para0(1,3)=rho_m;

grid = 0:885.4343;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on


para=zeros(m,1);

x01 = [8855.3];
eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);

for iter=1:m
    iter
        
        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
        [x,fval] = fmincon(@(x)beta(iter)*norm(max((func_flux_phase_transition(density,x,80.7027,545.8773,885.4343)...
            -flux),0*flux))^2+(1-beta(iter))*norm(max(-(func_flux_phase_transition(density,x,80.7027,545.8773,885.4343)...
            -flux),0*flux))^2,x01,[],[],[],[],[],[],[],options);
        P_rm=P_rm_0;
        para(iter,1)=x;
%         for j=1:1
%             P_sort=sortrows(P_rm,j);
%             P_sort_1=zeros(length(P_rm(:,1)),2);
%             count=1;
%             for i =1:(length(P_rm(:,1))-1)
%                 if P_sort(i,iter)>P_sort(i+1,iter)-0.1 && P_sort(i,iter)<P_sort(i+1,iter)+0.1
%                     count=count+1;
%                     P_sort_1(i,1)=P_sort(i,iter);
%                     P_sort_1(i,2)=count;
%                 else
%                     P_sort_1(i,1)=P_sort(i,iter);
%                     P_sort_1(i,2)=count;
%                     count=1;
%                 end
%                 
%             end
%             P_sort_1_end=sortrows(P_sort_1,2);
%             para(iter,1)=P_sort_1_end(length(P_rm(:,1)),1);
%             num=P_sort_1_end(length(P_rm(:,1)),2);
%             
%             index0 = find(P_sort_1(:,2)==num)+1;
%             index=max(index0);
%             P_rm=P_sort([index-num+1:index],:);
%             
%         end
if iter==eq_id
    plot(func_flux_phase_transition(grid,para(iter,1),80.7027,545.8773,885.4343),'m-','linewidth',0.05),hold on
else
    plot(func_flux_phase_transition(grid,para(iter,1),80.7027,545.8773,885.4343),'b-','linewidth',0.05),hold on
end
        
 
    
end
hold off

savefile = sprintf('GPT_q_0411-ngsim_newptfd');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'para');    % P is the parameters, n cross 3

savefile = sprintf('GPT_q_0411-ngsim-paraeq');
% savefile = sprintf('data_fd_minn_garz');
% 
%save(savefile,'para0'); 



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