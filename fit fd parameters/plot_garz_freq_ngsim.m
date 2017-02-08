function plot_garz_freq_ngsim
clear;
load GARZ_varyingIC-realistic_rhom0_fval_ngsim.mat
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
q_max = max(max(flux));           % ???
figure(1)

P_sort=sortrows(P,32);
P_sort_rm=zeros(length(P),2);
count=1;
for i =1:(length(P(:,1))-1)
    if P_sort(i,32)>P_sort(i+1,32)-0.5 && P_sort(i,32)<P_sort(i+1,32)+0.5
        count=count+1;
        P_sort_rm(i,1)=P_sort(i,32);
        P_sort_rm(i,2)=count;
    else
        P_sort_rm(i,1)=P_sort(i,32);
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

grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on

para=zeros(9,3);

for iter=1:9
    P_rm=P_rm_0;
    for j=1:3
        P_sort=sortrows(P_rm,(iter-1)*3+j);
        P_sort_1=zeros(length(P_rm(:,1)),2);
        count=1;
        for i =1:(length(P_rm(:,1))-1)
            if P_sort(i,(iter-1)*3+j)>P_sort(i+1,(iter-1)*3+j)-0.1 && P_sort(i,(iter-1)*3+j)<P_sort(i+1,(iter-1)*3+j)+0.1
                count=count+1;
                P_sort_1(i,1)=P_sort(i,(iter-1)*3+j);
                P_sort_1(i,2)=count;
            else
                P_sort_1(i,1)=P_sort(i,(iter-1)*3+j);
                P_sort_1(i,2)=count;
                count=1;
            end
            
        end
        P_sort_1_end=sortrows(P_sort_1,2);
        para(iter,j)=P_sort_1_end(length(P_rm(:,1)),1);
        num=P_sort_1_end(length(P_rm(:,1)),2);
        
        index0 = find(P_sort_1(:,2)==num)+1;
        index=max(index0);
        P_rm=P_sort([index-num+1:index],:);
        
    end
    if iter~=7
        plot(grid,traffic_flux_smooth(grid,para(iter,1),para(iter,2),para(iter,3),rho_m),'b-','linewidth',2.5),hold on
    else
        plot(grid,traffic_flux_smooth(grid,para(iter,1),para(iter,2),para(iter,3),rho_m),'m-','linewidth',2.5),hold on
    end
    
end
hold off


test=1;
for i=i:length(grid)-1
    for j=1:8
        if traffic_flux_smooth(i,para(j,1),para(j,2),para(j,3),rho_m)>=traffic_flux_smooth(i,para(j+1,1),para(j+1,2),para(j+1,3),rho_m)
            test=test*1
        else
            i
            j
            test=test*0
        end
    end
end
test
rho_m
para







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
