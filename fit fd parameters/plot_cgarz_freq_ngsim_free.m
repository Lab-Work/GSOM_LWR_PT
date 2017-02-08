function plot_cgarz_freq_ngsim_free
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

P_sort=sortrows(P,22);
P_sort_rm=zeros(length(P(:,1)),2);
count=1;
for i =1:(length(P(:,1))-1)
    if P_sort(i,22)>P_sort(i+1,22)-0.1 && P_sort(i,22)<P_sort(i+1,22)+0.1
        count=count+1;
        P_sort_rm(i,1)=P_sort(i,22);
        P_sort_rm(i,2)=count;
    else
        P_sort_rm(i,1)=P_sort(i,22);
        P_sort_rm(i,2)=count;
        count=1;
    end
        
end
P_sort_rm_end=sortrows(P_sort_rm,2);
rho_m=P_sort_rm_end(length(P(:,1)),1);
num=P_sort_rm_end(length(P(:,1)),2);

index0 = find(P_sort_rm(:,1)==rho_m)+1;
index=max(index0);
P_rm_0=P_sort( [index-num+1:index] , : );

para0=zeros(1,4);
for j=1:4
    P_sort=sortrows(P_rm_0,18+j);
    P_sort_1=zeros(length(P_rm_0(:,1)),2);
    count=1;
    for i =1:(length(P_rm_0(:,1))-1)
        if P_sort(i,18+j)>P_sort(i+1,18+j)-0.1 && P_sort(i,18+j)<P_sort(i+1,18+j)+0.1
            count=count+1;
            P_sort_1(i,1)=P_sort(i,18+j);
            P_sort_1(i,2)=count;
        else
            P_sort_1(i,1)=P_sort(i,18+j);
            P_sort_1(i,2)=count;
            count=1;
        end
        
    end
    P_sort_1_end=sortrows(P_sort_1,2);
    para0(1,j)=P_sort_1_end(length(P_rm_0(:,1)),1);
    num=P_sort_1_end(length(P_rm_0(:,1)),2);
    
    index0 = find(P_sort_1(:,1)==para0(1,j))+1;
    index=max(index0);
    P_rm_0=P_sort([index-num+1:index],:);
    
end

rho0=para0(1,1);
uf=para0(1,2);
rhof=para0(1,3);

grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on

para0


para=zeros(9,2);

for iter=1:9
    P_rm=P_rm_0;
    for j=1:2
        P_sort=sortrows(P_rm,(iter-1)*2+j);
        P_sort_1=zeros(length(P_rm(:,1)),2);
        count=1;
        for i =1:(length(P_sort(:,1))-1)
            if P_sort(i,(iter-1)*2+j)>P_sort(i+1,(iter-1)*2+j)-0.1 && P_sort(i,(iter-1)*2+j)<P_sort(i+1,(iter-1)*2+j)+0.1
                count=count+1;
                P_sort_1(i,1)=P_sort(i,(iter-1)*2+j);
                P_sort_1(i,2)=count;
            else
                P_sort_1(i,1)=P_sort(i,(iter-1)*2+j);
                P_sort_1(i,2)=count;
                count=1;
            end
            
        end
        P_sort_1_end=sortrows(P_sort_1,2);
        para(iter,j)=P_sort_1_end(length(P_rm(:,1)),1);
        num=P_sort_1_end(length(P_rm(:,1)),2);
        
        index0 = find(P_sort_1(:,1)==para(iter,j))+1;
        index=max(index0);
        P_rm=P_sort([index-num+1:index-1],:);
        
    end
    if iter==1 || iter==9
        plot(grid,func_fd_seibold_2p_updated(grid,para(iter,1),para(iter,2),rho0,uf,rhof,rho_m),'b-','linewidth',2.5),hold on
    elseif iter==7
            plot(grid,func_fd_seibold_2p_updated(grid,para(iter,1),para(iter,2),rho0,uf,rhof,rho_m),'m-','linewidth',2.5),hold on
    else
    end
    
end




load CGARZ_varyingIC-realistic_rhom0_fval-temp.mat
para=zeros(9,2);
size(P1)

for iter=1:9
    P_rm=P1;
    for j=1:2
        P_sort=sortrows(P_rm,(iter-1)*2+j);
        size(P_sort)
        P_sort_1=zeros(length(P_rm(:,1)),2);
        count=1;
        for i =1:(length(P_sort(:,1))-1)
            if P_sort(i,(iter-1)*2+j)>P_sort(i+1,(iter-1)*2+j)-0.1 && P_sort(i,(iter-1)*2+j)<P_sort(i+1,(iter-1)*2+j)+0.1
                count=count+1;
                P_sort_1(i,1)=P_sort(i,(iter-1)*2+j);
                P_sort_1(i,2)=count;
            else
                P_sort_1(i,1)=P_sort(i,(iter-1)*2+j);
                P_sort_1(i,2)=count;
                count=1;
            end
            
        end
        P_sort_1_end=sortrows(P_sort_1,2);
        para(iter,j)=P_sort_1_end(length(P_rm(:,1)),1);
        num=P_sort_1_end(length(P_rm(:,1)),2);
        
        index0 = find(P_sort_1(:,1)==para(iter,j))+1;
        index=max(index0);
        P_rm=P_sort([index-num+1:index-1],:);
        
    end
    if iter~=1 && iter~=9 && iter~=7
        plot(grid,func_fd_seibold_2p_updated(grid,para(iter,1),para(iter,2),rho0,uf,rhof,rho_m),'b-','linewidth',2.5),hold on
    end
    
end
hold off




% test=1;
% for i=i:length(grid)-1
%     for j=1:8
%         if func_fd_seibold_2p(i,para(j,1),para(j,2),rho0,v0,f0,uf,rhof,rho_m)>=func_fd_seibold_2p(i,para(j+1,1),para(j+1,2),rho0,v0,f0,uf,rhof,rho_m)
%             test=test*1;
%         else
%             i
%             j
%             test=test*0;
%         end
%     end
% end
% test
% rho_m
% para

% for iter =1:10;
%     scatter3(iter,iter,iter),hold on;
% end
% hold off

%         xlabel('initial condition index','fontsize',14);
%         ylabel('\rho_{max}','fontsize',14);
%         strt = ['Parameter-GARZ-rhom ',' \beta = ',num2str(beta(i)),' count=',num2str(count),'/',num2str(882)];
%         strts = ['Parameter-GARZ-rhom-log',num2str(i),'.fig'];
%         title(strt,'fontsize',14);


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



