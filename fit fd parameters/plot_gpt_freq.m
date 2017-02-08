clear;
load GPT_varyingIC_real_fval-NonUb.mat
beta = [.001,.01,.03,.05,.1,.3,.5,0.9,0.99];


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

para0=zeros(1,2);
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

grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on


para=zeros(9,1);

for iter=1:9
    if iter~=7
        P_rm=P_rm_0;
        for j=1:1
            P_sort=sortrows(P_rm,j);
            P_sort_1=zeros(length(P_rm(:,1)),2);
            count=1;
            for i =1:(length(P_rm(:,1))-1)
                if P_sort(i,iter)>P_sort(i+1,iter)-0.1 && P_sort(i,iter)<P_sort(i+1,iter)+0.1
                    count=count+1;
                    P_sort_1(i,1)=P_sort(i,iter);
                    P_sort_1(i,2)=count;
                else
                    P_sort_1(i,1)=P_sort(i,iter);
                    P_sort_1(i,2)=count;
                    count=1;
                end
                
            end
            P_sort_1_end=sortrows(P_sort_1,2);
            para(iter,1)=P_sort_1_end(length(P_rm(:,1)),1);
            num=P_sort_1_end(length(P_rm(:,1)),2);
            
            index0 = find(P_sort_1(:,2)==num)+1;
            index=max(index0);
            P_rm=P_sort([index-num+1:index],:);
            
        end
        plot(grid,func_flux_creen_pt(grid,para(iter,1),rho_c,v_m,rho_m),'b-','linewidth',2.5),hold on
    else
        plot(grid,func_flux_creen_pt(grid,0,rho_c,v_m,rho_m),'m-','linewidth',2.5),hold on
    end
    
end
hold off

test=1;
for i=i:length(grid)-1
    for j=1:8
        if func_flux_creen_pt(grid,para(j,1),rho_c,v_m,rho_m)>=func_flux_creen_pt(grid,para(j+1,1),rho_c,v_m,rho_m)
            test=test*1;
        else
            i
            j
            test=test*0;
        end
    end
end
test
rho_m