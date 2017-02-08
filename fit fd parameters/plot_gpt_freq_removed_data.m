clear;
load GPT_varyingIC_real_fval-NonUb-0406.mat

h = 1.0e-03;
beta = linspace(h,1-h,100);
m=length(beta);
% h = 1.0e-03;
% h = .01;
% beta = h/2:h:1-h/2;
% beta = .5;
%m = 100;
% beta = linspace(h,1-h,2);
% beta = [h,.5,1-h];
% % test beta data
%beta = [.001,.01,.03,.05,.1,.3,.5,0.9,0.99];

eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);


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
size(density)

figure(1)

P_sort=sortrows(P,3);
P_sort_rm=zeros(length(P(:,1)),2);
count=1;
for i =1:(length(P(:,1))-1)
    if P_sort(i,3)>P_sort(i+1,3)-0.1 && P_sort(i,3)<P_sort(i+1,3)+0.1
        count=count+1;
        P_sort_rm(i,1)=P_sort(i,3);
        P_sort_rm(i,2)=count;
    else
        P_sort_rm(i,1)=P_sort(i,3);
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
    P_sort=sortrows(P_rm_0,j);
    P_sort_1=zeros(length(P_rm_0(:,1)),2);
    count=1;
    for i =1:(length(P_rm_0(:,1))-1)
        if P_sort(i,j)>P_sort(i+1,j)-0.1 && P_sort(i,j)<P_sort(i+1,j)+0.1
            count=count+1;
            P_sort_1(i,1)=P_sort(i,j);
            P_sort_1(i,2)=count;
        else
            P_sort_1(i,1)=P_sort(i,j);
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


para=zeros(m,1);

x01 = [0];

for iter=1:m
    iter
    if iter~=eq_id
        
        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
        [x,fval] = fmincon(@(x)beta(iter)*norm(max((func_flux_creen_pt(density,x,rho_c,v_m,rho_m)...
            -flux),0*flux))^2+(1-beta(iter))*norm(max(-(func_flux_creen_pt(density,x,rho_c,v_m,rho_m)...
            -flux),0*flux))^2,x01,[],[],[],[],-1,[],[],options);
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
        plot(grid,func_flux_creen_pt(grid,para(iter,1),rho_c,v_m,rho_m),'b-','linewidth',0.05),hold on
    else
        plot(grid,func_flux_creen_pt(grid,0,rho_c,v_m,rho_m),'m-','linewidth',0.05),hold on
    end
    
end
hold off

savefile = sprintf('GPT_q_0406');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'para');    % P is the parameters, n cross 3

% test=1;
% for i=i:length(grid)-1
%     for j=1:8
%         if func_flux_creen_pt(grid,para(j,1),rho_c,v_m,rho_m)>=func_flux_creen_pt(grid,para(j+1,1),rho_c,v_m,rho_m)
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