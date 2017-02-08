clear;
load GARZ_varyingIC-realistic_rhom0_fval.mat
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

P_sort=sortrows(P,32+7);
rho_m=P_sort(1,32);

grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
para=zeros(9,3);
for i =1:9
    P_sort=sortrows(P,32+i);
    para(i,1)=P_sort(1,(i-1)*3+1);
    para(i,2)=P_sort(1,(i-1)*3+2);
    para(i,3)=P_sort(1,(i-1)*3+3);
    if i~=7
        plot(grid,traffic_flux_smooth(grid,P_sort(1,(i-1)*3+1),P_sort(1,(i-1)*3+2),P_sort(1,(i-1)*3+3),rho_m),'b-','linewidth',2.5),hold on
    else
        plot(grid,traffic_flux_smooth(grid,P_sort(1,(i-1)*3+1),P_sort(1,(i-1)*3+2),P_sort(1,(i-1)*3+3),rho_m),'m-','linewidth',2.5),hold on
    end
end
hold off
test=1;
for i=i:length(grid)-1
    for j=1:8
        if traffic_flux_smooth(i,para(j,1),para(j,2),para(j,3),rho_m)>=traffic_flux_smooth(i,para(j+1,1),para(j+1,2),para(j+1,3),rho_m)
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
para


% for iter =1:10;
%     scatter3(iter,iter,iter),hold on;
% end
% hold off

%         xlabel('initial condition index','fontsize',14);
%         ylabel('\rho_{max}','fontsize',14);
%         strt = ['Parameter-GARZ-rhom ',' \beta = ',num2str(beta(i)),' count=',num2str(count),'/',num2str(882)];
%         strts = ['Parameter-GARZ-rhom-log',num2str(i),'.fig'];
%         title(strt,'fontsize',14);




