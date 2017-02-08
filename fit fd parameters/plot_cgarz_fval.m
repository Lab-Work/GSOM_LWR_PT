function plot_cgarz_fval
clear;
load CGARZ_varyingIC-realistic_rhom0_fval.mat
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

P_sort=sortrows(P,29);
rho_m=P_sort(1,22);

grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
rho0 = 36.2; uf = 93.7124;  rhof = 529.9549;
v0 = diff_fd_free(rho0,uf,rhof);
f0 = fd_free(rho0,uf,rhof);
para=zeros(9,2);

for i =1:9
    P_sort=sortrows(P,22+i);
    para(i,1)=P_sort(1,(i-1)*2+1);
    para(i,2)=P_sort(1,(i-1)*2+2);
    if i~=7
        plot(grid,func_fd_seibold_2p(grid,P_sort(1,(i-1)*2+1),P_sort(1,(i-1)*2+2),rho0,v0,f0,uf,rhof,rho_m),'b-','linewidth',2.5),hold on
    else
        plot(grid,func_fd_seibold_2p(grid,P_sort(1,(i-1)*2+1),P_sort(1,(i-1)*2+2),rho0,v0,f0,uf,rhof,rho_m),'m-','linewidth',2.5),hold on
    end
end
hold off

test=1;
for i=i:length(grid)-1
    for j=1:8
        if func_fd_seibold_2p(i,para(j,1),para(j,2),rho0,v0,f0,uf,rhof,rho_m)>=func_fd_seibold_2p(i,para(j+1,1),para(j+1,2),rho0,v0,f0,uf,rhof,rho_m)
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

function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);

function y = diff_fd_free(rho,uf,rhof)

y = uf*(1-2*rho/rhof);




