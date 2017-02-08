
beta = [.001,.01,.03,.05,.1,.3,.5,0.9,0.99];
for i =1:9
    load GPT_varyingIC521.mat
        figure(i)
        rhoc = P_current(1,10);
        vm = P_current(1,11);
        rhom = P_current(1,12);
        count=0;
        for j=1:526;
            strld = ['GPT_varyingIC',num2str(j),'.mat']
            load(strld);
            
            if (P_current(1,10)<rhoc+5 && P_current(1,10)>rhoc-5 && P_current(1,11)<vm+5 && P_current(1,11)>vm-5 && P_current(1,12)<rhom+5 && P_current(1,12)>rhom-5 && P_current(1,i)>-1)
                count=count+1;
%                scatter3(P_current(1,10),P_current(1,11),P_current(1,12),'m'),hold on 
%                set(gca,'Xscale','log','Zscale','log','Yscale','log')
                   scatter(j,P_current(1,i),'m'),hold on
%                   set(gca,'Yscale','log')
            else
%               if (P_current(1,(i-1)*3+1)<10000 && P_current(1,(i-1)*3+2)<1 && P_current(1,(i-1)*3+3)<10000)
%                   scatter3(P_current(1,10),P_current(1,11),P_current(1,12),'b'),hold on
%                   set(gca,'Xscale','log','Zscale','log','Yscale','log')
%               end
% %                  if (P_current(1,32)<10000)
                  scatter(j,P_current(1,i),'b'),hold on
%                   set(gca,'Yscale','log')
% %                  end
            end
        end
        hold off
%         xlabel('\rho_{c}','fontsize',14);
%         ylabel('v_{m}','fontsize',14);
%         zlabel('\rho_{max}','fontsize',14);
%         strt = ['Parameter-GPT ',' \beta = ',num2str(beta(i)),' count=',num2str(count),'/',num2str(526)];
%         strts = ['Parameter-GPT-',num2str(i),'.fig'];
%         title(strt,'fontsize',14);
%         savefig(strts);


        xlabel('initial condition index','fontsize',14);
        ylabel('q','fontsize',14);
        strt = ['Parameter-GPT-q ',' \beta = ',num2str(beta(i)),' count=',num2str(count),'/',num2str(526)];
        strts = ['Parameter-GPT-q',num2str(i),'.fig'];
        title(strt,'fontsize',14);
        savefig(strts);
end

% for iter =1:10;
%     scatter3(iter,iter,iter),hold on;
% end
% hold off

% load Minneapolis_data.mat   % every 30 seconds, V--volume, S---speed
% % %---------------------------------------------------
% % % deal with data and find out density and flux data D and Q
% k = 2;
% Q1 = squeeze(V(:,4*(k-1)+1,:));
% V1 = squeeze(S(:,4*(k-1)+1,:));
% for i = 2:4              % 4 lanes add up
%   Q1 = Q1+squeeze(V(:,4*(k-1)+i,:));
%   V1 = V1+squeeze(S(:,4*(k-1)+i,:)); 
% end
% % change units
% V1 = 1.609344*V1;        % change into km/hour
% Q1 = 120*Q1;             % change into #ofvehicle/hour
% V1 = V1/4;               % average velocity of 4 lanes
% flux = Q1; vel = V1;
% % remove NAN number
% density = flux./vel;
% flux = flux(~isnan(density));
% density = density(~isnan(density));
% figure(1)    % 
% 
% grid = 0:P_current(1,32);    % plot of the sequence of flow-density curves
% plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on
% 
%  plot(grid,traffic_flux_smooth(grid,P_current(1,1),P_current(1,2),P_current(1,3),P_current(1,32)),'b-','linewidth',2.5),hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,4),P_current(1,5),P_current(1,6),P_current(1,32)),'m-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,7),P_current(1,8),P_current(1,9),P_current(1,32)),'r-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,10),P_current(1,11),P_current(1,12),P_current(1,32)),'k-','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,13),P_current(1,14),P_current(1,15),P_current(1,32)),'r-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,16),P_current(1,17),P_current(1,18),P_current(1,32)),'m-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,19),P_current(1,20),P_current(1,21),P_current(1,32)),'b-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,22),P_current(1,23),P_current(1,24),P_current(1,32)),'k-.','linewidth',2.5), hold on
%  plot(grid,traffic_flux_smooth(grid,P_current(1,25),P_current(1,26),P_current(1,27),P_current(1,32)),'b--','linewidth',2.5), hold on 


