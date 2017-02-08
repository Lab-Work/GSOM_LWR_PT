load CGARZ_varyingIC_eq.mat;
beta = [.001,.01,.03,.05,.1,.3,.5,0.9,0.99];

        figure(1)
        sigma = P(169,1);
        mu = P(169,2);
        rhom = P(169,3);
        count=0;
        for j=1:length(P);
            
            if (P(j,1)<sigma+50 && P(j,1)>sigma-50 && P(j,2)<mu+50 && P(j,2)>mu-50 && P(j,3)<rhom+50 && P(j,3)>rhom-50)
                count=count+1;
                scatter3(P(j,1),P(j,2),P(j,3),'m'),hold on             
            else
                scatter3(P(j,1),P(j,2),P(j,3),'b'),hold on
            end
        end
        hold off
        xlabel('\sigma','fontsize',14);
        ylabel('\mu','fontsize',14);
        zlabel('\rho_{max}','fontsize',14);
        strt = ['Parameter-CGARZ ',' \beta = ',num2str(0.5),' count=',num2str(count),'/',num2str(length(P))];
        strts = ['Parameter-CGARZ-eq','.fig'];
        title(strt,'fontsize',14);
        savefig(strts);





