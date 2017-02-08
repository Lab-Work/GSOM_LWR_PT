load CGARZ_varyingIC.mat;
beta = [.001,.01,.03,.05,.1,.3,.5,0.9,0.99];
for i =1:9
    if i~=7
        figure(i)
        sigma = P(169,(i-1)*2+1);
        mu = P(169,(i-1)*2+2);
        rhom = P(169,22);
        count=0;
        for j=1:length(P);
            
            if (P(j,(i-1)*2+1)<sigma+50 && P(j,(i-1)*2+1)>sigma-50 && P(j,(i-1)*2+2)<mu+50 && P(j,(i-1)*2+2)>mu-50 && P(j,22)<rhom+50 && P(j,22)>rhom-50)
                count=count+1;
                scatter3(P(j,(i-1)*2+1),P(j,(i-1)*2+2),P(j,22),'m'),hold on             
            else
                scatter3(P(j,(i-1)*2+1),P(j,(i-1)*2+2),P(j,22),'b'),hold on
            end
        end
        hold off
        xlabel('\sigma','fontsize',14);
        ylabel('\mu','fontsize',14);
        zlabel('\rho_{max}','fontsize',14);
        strt = ['Parameter-CGARZ ',' \beta = ',num2str(beta(i)),' count=',num2str(count),'/',num2str(length(P))];
        strts = ['Parameter-CGARZ-',num2str(i),'.fig'];
        title(strt,'fontsize',14);
        savefig(strts);
    end
end

for iter =1:10;
    scatter3(iter,iter,iter),hold on;
end
hold off


