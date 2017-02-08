function P = gen_data_fd_minn_entire_cgarz_fval_shrink_removed_data


%//////////////////////////////////////////////////////////////////////////
% This code is used to generate FD data from wrighted least squre fitting
% processes. This code can generate fd parameter either with 2 free
% parameters and 3 free parameters in the order [lambda, p, alpha].
% Shimao Fan, Temple University
% Jan 23 2012

% calculate fd for different days for one given sensor station
% k denotes sensor station
% sensor station 1~6 have 4 lanes
% sensor station 7~10 have 3 lanes
% 3 on ramp sensors and 3 off ramp sensor
% for 4 lanes piece of road
%//////////////////////////////////////////////////////////////////////////
%Ye Sun
%RELAX RHOM
%//////////////////////////////////////////////////////////////////////////
% close all;
% clc;
% clear all;
%---------------------------------------------------
% The maximum traffic density (maximum possible vel can fit into given 
%---------------------------------------------------
% load data
close all;
% Loading detector data part
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

% load data_minn_rho.mat   % density data
% load data_minn_flux.mat  % flux data
%==========================================
% The weighted least squre fitting process
%==========================================

h = 1.0e-02;
% h = 1.0e-03;
% h = .01;
% beta = h/2:h:1-h/2;
% beta = .5;
%m = 100;
% beta = linspace(h,1-h,2);
% beta = [h,.5,1-h];
% % test beta data
beta = [.2,.001,.7,.01,.3,.99,0.5,0.999,0.8];
% beta = [.000015,0.9,0.99,0.999,.9999];
%------------------------------------------
%YE GARZ-------------------------------------------------------
%  3 parameters square root fd (original one)
% lambda0=20;
% p0=0.5;
% alpha0=1000;
% rhom0=10;
% x0 = [lambda0,p0,alpha0,rhom0];   % (lambda, p, alpha,rhom);
% x1 = [lambda0,p0,alpha0];   % (lambda, p, alpha);
% beta = linspace(h,1-h,100);
%YE GARZ-------------------------------------------------------

%---------------------------------
% 3 parameter arctangent fd
% x0 = [10,100,10];
% P = zeros(m,3);  % initialize the matrix
% lb = [.1, 50 ,1];
% ub = [1000,rhom, 100]; 
%---------------------------------
% 2-p collapsed model, arctangent
% beta = linspace(h,1-h,100);
% 
% % uf = 114; rhof = 1026; 
% % uf = 104.4861; rhof = 566.9842;
% % rho0 = 48.5;
% % rho0 = 33;
%YE CGARZ-------------------------------------------------------


lb = [1.0e-03,0,1.0e-03,0,1.0e-03,0,0,50,0,0];
ub = [1000,1000,1000,1000,1000,1000,1000,150,2000,1000]; 

%YE CGARZ-------------------------------------------------------
%--------------------------------------
% newell-danganzo model
% x0 = [100 rhom/4];
% beta = 0.5;
%--------------------------------------
% free flow branch of cgarz model
% vm = 93.6380;   rhoc = 81.6558;
% x0 = [vm rhom];
% beta = 0.5;
% index = find(density<=rhoc);
% D_free = density(index);
% Q_free = flux(index);
%---------------------------------------
% determine standard curve in gpt
% x0 = [rhom/4 95];
% beta = 0.5;
%-------
% a family of perturbed curves
% rhoc = 80.6972;   vm = 93.7165;
% beta = [h 1-h];
% x0=0;
%-----------------------------------
% a quadratic flux curve
%x0 = 90;
%beta = 0.5;
%---------------------------------
%m = length(beta);
%P = zeros(m,1);  % initialize the matrix


%---------------------------------

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
%YE GARZ-------------------------------------------------------
% options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
% %options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
%  lbeq = [0,0,0,0];
% %options = optimset('MaxFunEvals', 1*10^10);
%  [x,fval] = fmincon(@(x)beta(7)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
%     -flux),0*flux))^2+(1-beta(7))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
%     -flux),0*flux))^2,x0,[],[],[],[],lbeq,[],[],options);
% %  [x,fval] = fminsearch(@(x)beta(7)*norm(max((traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
% %     -flux),0*flux))^2+(1-beta(7))*norm(max(-(traffic_flux_smooth(density,x(1),x(2),x(3),x(4))...
% %     -flux),0*flux))^2,x0,options);
% 
%   P(7,1) = x(1);
%   P(7,2) = x(2);
%   P(7,3) = x(3);
%   
%   rho_m = x(4);
%   P(7,1)
%   P(7,2)
%   P(7,3)
%   rho_m
%YE GARZ-------------------------------------------------------

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

count=0;

%YE CGARZ-------------------------------------------------------
initials = [5,10];
initialm = [50,100];
initialr0 = [25,50];
initialuf=[50,100];
initialrf=[800,1200];
initialr = [500,800];
% initials = [50];
% initialm = [100];
% initialr0 = [25];
% initialuf=[50];
% initialrf=[800];
% initialr = [500];
N1 = length(initials);
N2=length(initialm);
N3=length(initialr0);
N4 = length(initialuf);
N5=length(initialrf);
N6=length(initialr);
N=N1*N2*N3*N4*N5*N6;
m = length(beta);
P = zeros(N,2*(m+2)+6+1);  % initialize the matrix
% P = zeros(m,3);  % initialize the matrix
%YE CGARZ-------------------------------------------------------

%YE GARZ-----------------------------------------------------------------
% initiall = [1,50,100,500,1000,4000];
% initialp = [0.1,0.5,0.9];
% initiala = [1,50,100,500,1000,4000];
% initialr = [400,600,800];
% % initiall = [50];
% % initialp = [0.5];
% % initiala = [1000];
% % initialr = [600,400];
% lbeq = [0,0,0,0];
% lb1 = [0,0,0];
% N1=length(initiall);
% N2=length(initialp);
% N3=length(initiala);
% N4=length(initialr);
% N=N1*N2*N3*N4;
% m = length(beta);
% P = zeros(N,3*(m)+5+m);  % initialize the matrix
% P_current = zeros(1,3*(m)+5+m)
% % P = zeros(m,3);  % initialize the matrix
%YE GARZ-----------------------------------------------------------------
tau=500;%%currently tau=300 is the best
%y1=shrinkage(3,5,6,10)
% YE CGARZ----------------------------------------------------------------
for iter1 = 1:N1
    for iter2 = 1:N2
        for iter3 = 1:N3
            for iter4 =1:N4
                for iter5=1:N5
                    for iter6=1:N6
                        count=count+1;
                        count
                        
                        %             x0 = [initials(iter1) initialm(iter2) initialr(iter3)];
                        x0=[initials(iter1) initialm(iter2) initials(iter1) initialm(iter2) initials(iter1) initialm(iter2) initialr0(iter3) initialuf(iter4) initialrf(iter5) initialr(iter6)];
                        %            x01 = [initials(iter1) initialm(iter2)];
                        options = optimoptions('fmincon','TolFun', 1*10^(-15),'TolX',1*10^(-15));
                        %options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
                        %options = optimset('MaxFunEvals', 1*10^10);
                        %             [x,fval] = fmincon(@(x) beta(9)*(norm(max((func_fd_seibold_2p_updated(density,x(1),x(2),x(3),x(4),x(5),x(6))...
                        %                 -flux),0*flux)))^2+(1-beta(9))*(norm(max(-(func_fd_seibold_2p_updated(density,x(1),x(2),x(3),x(4),x(5),x(6))...
                        %                 -flux),0*flux)))^2,x0,[],[],[],[],lb,ub);%CGARZ!!!
                        
                        
                        [x,fval] = fmincon(@(x) beta(9)*(norm(max(shrinkage(density,(func_fd_seibold_2p_updated(density,x(1),x(2),x(7),x(8),x(9),x(10))...
                            -flux),x(7),tau),0*flux)))^2+(1-beta(9))*(norm(max(shrinkage(density,-(func_fd_seibold_2p_updated(density,x(1),x(2),x(7),x(8),x(9),x(10))...
                            -flux),x(7),tau),0*flux)))^2+beta(7)*(norm(max(shrinkage(density,(func_fd_seibold_2p_updated(density,x(3),x(4),x(7),x(8),x(9),x(10))...
                            -flux),x(7),tau),0*flux)))^2+(1-beta(7))*(norm(max(shrinkage(density,-(func_fd_seibold_2p_updated(density,x(3),x(4),x(7),x(8),x(9),x(10))...
                            -flux),x(7),tau),0*flux)))^2+beta(1)*(norm(max(shrinkage(density,(func_fd_seibold_2p_updated(density,x(5),x(6),x(7),x(8),x(9),x(10))...
                            -flux),x(7),tau),0*flux)))^2+(1-beta(1))*(norm(max(shrinkage(density,-(func_fd_seibold_2p_updated(density,x(5),x(6),x(7),x(8),x(9),x(10))...
                            -flux),x(7),tau),0*flux)))^2,x0,[],[],[],[],lb,ub,[],options);%CGARZ!!!
                        
                        P(count,m*2+4+1) = initials(iter1);
                        P(count,m*2+4+2) = initialm(iter2);
                        P(count,m*2+4+3) = initialr0(iter3);                                               
                        P(count,m*2+4+4) = initialuf(iter4);
                        P(count,m*2+4+5) = initialrf(iter5);
                        P(count,m*2+4+6) = initialr(iter6);
                        P(count,m*2+4+7) = fval;
                        
                        P(count,17) = x(1);
                        P(count,18) = x(2);
                        P(count,19) = x(7);
                        P(count,20) = x(8);
                        P(count,21) = x(9);
                        P(count,22) = x(10);
                        P(count,1) = x(5);
                        P(count,2) = x(6);
                        P(count,13) = x(3);
                        P(count,14) = x(4);

                    end
                end
            end
        end
    end   
end
%YE CGARZ-------------------------------------------------------------



% fd_equi = P;
% rhoc = rho_critical(fd_equi(1),fd_equi(2),rhom)

% data_fd_ngsim_p3 have 3 parameters
% savefile = sprintf('data_fd_pt_minn_%1.0f',rhom);
%strsv = ['rhom0_',num2str(rhom0),'sigma0_',num2str(sigma0),'mu0_',num2str(mu0)]

% 
% % load data_fd_minn_arctan.mat
% 
% 
% %\\\\\\\\\\\\\\\\\\
% %  The plots part
% %//////////////////
% 
% %figure(1)    % 
% 
% 
% % plot(d_s1,f_s1,'r.'), hold on
% % plot(d_s3,f_s3,'b*'), hold on

P_sort=sortrows(P,29);

savefile = sprintf('CGARZ_entire_fval-rhof-1-minn-0126-tau500');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'P_sort');    % P is the parameters, n cross 3
rho_m=P_sort(1,22);

grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on

rho0 = P_sort(1,19);; uf = P_sort(1,20);;  rhof = P_sort(1,21);
v0 = diff_fd_free(rho0,uf,rhof);
f0 = fd_free(rho0,uf,rhof);


plot(grid,func_fd_seibold_2p_updated(grid,P_sort(1,1),P_sort(1,2),rho0,uf,rhof,rho_m),'r-','linewidth',2.5),hold on
plot(grid,func_fd_seibold_2p_updated(grid,P_sort(1,17),P_sort(1,18),rho0,uf,rhof,rho_m),'r-','linewidth',2.5),hold on
plot(grid,func_fd_seibold_2p_updated(grid,P_sort(1,13),P_sort(1,14),rho0,uf,rhof,rho_m),'m-','linewidth',2.5),hold on

lb = [1.0e-03,rho0];
ub = [1000,2000]; 
lb1 = [1.0e-03, rho0];
ub1 = [1000,2000]; 

count=0;

%YE CGARZ-------------------------------------------------------
initials = [50,100];
initialm = [50,100];
% initials = [100];
% initialm = [100];
% initialr = [500,600];
N1 = length(initials);
N2=length(initialm);
N=N1*N2;
m = length(beta);
P2 = zeros(N,3*9+2);  % initialize the matrix

for iter1 = 1:N1
    for iter2 = 1:N2       
        count=count+1;       
        x01 = [initials(iter1) initialm(iter2)];     
        for j = 1:m         
            if j~=7 && j~=1 && j~=9
                %        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
                %        options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
                %options = optimset('MaxFunEvals', 1*10^10);
                
                [x,fval] = fmincon(@(x) beta(j)*(norm(max(shrinkage(density,(func_fd_seibold_2p_updated(density,x(1),x(2),rho0,uf,rhof,rho_m)...
                    -flux),rho0,tau),0*flux)))^2+(1-beta(j))*(norm(max(shrinkage(density,-(func_fd_seibold_2p_updated(density,x(1),x(2),rho0,uf,rhof,rho_m)...
                    -flux),rho0,tau),0*flux)))^2,x01,[],[],[],[],lb1,ub1);%CGARZ!!!
               
                
                
                               
                P2(count,(j-1)*3+1) = x(1);
                P2(count,(j-1)*3+2) = x(2);
                P2(count,(j-1)*3+3)=fval;
                P2(count,28)=initials(iter1);
                P2(count,29)=initialm(iter2);
                count
                j
            end           
        end
    end    
end

savefile = sprintf('CGARZ_entire_fval-rhof-2-minn-0126-tau500');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'P2');    % P is the parameters, n cross 3

para=zeros(9,2);

for i =1:9
    P_sort=sortrows(P2,3*i);
    para(i,1)=P_sort(1,(i-1)*3+1);
    para(i,2)=P_sort(1,(i-1)*3+2);
    if i~=7 && i~=1 && i~=9
        plot(grid,func_fd_seibold_2p_updated(grid,P_sort(1,(i-1)*3+1),P_sort(1,(i-1)*3+2),rho0,uf,rhof,rho_m),'b-','linewidth',2.5),hold on
    end
    
end
hold off
para
a1=max(size(density))






strr = ['\rho_{f} = ',num2str(rho0)];
text(rho_m-100,10000,strr,'fontsize',14);

function y = fd_free(rho,uf,rhof)  % quadratic form

y = uf*rho.*(1-rho/rhof);

function y = diff_fd_free(rho,uf,rhof)

y = uf*(1-2*rho/rhof);

function y = linear_bound(rho,v0,rho0,f0)

y = f0+v0*(rho-rho0);

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


function y = shrinkage(rho,d,rho0,tau)
case1 = rho<=rho0 & d<=tau & d>=-tau;
case2 = ~case1;
y = case1.*0+case2.*d;