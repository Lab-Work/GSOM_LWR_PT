function P =plot_cgarz_freq_shrink_removed_data


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

h = 1.0e-03;
%h = 1.0e-03;
beta = linspace(h,1-h,100);
m=length(beta);
eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);
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


load CGARZ_entire_varyingIC-0406-tau300.mat;

index_rf = find(P(:,7)>1);
P=P(index_rf,:);

para_eq=zeros(6,1);


P_sort=sortrows(P,10);
P_sort_rm=zeros(length(P(:,1)),2);
count=1;
for i =1:(length(P(:,1))-1)
    if P_sort(i,10)>P_sort(i+1,10)-2 && P_sort(i,10)<P_sort(i+1,10)+2
        count=count+1;
        P_sort_rm(i,1)=P_sort(i,10);
        P_sort_rm(i,2)=count;
    else
        P_sort_rm(i,1)=P_sort(i,10);
        P_sort_rm(i,2)=count;
        count=1;
    end
        
end
P_sort_rm_end=sortrows(P_sort_rm,2);
rho_m=490.5775;
% rho_m=P_sort_rm_end(length(P(:,1)),1);
num=P_sort_rm_end(length(P(:,1)),2);

index0 = find(P_sort_rm(:,1)==rho_m)+1;
index=max(index0);
P_rm_0=P_sort( [index-num+1:index] , : );

para_eq(6,1)=rho_m;

para0=zeros(1,3);
for j=3:1
    P_sort=sortrows(P_rm_0,6+j);
    P_sort_1=zeros(length(P_rm_0(:,1)),2);
    count=1;
    for i =1:(length(P_rm_0(:,1))-1)
        if P_sort(i,6+j)>P_sort(i+1,6+j)-3 && P_sort(i,6+j)<P_sort(i+1,6+j)+3
            count=count+1;
            P_sort_1(i,1)=P_sort(i,6+j);
            P_sort_1(i,2)=count;
        else
            P_sort_1(i,1)=P_sort(i,6+j);
            P_sort_1(i,2)=count;
            count=1;
        end
        
    end
   j
    P_sort_1_end=sortrows(P_sort_1,2);
    para0(1,j)=P_sort_1_end(length(P_rm_0(:,1)),1);
    num=P_sort_1_end(length(P_rm_0(:,1)),2)
    
    index0 = find(P_sort_1(:,2)==num)+1
    index=max(index0);
    P_rm_0=P_sort([index-num+1:index],:);
    
end


load CGARZ_entire_varyingIC-0406-tau300.mat;

index_rf = find(P(:,7)>1);
P=P(index_rf,:);

para_eq=zeros(6,1);


P_sort=sortrows(P,10);
P_sort_rm=zeros(length(P(:,1)),2);
count=1;
for i =1:(length(P(:,1))-1)
    if P_sort(i,10)>P_sort(i+1,10)-2 && P_sort(i,10)<P_sort(i+1,10)+2
        count=count+1;
        P_sort_rm(i,1)=P_sort(i,10);
        P_sort_rm(i,2)=count;
    else
        P_sort_rm(i,1)=P_sort(i,10);
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

para_eq(6,1)=rho_m;

para0=zeros(1,3);
for j=3:1
    P_sort=sortrows(P_rm_0,6+j);
    P_sort_1=zeros(length(P_rm_0(:,1)),2);
    count=1;
    for i =1:(length(P_rm_0(:,1))-1)
        if P_sort(i,6+j)>P_sort(i+1,6+j)-3 && P_sort(i,6+j)<P_sort(i+1,6+j)+3
            count=count+1;
            P_sort_1(i,1)=P_sort(i,6+j);
            P_sort_1(i,2)=count;
        else
            P_sort_1(i,1)=P_sort(i,6+j);
            P_sort_1(i,2)=count;
            count=1;
        end
        
    end
   j
    P_sort_1_end=sortrows(P_sort_1,2);
    para0(1,j)=P_sort_1_end(length(P_rm_0(:,1)),1);
    num=P_sort_1_end(length(P_rm_0(:,1)),2)
    
    index0 = find(P_sort_1(:,2)==num)+1
    index=max(index0);
    P_rm_0=P_sort([index-num+1:index],:);
    
end

rho0=63.8473;
uf=100.8327;
rhof=901.6230;
% rho0=para0(1,1);
% uf=para0(1,2);
% rhof=para0(1,3);

para_eq(3,1)=rho0;
para_eq(4,1)=uf;
para_eq(5,1)=rhof;


grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on



paraeq=zeros(1,2);

for iter=2:2
    P_rm=P_rm_0;
    if length(P_rm(:,1))~=1
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
            j
            P_sort_1_end=sortrows(P_sort_1,2)
            paraeq(1,j)=P_sort_1_end(length(P_rm(:,1)),1);
            
        end

            plot(grid,func_fd_seibold_2p_updated(grid,paraeq(1,1),paraeq(1,2),rho0,uf,rhof,rho_m),'m-','linewidth',2.5),hold on
   
    else
        for j=1:2
            paraeq(1,j)=P_rm(1,(iter-1)*2+j);
        end

            plot(grid,func_fd_seibold_2p_updated(grid,paraeq(1,1),paraeq(1,2),rho0,uf,rhof,rho_m),'m-','linewidth',2.5),hold on

        
    end
end

para_eq(1,1)=6.4760;
para_eq(2,1)=83.0553;

% para_eq(1,1)=paraeq(1,1);
% para_eq(2,1)=paraeq(1,2);


lb = [1.0e-03,rho0];
ub = [1000,2000]; 
lb1 = [1.0e-03, rho0];
ub1 = [1000,2000]; 

count=0;

%YE CGARZ-------------------------------------------------------
% initials = [1,10,50,100,500,1000];
% initialm = [50,100,500,1000];
initials = [paraeq(1,1)];
initialm = [paraeq(1,2)];
% initials = [100];
% initialm = [100];
% initialr = [500,600];
N1 = length(initials);
N2=length(initialm);
N=N1*N2;
m = length(beta);
P2 = zeros(m,2);  % initialize the matrix
tau=300;

for iter1 = 1:N1
    for iter2 = 1:N2       
        count=count+1;       
        x01 = [initials(iter1) initialm(iter2)];     
        for j = 1:m
            j
            
                %        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
                %        options = optimset('MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5);
                %options = optimset('MaxFunEvals', 1*10^10);
                                
                [x,fval] = fmincon(@(x) beta(j)*(norm(max(shrinkage(density,(func_fd_seibold_2p_updated(density,x(1),x(2),rho0,uf,rhof,rho_m)...
                    -flux),rho0,tau),0*flux)))^2+(1-beta(j))*(norm(max(shrinkage(density,-(func_fd_seibold_2p_updated(density,x(1),x(2),rho0,uf,rhof,rho_m)...
                    -flux),rho0,tau),0*flux)))^2,x01,[],[],[],[],lb1,ub1);%CGARZ!!!

                               
                P2(j,1) = x(1);
                P2(j,2) = x(2);
                count;
                
                plot(grid,func_fd_seibold_2p_updated(grid,P2(j,1),P2(j,2),rho0,uf,rhof,rho_m),'b-','linewidth',0.05),hold on
                      
        end
        
        strr = ['\rho_{f} = ',num2str(rho0)];
        text(rho_m-100,10000,strr,'fontsize',14);
        strrm = ['\rho_{m} = ',num2str(rho_m)];
        text(rho_m-100,8000,strrm,'fontsize',14);
    end    
end

savefile = sprintf('CGARZ_100beta_2para_updated');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'P2');    % P is the parameters, n cross 3

savefile = sprintf('CGARZ_eq_6para_updated');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'para_eq');    % P is the parameters, n cross 3




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