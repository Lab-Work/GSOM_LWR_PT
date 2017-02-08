function plot_gpt_freq_ngsim_100beta_newptfd_fixedthotilde
clear;
load GPT_varyingIC_real_fval_ngsim_newptfd-fixrhotilde
h = 1.0e-03;
beta = linspace(h,1-h,100);
m=length(beta);

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
q_max = max(max(flux));           % ???
figure(1)

P_sort=sortrows(P,9);

rho_m=P_sort(1,3);
v_m=P_sort(1,1);
rho_tilde=P_sort(1,2);
q_eq=P_sort(1,4);

para0=zeros(4,1);

para0(3,1)=P_sort(1,3);
para0(1,1)=P_sort(1,1);
para0(2,1)=P_sort(1,2);
para0(4,1)=P_sort(1,4);


grid = 0:rho_m;    % plot of the sequence of flow-density curves
plot(density,flux,'.','color',[0.2,0.2,0.2],'markersize',6), hold on


para=zeros(m,1);

x01 = [q_eq];
eq_id0=find(beta(:)>0.49 & beta(:)<0.50);
eq_id=max(eq_id0);

for iter=1:m
        
        options = optimoptions('fmincon','MaxFunEvals', 1*10^10, 'MaxIter', 1*10^5, 'TolX',1*10^(-15), 'TolCon',1*10^(-15));
        [x,fval] = fmincon(@(x)beta(iter)*norm(max((func_flux_phase_transition(density,x,v_m,rho_tilde,rho_m)...
            -flux),0*flux))^2+(1-beta(iter))*norm(max(-(func_flux_phase_transition(density,x,v_m,rho_tilde,rho_m)...
            -flux),0*flux))^2,x01,[],[],[],[],[],[],[],options);
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
if iter~=eq_id
    plot(func_flux_phase_transition(grid,para(iter,1),v_m,rho_tilde,rho_m),'b-','linewidth',0.05),hold on
end   
end


plot(func_flux_phase_transition(grid,para(eq_id,1),v_m,rho_tilde,rho_m),'m-','linewidth',2.5),hold on
hold off

savefile = sprintf('GPT_q_0516-ngsim_newptfd');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'para');    % P is the parameters, n cross 3

savefile = sprintf('GPT_q_0516-ngsim-paraeq');
% savefile = sprintf('data_fd_minn_garz');
% 
save(savefile,'para0'); 



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