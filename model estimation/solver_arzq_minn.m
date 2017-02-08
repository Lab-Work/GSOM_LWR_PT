% function [rho_num,vel_num,erho,evel] = solver_arzq_minn(pos,vm,wm,rhom,h,D,Vel)

function [erho,evel] = solver_arzq_minn(tfinal,pos,vm,wm,rhom,h,D,Vel)

%==========================================================================
%  Simulate AW-Rascal model by Lax-friendrish and HLL approximate Riemann 
%  solver method, In this code, we take the data into classical AR model,
%  that is, the initial and boundary condition are from real traffic data.
%  U = [\rho,u], where rho is the density, and u is the velocity
%  Q = [rho;rho*(u+P(rho))], where P(rho) is pressure function
%  B = density data;
%  Vel = velocity data;
%  [a,b] = time intervel
%  gamma = see Ve function, Ve = vm*(1-(rho/rhom).^gamma);
%  vm, rhom = maximum velocity and density

%  Shimao Fan,
%  Modified in Jan 19th, 2012
%  I try to avoid the difficulty of changing time steps in previous version
%  Temple University, Math Department
%  Last modified on Jan 24 2013
%==========================================================================
% set up dt and dx

L = pos(end);        % length of road in KM
% tfinal = 1.0;     % final time
N = ceil(L/h);    % number of grid
dx = L/N;         % final grid size
% whos
%------------------------------------------
k = .95*h/wm;   % approximate time step
M = ceil(tfinal/k);         % number of time steps
dt = tfinal/M;              % here M is the number of time steps  
lambda = dt/dx;
%------------------------------------------
nofdata = tfinal*120;   % this is number of data points
inter = round(M/nofdata);   % how to pick the data
%------------------------------------------
%--------------
% make a grid
%--------------
x = transpose(dx/2:dx:L-dx/2)'; % stagged grids
t = 0:dt:tfinal;                % time vector, uniform
U=zeros(2,N);      % the solution [density, velocity]
Q=zeros(2,N);      % the unknowns in conserved form
%-----------------------------------
% find position of sensor 2
position = ceil(pos(2:end-1)/dx);     % position of sensor 2
%-----------------------------------
% initial condition of density rho
%-----------------------------------
U(1,:)= 10;     % initial density/ empty
U(2,:) = 100;   % the same initial velocity
%------------------------------------------
for i = 1:length(pos)   % number of sensor
    DensData(i,:) = intpol_bc(D(i,:),t);
    VelData(i,:) = intpol_bc(Vel(i,:),t);
end

% boundary condition
D_lbc = DensData(1,:); V_lbc = VelData(1,:);
D_rbc = DensData(end,:); V_rbc = VelData(end,:);
%*************************************************************************

%//////////////////////
%  The main algorithm
%\\\\\\\\\\\\\\\\\\\\\\
% rho_num = zeros(1,M);
% vel_num = zeros(1,M);
e_rho = zeros(1,M);
e_vel = zeros(1,M);
start = 5*inter;  % this means starting point to compare model
% prediction and real traffic data. Since we just given a initial condition
% which not necessary agree with reality.
Q = unknownQ(U,vm,rhom);

for n = 1:M      % loop over all time
%    Q = unknownQ(U,vm,rhom);
   %----------------------------------------------------------------
   %  Insert new boundary conditions on right hand side
   lbc = [D_lbc(n);V_lbc(n)];
   rbc = [D_rbc(n);V_rbc(n)];
   Up1=[U(:,2:end),rbc];    % right boundary
   Um1=[lbc,U(:,1:end-1)];  % left boundary 
   %------------------------------------------------
   % HLL Riemann solver
%    Q=Q+dt/dx*(HLL(Um1,U,lda,p)-HLL(U,Up1,lda,p));
   Q = Q+lambda*(HLL(Um1,U,vm,rhom)-HLL(U,Up1,vm,rhom));
   %-----------------------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   %-----------------------------------------------------------------
   U(1,:)=Q(1,:);      % traffic density
   U(2,:)=Q(2,:)./Q(1,:)+Ve(U(1,:),vm,rhom);             % traffic velocity
   %----------------------------------------
% we store velocity and density solution and the cooresponding empirical
% data
   rho = U(1,position);
   vel = U(2,position);   
   
%    rho_num(n) = rho;
%    vel_num(n) = vel;
      
   data_rho = DensData(2:end-1,n)';
   data_vel = VelData(2:end-1,n)';
   diff_rho = data_rho-rho;
   diff_vel = data_vel-vel;
   
   e_rho(n) = norm(diff_rho,1);
   e_vel(n) = norm(diff_vel,1);   
%   w = NumVel(n)-Ve(NumDens(n),vm,rhom,gamma)+vm;

%    %----------------------
%    %   The plots part
%    %----------------------
%    inter = 20;
% %    figure(1)
%   if mod(n,inter) == 0
%    subplot(2,1,1)
%    plot(x,U(1,:),'r-','linewidth',3), hold on
%         plot(pos,DensData(n),'bs','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold off
%    title(sprintf('The traffic density at time t=%4.3f, w=%3.1f',n*dt,w),'fontsize',14)
%    axis([0,L,0,400])
%    ylabel('Density','fontsize',14)
% %    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
%    set(gca,'position',[.04 .55 .95 .4])
%    %----------------------------------------------------------------
%    subplot(2,1,2)
%    plot(x,U(2,:),'b-','linewidth',3),hold on
%         plot(pos,VelData(n),'bo','MarkerEdgeColor','k',...
%                 'MarkerFaceColor','r',...
%                 'MarkerSize',10,'linewidth',1), hold off
%    title(sprintf('The traffic velocity at time t=%4.3f, w=%3.1f',n*dt,w),'fontsize',14)
%    axis([0,L,0,120])
%    ylabel('Velocity','fontsize',14)
%    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
%    res = 800;
%    set(gcf,'paperpositionmode','auto')
%    set(gcf,'position',[10  50 res res*.8])
%    set(gca,'position',[.04 .08 .95 .4])
%    drawnow
%    % Print figures
%    filename_save = sprintf('fig_arzq_minn_day17_%06d',n/inter);
%    print(gcf,'-dpng',filename_save,'-r290')
%   end

end

% erho = mean(abs(DensData(start:M)-NumDens(start:M)))/400;
% evel = mean(abs(VelData(start:M)-NumVel(start:M)))/120;
% erho = mean(abs(DensData(start:M)-NumDens(start:M)))/rhom;
% evel = mean(abs(VelData(start:M)-NumVel(start:M)))/vm;
% rho_max = max(maxdens);

nn = length(pos)-2;    % number of sensors validated; the interior ones

erho = mean(e_rho(start:M))/nn;
evel = mean(e_vel(start:M))/nn;

e_rho = e_rho(start:inter:end);
e_vel = e_vel(start:inter:end);
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%---------------------------------------------------------------
function F = HLL(Ul,Ur,vm,rhom) % the HLL Riemann solver
    % where v contains boundary condition of velocity
    % eigenvalues:
    eigns1 = lambda(Ul,vm,rhom);  % left bound
    eigns2 = lambda(Ur,vm,rhom);  % right bound
    % another way to define c1 and c2
    c1 = minimum(eigns1(1,:),eigns2(1,:));
    c2 = maximum(eigns1(2,:),eigns2(2,:));
    case1 = c1>=0;
    case2 = c1<0 & c2>=0;
    case3 = c2<0;
    Ql = unknownQ(Ul,vm,rhom);
    Qr = unknownQ(Ur,vm,rhom);
    % left flux
    Fl = flux(Ul,vm,rhom);  % plug into left boundary velocity information
    % right flux
    Fr = flux(Ur,vm,rhom);
    %-----------------------------------------------
    % intermediate flux
    F_hll(1,:) = (c2.*Fl(1,:)-c1.*Fr(1,:)+(Qr(1,:)-Ql(1,:)).*c1.*c2)./...
            (c2-c1);
    F_hll(2,:) = (c2.*Fl(2,:)-c1.*Fr(2,:)+(Qr(2,:)-Ql(2,:)).*c1.*c2)./...
            (c2-c1);
    F(1,:) = case1.*Fl(1,:)+case2.*F_hll(1,:)+case3.*Fr(1,:);
    F(2,:) = case1.*Fl(2,:)+case2.*F_hll(2,:)+case3.*Fr(2,:);

%--------------------------------------------------------------------------   
function y = lambda(U,vm,rhom)  % U = (density, w), calculate eigenvalues
       diffv = DiffVel(vm,rhom);
       y(1,:) = U(2,:) + U(1,:).*diffv;   % modified here
       y(2,:) = U(2,:);
%----------------------------------------------
function y = flux(U,vm,rhom)  % U = [rho,w], the numerical flux
        ve = Ve(U(1,:),vm,rhom);
        y(1,:) = U(1,:).*U(2,:);
        y(2,:) = U(1,:).*U(2,:).*(U(2,:)-ve);
%--------------------------------------------------------
function y = unknownQ(U,vm,rhom)    % the conserved unknowns [rho, rho*(v-Ve)]
        y(1,:) = U(1,:);     % density
        y(2,:) = U(1,:).*(U(2,:)-Ve(U(1,:),vm,rhom));
%--------------------------------------------------------   
function y = DiffVel(vm,rhom)% define derivative of velocity
        y = -vm/rhom;
%-------------------------------------------------------
function y = Ve(rho,vm,rhom)  % density and 3 parameters 
     y = vm*(1-(rho/rhom));   % which is the quadratic form
%---------------------------------------------------------------- 
function y = minimum(a,b) % here a and b could be vectors
    temp1 = a<=b;
    temp2 = a>b;
    y = temp1.*a+temp2.*b;    % where y is the minimum
  
function y = maximum(a,b) % so as above function
    temp1 = a<=b;
    temp2 = a>b;
    y = temp1.*b+temp2.*a;
%-------------------------------------------------------------   
function y = intpol_bc(data,t)   % cubic
    % data--- traffic data, traffic density
    % t --- time vector
    % assume picewise constant boundary conditions that came from real data
    tfinal = ceil(t(end)*60); % convert hour to minutes
    t = t*60;           % convert hour to minutes
    dt = 1/2;           % average avery 30 seconds, or half minute
    pl = 0:dt:tfinal;   % time interval we interested in for calculation
    n = length(pl);     % number of data we use
    y = spline(pl,data(1:n),t);
%-------------------------------------------------------------
function y = bc(data,t)     % piecewise constant
    % data--- traffic data, traffic density
    % t --- time vector [0:dt:tfinal], time start from zero
    % assume picewise constant boundary conditions that came from real data
    tfinal = t(end)*60; % convert hour to minutes
    t = t*60;           % convert hour to minutes
    dt = 1/2;      % average avery 30 seconds, or half minute
    pl = 0:dt:tfinal;
%     pl = -dt/2:dt:tfinal+dt;  % shift to left by dt/2, that is 1/4 min
    for i = 1:length(pl)-1   % loop over all pices
        c = pl(i)<=t & t<pl(i+1);
        B(i,:) = c*data(i);
    end
%     B
    for j = 1:length(t)
        y(j) = sum(B(:,j));
    end

%%%%%%%%%%%%%%%%%%%%%%End of Code%%%%%%%%%%%%%%%%%%%%%%%

