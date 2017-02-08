function [num_rho,num_vel,erho,evel] = solver_gpt_green_ngsim_ctm(tfinal,D,V,BD,BV,...
    x,ref,q_min,q_max,rhoc,vm,rhom)
%==========================================================================
% Generalized Aw-Rascal(GAR) model, with fundamental diagrams families has 3
% free parameters(TP), in data from California(CA)
% Riemann solver, HLL;
% fundamental diagram fitted from macroscopic data;
% initial density and boundary conditions are formed from NGSIM data
% Shimao Fan,
% Temple University, Math department,
% Fet 7th, 2012
%***********************************
% D === NGSIM_DensData
% V ==== NGSIM_VelData
% cl, cp and ca are parameters of fundamental diagrams
% h is size of grid point
% rhom is the maximum density
% n1 is average of time for interior data
% n2 is for boundary
% x0 is space vector choose h = 32 feet
% hom = 0 ..... homogeneous model
% hom = 1 ..... with stiff relaxiation terms
%***********************************
% Add relaxiation term on right hand side.
% Shimao Fan, March 15 2012
% Temple University.
%==========================================================================
% close all;
% some parameters
% pick the middle curve/ NGSIM data
% lda = P(1);
% p = P(2);
% alpha = P(3);
% rhovec = .1:.1:rhom;           % traffic density vector       
%------------------------------------------------
% determine dx and dt
% tfinal = 14;         % mins
h = x(2)-x(1);
N = length(x);       % length of x vector
k = 1/36000;         % length of time step, does not work start from ref = 8
t0 = 0:k:(tfinal/60);
% m1 = length(t0);
%-------------------------------------------
    dt = 2*k/ref; 
    t = 0:dt:(tfinal/60);
    D_lbc = spline(t0,BD(:,1),t);
    D_rbc = spline(t0,BD(:,2),t);
    V_lbc = spline(t0,BV(:,1),t);
    V_rbc = spline(t0,BV(:,2),t);
    % perform projection/modification here.
    index_l = find(D_lbc>rhom);
    index_r = find(D_rbc>rhom);
    D_lbc(index_l) = .9*rhom;
    V_lbc(index_l) = 1.0;
    D_rbc(index_r) = .9*rhom;
    V_rbc(index_r) = 1.0;
    % interpol respect to time for all the data
    for i = 1:size(D,2)   % loop respect space grids
        D_data(:,i) = spline(t0,D(:,i),t);
        V_data(:,i) = spline(t0,V(:,i),t);
    end
%--------------------------------------------    
% at this end, we find time step satisfy CFL condition
M = ceil(tfinal/(60*dt));
lambda = dt/h;        % the length of road
% L = x(end);      
% pos = ceil(.224/h);
%------------------------------------------------
% initial condition
icv = V(1,:);
icd = D(1,:);
index = find(icd>rhom);
icd(index) = .9*rhom;   % also do projection to initial data
icv(index) = 1.0;
U(1,:) = icd;
% calculate initial w vector
for i = 1:N
   U(2,i) = func_inverse_v_gpt(icd(i),icv(i),q_min,q_max,rhoc,vm,rhom);
end

%*************************************************************************
%*******************  ----------------------  ****************************
%                       The main algorithm
%*******************  ----------------------  ****************************
%*************************************************************************

erho = zeros(1,M);
evel = zeros(1,M);
% Q = unknownQ(U);

% w_min = 55;
% w_max = 80;
% num_rho = zeros(M,N);
% num_vel = zeros(M,N);
midd = round(N/2);   % middle position of study area
%=================
%U = [rho,q], the two conserved quantities
%=================
for n = 1:M      % loop over all time U = [density, w]; as unknowns
   % represent Q as function of U, the conservative unknowns
%    Q(1,:)=U(1,:);
%    Q(2,:)=Q(1,:).*U(2,:);
   %----------------------------------------------------------------------
   %  Insert new boundary conditions on right hand side
   rhol = D_lbc(n); rhor = D_rbc(n); % density of left and right bound
   vell = V_lbc(n); velr = V_rbc(n);
%    diffl = abs(rhol-rhovec); diffr = abs(rhor-rhovec);
%    [C,indl] = min(diffl); [C,indr] = min(diffr);  % index of the density
%    uvecl = squeeze(Mw(indl,2,:));   
%    uvecr = squeeze(Mw(indr,2,:)); % velocity vector
%    wvec = squeeze(Mw(indl,1,:));                          % w vector
%    dvl = abs(V_lbc(n)-uvecl); dvr = abs(V_rbc(n)-uvecr);
%    [C,indl2] = min(dvl);    [C,indr2] = min(dvr);  % row of the velocity
%    wl = wvec(indl2);    wr = wvec(indr2);          % find correaponding w values

   % apply bisection root finder
   wl = func_inverse_v_gpt(rhol,vell,q_min,q_max,rhoc,vm,rhom);
   wr = func_inverse_v_gpt(rhor,velr,q_min,q_max,rhoc,vm,rhom);
   %--------------------------------------
% %    % apply newton's iteration
%     wl = w_max;
%     wr = wl;    % initial guess    
%     for i =1 : 2   % something wrong with this root finding process        
%        wl = wl+(vell-Velocity(rhol,wl,cl,cp,ca,rhom))./...
%            func_diff_v_w(rhol,wl,cl,cp,ca,rhom);
%        wr = wr+(velr-Velocity(rhor,wr,cl,cp,ca,rhom))./...
%           func_diff_v_w(rhor,wr,cl,cp,ca,rhom);       
%     end  
%     
%    wl = min(max(wl,w_min),w_max);
%    wr = min(max(wr,w_min),w_max); 
   %--------------------------------------
   lbc = [D_lbc(n);wl];
   rbc = [D_rbc(n);wr];
   Up1=[U(:,2:end),rbc];             % right boundary
   Um1=[lbc,U(:,1:end-1)];           % left boundary 

   U = U+ lambda*(Num_Flux(Um1,U,rhoc,vm,rhom)-Num_Flux(U,Up1,rhoc,vm,rhom));
   %----------------------------------------------------------------------
   %  find [rho,u] from conserved form, solution of next time step
   %-----------------------------------------------------------------
%    U(1,:)=Q(1,:);                     % traffic density
%    U(2,:)=Q(2,:)./Q(1,:);             % traffic velocity
   %----------------------------------------
   rho = U(1,:);
   q = U(2,:);
%    w = projection(w,w_min,w_max);
   vel = Velocity(rho,q,rhoc,vm,rhom);
   %----------------------------------------
   num_rho(n) = rho(midd);
   num_vel(n) = vel(midd);
   
   %////////////////////
   %   The plots part
   %\\\\\\\\\\\\\\\\\\\\
%   if mod(n,10) == 0
%    subplot(2,1,1)
%    plot(x,U(1,:),'m-','linewidth',2.5), hold on
%    plot(x,D(n,:),'b-','linewidth',2.5), hold off
%    title(sprintf('The traffic density at time t=%4.3f',n*dt))
%    axis([0,x(end),0,600])
%    ylabel('Density','fontsize',14)
%    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
% %    %----------------------------------------------------------------
%    subplot(2,1,2)
%    plot(x,vel,'m-','linewidth',2.5),hold on
%    plot(x,V(n,:),'b-','linewidth',2.5), hold off
% %    title(sprintf('The traffic velocity at time t=%4.3f',n*dt))
%    axis([0,L,0,100])
%    ylabel('Velocity','fontsize',14)
%    xlabel('Space','fontsize',14)
%    set(gca,'linewidth',2)
%    set(gca,'fontsize',14)
%    res = [971 600];
%    set(gcf,'paperpositionmode','auto')
%    set(gcf,'position',[10 50 res(1) res(2)])
%    drawnow
% %    Print figures
% %    filename_save = sprintf('AwRscal_%06d',n/40);
% %    print(gcf,'-dpng',filename_save,'-r200')
%   end

%------------------------------------------------
% Calculate error 
%    [e1,e2] = comput_error_nagim_time(D_data(n,:),...
%              V_data(n,:),rho,vel);
%    erho(n) = e1;
%    evel(n) = e2;
   rho_diff = D_data(n,:)-rho;
   vel_diff = V_data(n,:)-vel;
   e1 = norm(rho_diff(1:end),1);
   e2 = norm(vel_diff(1:end),1);
   erho(n) = e1;
   evel(n) = e2;
   
%    U(:,1) = lbc;
%    U(:,end) = rbc;
%    
   % boundary adjustment
%    % rescale with data
%    es1 = norm((D_data(n,:)-rho)./D_data(n,:),1);
%    es2 = norm((V_data(n,:)-vel)./V_data(n,:),1);
%    
%    erhos(n) = es1;
%    evels(n) = es2;
%    n
end

%---------------------------------
% rhom = 700;
% vm = 75;
% erho = erho/ref/rhom;
% evel = evel/ref/vm;
% erho = erho/ref/700;
% evel = evel/ref/75;
erho = erho/ref;
evel = evel/ref;
%\\\
% erhos = erhos/ref;
% evels = evels/ref;

%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================
% function y = projection(w,w_min,w_max)  % given w vector, this function 
%         % will perform projection. Which means for w outside of w intervel,
%         % that is w below w_min and above w_max, project the value to the
%         % upper and lower bounds.
%         case1 = w<w_min;
%         case2 = w>=w_min & w<=w_max;
%         case3 = w>w_max;
%         y = w_min*case1+w.*case2+w_max*case3;
%---------------------------------------------------------
function y = func_diff_v_w(rho,w,cl,cp,ca,rhom)
     eps = 1.0e-06;
     dy = Velocity(rho,w+eps,cl,cp,ca,rhom) - ...
            Velocity(rho,w,cl,cp,ca,rhom);
     y = dy/eps;
%---------------------------------------------------------
function  F = Num_Flux(Ul,Ur,rhoc,vm,rhom) % numerical flux
   % find sending function
   s_flux = sending(Ul(1,:),Ul(2,:),rhoc,vm,rhom);
   % receiving of flux
   r_flux = receiving(Ul,Ur,rhoc,vm,rhom);
   
   w = Ul(2,:)./Ul(1,:);
   % conservation of mass
   F(1,:) = min(s_flux,r_flux);   % take the minimum
   % conservation of momentum
   F(2,:) = w.*F(1,:);
%----------------------------------------------
% define the receiving function
function y = receiving(Ul,Ur,rhoc,vm,rhom)
      % determined by Ul
      w = Ul(2,:)./Ul(1,:);
      % velocity on the right side
      ur = Velocity(Ur(1,:),Ur(2,:),rhoc,vm,rhom);
      % newton iteration is faster/less costy
      rho_middle = Riemann_G(Ul,Ur,rhoc,vm,rhom); % middle state
      q = w.*rho_middle;
      sigma = func_rhoc(rhoc,q,vm,rhom);  % a vector of rhoc
      case1 = rho_middle<=sigma;  % maximum 
      case2 = rho_middle>sigma;   % maximum potential
      q_max = traffic_flux_smooth(sigma,q,rhoc,vm,rhom); % maximum flow rate
      y = case2.*rho_middle.*ur + case1.*q_max;  
%----------------------------------------------
% define sending functions
function y = sending(rho,q,rhoc,vm,rhom)% U = [\rho, w] two state variables

      sigma = func_rhoc(rhoc,q,vm,rhom);  % a vector of rhoc/3-params fd
      case1 = rho<=sigma;  % maximum possible
      case2 = rho>sigma;   % maximum
      q_max = traffic_flux_smooth(sigma,q,rhoc,vm,rhom); % maximum flow rate
      q = traffic_flux_smooth(rho,q,rhoc,vm,rhom);
      y = case1.*q + case2.*q_max;
%-----------------------------------------------
% calculate the middle density
function rho = Riemann_G(Ul,Ur,rhoc,vm,rhom) % calculate the riemann solution/intermediate state
 % Ul and Ur are initial data of Riemann problem
  u = Velocity(Ur(1,:),Ur(2,:),rhoc,vm,rhom); % right velocity
  w = Ul(2,:)./Ul(1,:);     % left propoty
  rho = Ul(1,:);   % not sure which one is a good pick
  for i =1 : 2   % something wrong with this root finding process
       rho = rho+(u-Velocity(rho,rho.*w,rhoc,vm,rhom))./...
           DiffVel3(rho,rho.*w,rhoc,vm,rhom);     
  end
%--------------------------------------------------------
% function y = unknownQ(U)    % the conserved unknowns [rho, rho*w]
%         y(1,:) = U(1,:);
%         y(2,:) = U(1,:).*U(2,:);
%----------------------------------------------------   
function y = DiffVel3(rho,q,rhoc,vm,rhom)% define derivative of velocity
   eps = 1.0e-3;
   dy = func_vel_creen_pt(rho+eps,q,rhoc,vm,rhom)-...
       func_vel_creen_pt(rho,q,rhoc,vm,rhom);
   y = dy/eps;
%----------------------------------------------------
function f = traffic_flux_smooth(rho,q,rhoc,vm,rhom)
         v = func_vel_creen_pt(rho,q,rhoc,vm,rhom);
         f = rho.*v;
%----------------------------------------------------         
function y = Velocity(rho,q,rhoc,vm,rhom) % velocity function
       % x is density;   
   y = func_vel_creen_pt(rho,q,rhoc,vm,rhom); 
%--------------------------------------------------------------
% calculate the critical denisty
function y = func_rhoc(rhoc,q,vm,rhom)
        a = -rhoc.*vm./((2*rhoc-3*rhom).^2);
        ca = a;
        cb = -a.*(rhom+rhoc) + rhoc.*vm./(rhoc-rhom)-vm./(1+q);
        cc = a.*rhom.*rhoc -rhoc.*vm.*rhom./(rhoc-rhom);
             
        y = (-cb - sqrt(cb.^2-4*ca.*cc))./(2*ca); 
%---------------------------------------------------------------- 
% function y = minimum(a,b) % here a and b could be vectors
%     temp1 = a<=b;
%     temp2 = a>b;
%     y = temp1.*a+temp2.*b;    % where y is the minimum
%   
% function y = maximum(a,b) % so as above function
%     temp1 = a<=b;
%     temp2 = a>b;
%     y = temp1.*b+temp2.*a;

%%%%%%%%%%%%%%%%%%%%%%End of Code%%%%%%%%%%

