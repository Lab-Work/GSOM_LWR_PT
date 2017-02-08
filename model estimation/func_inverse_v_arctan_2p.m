function w = func_inverse_v_arctan_2p(rho,u,q_min,q_midd,q_max,ca,cb,...
    rho0,v0,f0,uf,rhof,rhom,proj)

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
% For given rho and u, this function try to identify which curve this point
% belongs to. I would like to solve the problem in a flow-density plane,
% that is (rho, u*rho), and find the corresponding 'w' by the Newton's
% Iteration.
% Shimao Fan
% Feb 15 2012
% Temple University
%*********************************************
% this is for 3-parameter arctangent function
% modified May 28 2013
% Proj counts for different projection techniques
%*********************************************
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

tol = 1.0e-02;     % creteria for bisection search

q = u.*rho;   % flux
v = [max(q_min,q),q_max]; % bound for bisection search
f = @(w) q - traffic_flux_fitted_arctan_2p(rho,w,ca,cb,rho0,v0,f0,uf,rhof,rhom);

%\\\\\\\\\\\\\\\\\\\\\\\\
% apply projection here
%////////////////////////

% if data lay above the upper curve
case1 = f(q_min)<=0;   % data below lower bound
case2 = f(q_max)>0;   % data above upper bound
case3 = rho<=rho0;    % data at free regime 
%case4 = rho>rho0;    % data at free regime 

% 
% if case3
% %         w = q_midd;
% %         w = q_min;
%         w = q_max;
% else if case2   % if rho>rho0
%          w = q_max;
%     else if case1
%             w = q_min;
%         else
%            w = bisection(f,v,tol);
%         end
%     end
% end

%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
switch proj
    case 1
      if case3
          w = q_midd;
      else if case2
          w = q_max;
          else if case1
                  w = q_min;
              else
                  w = bisection(f,v,tol);
              end
          end
      end

 %     if case3
 %        w = q_midd;
 %     end
    case 2  
        if case1 
            w = q_min;
        else if case2
            w = q_max;
            else
             w = bisection(f,v,tol);
            end
        end

end
%\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\

%-----------------------------

% if case3    % on the left of rho_f
%     w = q_midd;
% else if case2
%         w = q_max;
%     else if case1
%             w = q_min;         
%     else
%         w = bisection(f,v,tol);
%         end
%     end
% end

%-------------------------------

% 
% if case1
%     w = q_min;
% else if case2
%         w = q_max;
%     else
%         w = bisection(f,v,tol);
%     end
% end
% 
% if case3
%    w = q_midd;
% %    w = q_min;
% %    w = q_max; 
% end

