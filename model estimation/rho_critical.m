function y = rho_critical(lambda,p,rhom)


% find critical density of smooth flux function or FDs

      a = sqrt(1+(p.*lambda).^2);
      b = sqrt(1+((1-p).*lambda).^2);
      y = rhom*((b-a)./(lambda.*(sqrt(lambda.^2-(b-a).^2)))+p);