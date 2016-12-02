function v = varsolv(p, n, y, h, g, prior, ey, sy, my, freqnum)
%% Runs a VAR(p) with n variables using data y, 
% forecasts with max horizon h periods and g periods (1), 
% under possible priors for sy = startyear, ey=endyear, 
% my = end of initial sample year at frequency freqnum
% Prior 1: Flat (ols)
% Prior 2: MN prior
% Prior 3: SOC

k = 1+n*p;
T = size(y,1)-p;
initT = size(y,1);

% Form the y and x matrices for regression
yreg = y(p+1:end,:);

x = ones(T,1);
for i = 1:p
    temp = y((p+1-i):(end-i), :);
    x = [x, temp];
end

% Define variables for plots
    gdp      = y(:,1);
    infl     = diff(y(:,2));
    ff       = y(:,3);
    nwage    = diff(y(:,7)+ y(:,2));
    lshare   = y(:,7)+y(:,6)-y(:,1);
    cshare   = y(:,4)-y(:,1);  

%% Prior 1: Get coefficients B
if prior == 1
    v.B = inv(x'*x) * x'*yreg;
    
   % Forecast GDP and GDP deflator one and four periods ahead (ex 1)
   xhat = ones(h, k);
   for i = 1:h
      for l = 1:p
         sel = (1 + (l-1)*n) : (l*n);
         xhat(i, 1+sel) = y(end+1-l, :);
      end
      yhat = xhat(i, :) * v.B;
      y = [y; yhat];
   end
   v.f1 = y(end-h+1,:);
   v.f4 = y(end,:);
   
   % Forecast the DC (ex 2)
    v.ydc        = y(1:p,:);
    xdc        = x(1,:);
    ydcf       = xdc(1,:) * v.B; 
    v.ydc(p+1,:) = ydcf;
    xhat       = ones(T,k);
    for t = 1:T-1
        for l = 1:p
            sel = (1 + (l-1)*n) : (l*n);
            xhat(t, 1+sel) = v.ydc(end+1-l, :);
        end
        yhat = xhat(t, :) * v.B;
        v.ydc = [v.ydc; yhat];
    end
    
    % Define variables for plots and regressions
    v.gdp  = v.ydc(:,1); % GDP
    v.infl = diff(v.ydc(:,2)); % infl
    v.ff   = v.ydc(:,3); % FF
    v.nwage  = diff(v.ydc(:,7)+ v.ydc(:,2)); % nom wage infl
    v.lshare = v.ydc(:,7)+v.ydc(:,6)-v.ydc(:,1); % labor share
    v.cshare = v.ydc(:,4)-v.ydc(:,1); % cons share
    
    % Regressions (ex 2)
    % Initialize vector of R-sq then compute and fill
    v.Rsq = zeros(1,6)
    % GDP
    i       = v.gdp;
    d       = gdp;
    dc         = regsolv(i,d);
    v.Rsq(1)   = dc.Rsq;
    % Infl
    i       = v.infl;
    d       = infl;
    dc         = regsolv(i,d);
    v.Rsq(2)   = dc.Rsq;
    % FF
    i       = v.ff;
    d       = ff;
    dc         = regsolv(i,d);
    v.Rsq(3)   = dc.Rsq;
    % Nwage
    i       = v.nwage;
    d       = nwage;
    dc         = regsolv(i,d);
    v.Rsq(4)   = dc.Rsq;
    % Lshare
    i       = v.lshare;
    d       = lshare;
    dc         = regsolv(i,d);
    v.Rsq(5)   = dc.Rsq;
    % Cshare
    i       = v.cshare;
    d       = cshare;
    dc         = regsolv(i,d);
    v.Rsq(6)   = dc.Rsq;
   

%% Prior 2: MN Prior: 
% beta|sigma ~ N(b, sigma kron omega)
% sigma      ~ IW(psi, d)
% psi (and thus sigma) determined within the loop
% omega is also formed within the loop using sigmahat
% centered on a random walk so b is:
elseif prior == 2
    d = n+2;
    lambda = 0.2;
    alfa   = 1;
    
    % Form the prior mean on beta
    B = zeros(k,n);
    B(2:n+1, :) = eye(n);

    omega  = zeros(k);
    sighat = zeros(n,1);
    omega(1,1) = lambda^2*10^6;
    
    % Get sigmahat, psi, omega
    for i = 1:n
        ary = y(2:end,i);
        arx = [ones(size(y,1)-1,1) y(1:end-1,i)];
        arb =  inv(arx'*arx)*(arx'*ary);
        sighat(i,1) = var(ary - arx*arb); 
    end    
    psi  = diag(sighat);
    isighat = 1./sighat;
    ipsi = diag(isighat);
    sigma  = iwishrnd(psi,n+2);
    for l = 1:p
        omega((l-1)*n+2:l*n+1,(l-1)*n+2:l*n+1) = lambda^2*ipsi*l^(-2); 
    end
    
    % Posterior mode
    v.B = (x'*x + inv(omega)) \ (x'*yreg+inv(omega)*B);
    
    % Forecast GDP and GDP deflator one and four periods ahead (ex 1)
    xhat = ones(h, k);
    for i = 1:h
      for l = 1:p
         sel = (1 + (l-1)*n) : (l*n);
         xhat(i, 1+sel) = y(end+1-l, :);
      end
      yhat = xhat(i, :) * v.B;
      y = [y; yhat];
    end
    v.f1 = y(end-h+1,:);
    v.f4 = y(end,:);
    
        % Forecast the DC (ex 2)
    v.ydc        = y(1:p,:);
    xdc        = x(1,:);
    ydcf       = xdc(1,:) * v.B; 
    v.ydc(p+1,:) = ydcf;
    xhat       = ones(T,k);
    for t = 1:T-1
        for l = 1:p
            sel = (1 + (l-1)*n) : (l*n);
            xhat(t, 1+sel) = v.ydc(end+1-l, :);
        end
        yhat = xhat(t, :) * v.B;
        v.ydc = [v.ydc; yhat];
    end
    
        % Define variables for plots and regressions
    v.gdp  = v.ydc(:,1); % GDP
    v.infl = diff(v.ydc(:,2)); % infl
    v.ff   = v.ydc(:,3); % FF
    v.nwage  = diff(v.ydc(:,7)+ v.ydc(:,2)); % nom wage infl
    v.lshare = v.ydc(:,7)+v.ydc(:,6)-v.ydc(:,1); % labor share
    v.cshare = v.ydc(:,4)-v.ydc(:,1); % cons share
    
    % Regressions (ex 2)
    % Initialize vector of R-sq then compute and fill
    v.Rsq = zeros(1,6)
    % GDP
    i       = v.gdp;
    d       = gdp;
    dc         = regsolv(i,d);
    v.Rsq(1)   = dc.Rsq;
    % Infl
    i       = v.infl;
    d       = infl;
    dc         = regsolv(i,d);
    v.Rsq(2)   = dc.Rsq;
    % FF
    i       = v.ff;
    d       = ff;
    dc         = regsolv(i,d);
    v.Rsq(3)   = dc.Rsq;
    % Nwage
    i       = v.nwage;
    d       = nwage;
    dc         = regsolv(i,d);
    v.Rsq(4)   = dc.Rsq;
    % Lshare
    i       = v.lshare;
    d       = lshare;
    dc         = regsolv(i,d);
    v.Rsq(5)   = dc.Rsq;
    % Cshare
    i       = v.cshare;
    d       = cshare;
    dc         = regsolv(i,d);
    v.Rsq(6)   = dc.Rsq;
    
%% Prior 3: SOC
% Add n rows with mean of prior obs
elseif prior == 3
    mu    = 1;
    ystar = diag(mean(y(1:p,:)))*mu;
    xstar = [zeros(n,1) repmat(ystar,[1 p])];
    
    y = [yreg;ystar];
    x = [x;xstar];
    
    v.B = inv(x'*x) * x'*y;
    
    % Forecast the DC (ex 2)
    v.ydc        = y(1:p,:);
    xdc        = x(1,:);
    ydcf       = xdc(1,:) * v.B; 
    v.ydc(p+1,:) = ydcf;
    xhat       = ones(T,k);
    for t = 1:T-1
        for l = 1:p
            sel = (1 + (l-1)*n) : (l*n);
            xhat(t, 1+sel) = v.ydc(end+1-l, :);
        end
        yhat = xhat(t, :) * v.B;
        v.ydc = [v.ydc; yhat];
    end
    
        % Define variables for plots and regressions
    v.gdp  = v.ydc(:,1); % GDP
    v.infl = diff(v.ydc(:,2)); % infl
    v.ff   = v.ydc(:,3); % FF
    v.nwage  = diff(v.ydc(:,7)+ v.ydc(:,2)); % nom wage infl
    v.lshare = v.ydc(:,7)+v.ydc(:,6)-v.ydc(:,1); % labor share
    v.cshare = v.ydc(:,4)-v.ydc(:,1); % cons share
    
    % Regressions (ex 2)
    % Initialize vector of R-sq then compute and fill
    v.Rsq = zeros(1,6)
    % GDP
    i       = v.gdp;
    d       = gdp;
    dc         = regsolv(i,d);
    v.Rsq(1)   = dc.Rsq;
    % Infl
    i       = v.infl;
    d       = infl;
    dc         = regsolv(i,d);
    v.Rsq(2)   = dc.Rsq;
    % FF
    i       = v.ff;
    d       = ff;
    dc         = regsolv(i,d);
    v.Rsq(3)   = dc.Rsq;
    % Nwage
    i       = v.nwage;
    d       = nwage;
    dc         = regsolv(i,d);
    v.Rsq(4)   = dc.Rsq;
    % Lshare
    i       = v.lshare;
    d       = lshare;
    dc         = regsolv(i,d);
    v.Rsq(5)   = dc.Rsq;
    % Cshare
    i       = v.cshare;
    d       = cshare;
    dc         = regsolv(i,d);
    v.Rsq(6)   = dc.Rsq;
    
%% Prior 4: MN Prior with optimal lambda (ex1):
% beta|sigma ~ N(b, sigma kron omega)
% sigma      ~ IW(psi, d)
% lambda is chosen to max ML with flat hyperprior
% psi (and thus sigma) determined within the loop
% omega is also formed within the loop using sigmahat
% centered on a random walk so b is:
elseif prior == 4
    d = n+2;
    % lambda = 0.2;
    alfa   = 1;
    
    % Form the prior mean on beta
    btilde = zeros(k,n);
    btilde(2:n+1, :) = eye(n);

    omega  = zeros(k);
    sighat = zeros(n,1);
    omega(1,1) = 10^6;
    
    % Get sigmahat, psi, omega
    for i = 1:n
        ary = y(2:end,i);
        arx = [ones(size(y,1)-1,1) y(1:end-1,i)];
        arb =  inv(arx'*arx)*(arx'*ary);
        sighat(i,1) = var(ary - arx*arb); 
    end    
    psi  = diag(sighat);
    isighat = 1./sighat;
    ipsi = diag(isighat);
    sigma  = iwishrnd(psi,d);
    for l = 1:p
        omega((l-1)*n+2:l*n+1,(l-1)*n+2:l*n+1) = ipsi*l^(-2); 
    end
    
    % Maximize log marginal likelihood wrt lambda to get opt
    ml = @(lam)(marglike(lam,n,T,d,omega,psi,x,yreg,btilde) );
    v.lambda_opt = fminbnd(ml,0.0001,0.2);
    omega_opt = omega * v.lambda_opt^2;
    
    % Posterior mode with opt omega
    v.B = (x'*x + inv(omega_opt)) \ (x'*yreg+inv(omega_opt)*btilde);
    
    % Forecast GDP and GDP deflator one and four periods ahead (ex 1)
    xhat = ones(h, k);
    for i = 1:h
      for l = 1:p
         sel = (1 + (l-1)*n) : (l*n);
         xhat(i, 1+sel) = y(end+1-l, :);
      end
      yhat = xhat(i, :) * v.B;
      y = [y; yhat];
    end
    v.f1 = y(end-h+1,:);
    v.f4 = y(end,:);
  
end

end


