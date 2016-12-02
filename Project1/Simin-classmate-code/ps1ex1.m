%% ~~~~~~~~~~~~~~~~~~~~ Exercise 1 ~~~~~~~~~~~~~~~~~~~~ %%
clear all
load dataVARmedium.mat
yfull = y;

% Set the prior (1=flat, 2=mn, 4= mn with opt lambda)
prior = 4;

% Feed into var function
p = 5;
n = size(y,2);
h = 4; % max per ahead forecast
g = 1; % other per ahead forecast
sy = 1959;
my = 1974;
ey = 2008;
freqnum = 4;

initT = (my-sy+1)*freqnum;
TT    = size(yfull,1);

% Variables for the forecast
ypred1 = zeros(TT+1,n);
ypred4 = zeros(TT+4,n);

for j = 0:(ey-my)*4 
   T = initT+j-p;
   y = yfull(1:(initT+j),:); 
   v = varsolv(p, n, y, h, g, prior, ey, sy, my, freqnum);
   % Collect the forecasts
   ypred1(initT+j+1,:) = v.f1;
   ypred4(initT+j+4,:) = v.f4;
   % Collect optimal lambdas
   if prior == 4
    lambda(initT+j+1,1) = v.lambda_opt;
   end
end

% Compute forecasted growth rate vs. actual
zhat1 = 1/g *(ypred1(initT+g :end-g,:)  - y(initT:end-g,:));
zhat4 = 1/h *(ypred4(initT+h :end-h,:)  - y(initT:end-h,:));

z1    = 1/g * diff(y(initT:end,:));
z4    = 1/h * (y(initT+h:end,:) - y(initT:end-h,:));

MSFE1 = mean((zhat1-z1).^2);
MSFE4 = mean((zhat4-z4).^2);

MSFE  = [MSFE1(1:2) MSFE4(1:2)];
% d1 = digits(3);
% tabMSFE = latex(vpa(sym(MSFE)));
% digits(d1);

% Plot the optimal lambda for the MN prior
if prior==4
    A = sy+1/freqnum:1/freqnum:(ey+1);
    t = A(initT+1:end)';
   
    figure
    plot(t,lambda(initT+1:end-1))
    title('Optimal Lambda for Minnesota Prior')
    saveas(gcf,'ps1ex1.jpg')
end
