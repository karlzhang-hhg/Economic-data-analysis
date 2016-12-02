function rt = hw1ex2bforecast(data,Tulim,p) % Tulim in no. of quarters

T = size(data,1);
n = size(data,2);

d = n+2;
lambda = 0.2;

%calculate the variance of residuals from AR(1) model, which is sigma^2 hat

xar = data(1:(Tulim-1),:);
yar = data(2:Tulim,:);

sigma = zeros(1,n);
for i = 1:1:n
    mdl = fitlm(xar(:,i),yar(:,i));
    r = mdl.Residuals.Raw; %resid
    sigma(i) = std(r)^2;
end;

% set up b and omega in Minnesota prior and other parameters

btemp = zeros(n, n*p+1);
btemp(1:n,2:(n+1)) = eye(n);
btilde = btemp';

% b = reshape(btilde, n*(n*p+1),1);
omega = zeros((n*p+1),(n*p+1));
for s = 1:1:p
    omega((1+n*(s-1)+1):(1+n*s),(1+n*(s-1)+1):(1+n*s))= ...
        diag(1./(sigma*(s^2)));
end;
omega(1,1) = 10^3; % chose a constant randomly
omega = lambda^2 * omega;

% psi = diag(sigma);

% set up y and x matrices

y = data(1:Tulim,:);
ytrim = data((p+1):Tulim,:);
Tlen = size(ytrim,1);

Y = reshape(ytrim,Tlen*n,1);

x = ones(Tlen,1);
for i = 1:1:p
    temp = y((p+1-i):(end-i), :);
    x = [x, temp];
end
In = eye(n);
X = kron(In,x);
% posterior distributions

% beta posterior
Bhat = inv(x'*x+ inv(omega))*...
    (x'*ytrim + inv(omega)*btilde);
betapos = reshape(Bhat,size(Bhat,1)*size(Bhat,2),1);

% predictions
h = Tulim-p-1;
ypred = y(1:p,:);
ypred(p+1,:) = x(1,:)* Bhat;

covmat = ones(h,1+n*p);
for i = 1:1:h
    for j = 1:1:p
        covmat(i,(1+n*(j-1)+1):(1+(n*j))) = ypred(end+1-j,:); 
    end;
    ypred= [ypred;covmat(i,:)*Bhat];
end;


% ypred = [y;zeros(h,n)];
% covmat = zeros(n,n*(n*p+1));
% for k = 1:1:h
%     for i = 1:1:n
%         covmat(i,(i-1)*(n*p+1)+1) = 1;
%         for j = 1:1:p
%         covmat(i,((i-1)*(n*p+1)+(1+n*(j-1)+1)):((i-1)*(n*p+1)+(1+(n*j)))) ...
%             = ypred(Tulim+k-j,:);        
%         end;
%     end;
%     ypred(Tulim+k,:) = covmat*betapos;
% end;

% predictions

rt.predgdp = ypred(:,1); % log real GDP
rt.predp = diff(ypred(:,2));   % price inflation
rt.predff = ypred(:,3); % federal funds rate

% ngdp = exp(col2)*exp(col1)/100
rt.prednwi = diff(ypred(:,7)+ypred(:,2)); % nominal wage inflation
rt.predls = ypred(:,7)+ ypred(:,6)- ypred(:,1); 
    % log labor share

rt.predlcr = ypred(:,4) - ypred(:,1);

% observations

rt.realgdp = data(:,1); % log real GDP
rt.realp = diff(data(:,2));   % price inflation
rt.realff = data(:,3); % federal funds rate
    % ngdp = exp(col2)*exp(col1)/100
rt.realnwi = diff(data(:,7)+data(:,2)); % nominal wage inflation

rt.realls = data(:,7)+ data(:,6)- data(:,1); 
        % log labor share
rt.reallcr =  data(:,4) - data(:,1);


% regression data~ ypred (dc) parameter by parameter
rsq = zeros(1,6);
rsq(1,1) = dcrsq(rt.realgdp, rt.predgdp);
rsq(1,2) = dcrsq(rt.realp, rt.predp);
rsq(1,3) = dcrsq(rt.realff, rt.predff);
rsq(1,4) = dcrsq(rt.realnwi, rt.prednwi);
rsq(1,5) = dcrsq(rt.realls, rt.predls);
rsq(1,6) = dcrsq(rt.reallcr, rt.predlcr);

rt.rsq =rsq;





