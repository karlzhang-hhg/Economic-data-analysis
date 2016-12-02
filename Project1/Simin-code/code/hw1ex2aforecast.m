function rt = hw1ex2aforecast(data,Tulim,p) % Tulim in no. of quarters

T = size(data,1);
n = size(data,2);

% condition on the first 5 periods. starting from the 6th row
y = data(1:Tulim,:);
ytrim = data((p+1):Tulim,:);
Tlen = size(ytrim,1);

Y = reshape(ytrim,Tlen*n,1);

x = ones(Tlen,1);
for i = 1:1:p
    temp = y((p+1-i):(end-i), :);
    x = [x, temp];
end
%x is a Tlen*(np+1) matrix
In = eye(n);

X = kron(In,x);

% mdl = LinearModel.fit(X,Y);
% mdl = fitlm(X,Y,'Intercept',false);
% betaols = mdl.Coefficients.Estimate; % return coeff


betaols = inv(x'*x)*(x'*ytrim);% betaols is a (np+1)*n matrix


% predictions
h = Tulim-p-1;
ypred = y(1:p,:);
ypred(p+1,:) = x(1,:)* betaols;


covmat = ones(h,1+n*p);
for i = 1:1:h
    for j = 1:1:p
        covmat(i,(1+n*(j-1)+1):(1+(n*j))) = ypred(end+1-j,:); 
    end;
    ypred= [ypred;covmat(i,:)*betaols];
end;

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
rt.rsq = rsq;


