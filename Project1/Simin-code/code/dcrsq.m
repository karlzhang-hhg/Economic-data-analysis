function rt = dcrsq(real, pred)
x = [ones(size(pred,1),1) pred];

betahat = inv(x'*x)*x'*real;
ssr = sum((real - x*betahat).^2);
sst =  sum((real- mean(real)).^2);
rsquare = 1-ssr/sst; 

rt = rsquare;