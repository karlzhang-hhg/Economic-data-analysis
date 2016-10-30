% Regression Function  
function r = regsolv(i,d)
         cons = ones(size(d,1),1);
         x    = [cons i];
         r.b    = inv(x'*x)*x'*d;
         r.ehat = d - x*r.b;
         r.that = d - mean(d);
         r.SST  = r.that'*r.that;
         r.SSR  = r.ehat'*r.ehat;
         r.Rsq= 1-r.SSR/r.SST;
end