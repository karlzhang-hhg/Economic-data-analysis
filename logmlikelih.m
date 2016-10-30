function val = logmlikelih(lam,X,Y,d_psi,d_phi,b_hat,M,T,d)
% Calculate the log marginal likelihood
% Input:
% lam: lambda value which we are optimizing in terms of
% X: matrix X in time series (lag);
% Y: matrix Y which are to be forecasted;
% d_psi: diagnoal elements of Psi matrix (which is diagonal);
% d_ome: diagonal elements of Omega matrix (which is diagonal);
% M: the dimension of y
% T: time horizon to be forecasted
% d: some constant specified in prior
% Output:
% val: the log marginal likelihood up to adding a constant
tempm1 = X*diag(d_phi.^0.5);    
eig1 = eig(tempm1'*tempm1);
B_hat = (lam^2*(X'*X)+diag(d_phi.^(-1)))\(lam^2*X'*Y + diag(d_phi.^(-1))*b_hat);
eps_hat = Y - X*B_hat;
tempm2 = B_hat-b_hat;
eig2 = eig(diag(d_psi.^(-0.5))*(lam^2*(eps_hat'*eps_hat)+tempm2'*diag(d_phi.^(-1))*tempm2)*diag(d_psi.^(-0.5)));
val = (-M/2)*sum(log(1+lam^2*eig1))+(-(T+d)/2)*sum(log(lam^2+eig2))+M*(T+d)*log(lam);
end