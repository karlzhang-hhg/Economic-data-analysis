function a0_0 = initPoint(X,Y,B_hat,b_hat,omega_d)
%Residual;
eps_hat = Y - X*B_hat;
%Another part;
temp1 = (diag(omega_d.^0.5)\(B_hat-b_hat));
%S_hat
S_hat = eps_hat'*eps_hat+temp1'*temp1;
%Cholesky decomposition on S_hat'*S_hat to get the initial A0 and then
%return the vector version of it.
A0_0 = chol(S_hat'*S_hat,'upper');
a0_0 = mat2vec(A0_0);
end