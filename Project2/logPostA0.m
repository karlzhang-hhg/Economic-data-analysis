function logfval = logPostA0(a0,X,Y,d_phi,b_hat,B_hat,M,T,lam)
%% Calculate the log posterior for A_0 and return the B_hat
%Input:
%lam: lambda;
%A0: The A_0 matrix (M-by-M);
%X: X matrix ((T+M)-by-K);
%Y: Y matrix ((T+M)-by-M);
%d_phi: Diagonal entries for square root of Phi matrix (K entries)
%b_hat: Matrix form for the mean of prior for beta (which is a vector form of B);
%B_hat: Matrix of coefficient (K-by-M);
%M: Dimension of vector variables;
%T: Time periods for forcasting;
%p: The maximum lag;
%Output:
%logfval: The log posterior function value;

%% ==========================================================================
%MAP for B;
% B_hat = (lam^2*(X'*X)+diag(d_phi.^(-1)))\(lam^2*X'*Y + diag(d_phi.^(-1))*b_hat);
%Residual;
eps_hat = Y - X*B_hat;
%Assemble the A0 from a0
A0 = vec2mat(a0,M);
%A simple way to calculate the arguments inside exponant;
% temp1 = eps_hat*A0';
% temp2 = lam^(-1)*(diag(d_phi.^0.5)\(B_hat-b_hat))*A0';
%The trace of a matrix times its transpose is just the Frobenius norm of
%that matrix;
%logfval = (T+M)*log(det(A0))-0.5*(sum(sum(temp1.^2))+sum(sum(temp2.^2))); %equivalent form
logfval = (T+M)*log(det(A0))-0.5*trace((eps_hat'*eps_hat+lam^(-2)*(B_hat-b_hat)'/diag(d_phi)*(B_hat-b_hat))*(A0'*A0));
end