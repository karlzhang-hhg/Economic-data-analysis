function var_imf = imf_var(M,p,A0,B_hat,shock,t_proj)
%% Use companion form to get the variance of impulse response function with certain length of projection.
%Input:
%   T: The time horizon to investigate the impulse response;
%   M: Dimension of vector variable;
%   p: The maximum time lag;
%   A0: The MAP of A0 matrix (M-by-M);
%   B_hat: The MAP of B matrix (K-by-M);
%   mp_shock: To specify which shock is in interest (M-by-1);
%Output:
%   imf: Impulse response function in terms of time periods (T-by-M);

% Phi is a K-by-K matrix
Phi = sparse([B_hat(2:end,:)';eye((p-1)*M,p*M)]);
% Big shock column: K-by-1
Shock = sparse([shock;zeros((p-1)*M,1)]);
% Big matrix pre-multiplied on shocks
G = sparse(zeros(p*M,p*M));
G(1:M,1:M) = A0\eye(M);

temp = G*diag(Shock)*G';
var_imf = temp;
for i = 1:t_proj-1
   temp = Phi*temp*Phi';
   var_imf = var_imf + temp;
end
var_imf = var_imf(1:M,1:M);
end