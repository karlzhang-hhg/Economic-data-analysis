function imf = IMF(T,M,p,A0,B_hat,mp_shock,time_range)
%% Use companion form to get the impulse response function interms of time.
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
MP_shock = sparse([mp_shock;zeros((p-1)*M,1)]);
% Big matrix pre-multiplied on shocks
G = sparse(zeros(p*M,p*M));
G(1:M,1:M) = A0\eye(M);

imf = zeros(T,M);
temp = G*MP_shock;
for i = 1:T
   temp = Phi*temp;
   imf(i,:) = temp(1:M)';
end
imf = [(time_range(1):(time_range(2)-time_range(1))/(T-1):time_range(2))',imf];
end