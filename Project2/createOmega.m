function omega_d = createOmega(Dat,p,T,M,K,vc,lam)
%% The first observation is given and the rest of observation up to the
% end of current window are for prediction as stated in the problem
% statement
% Follow Hamilton 1994 to get the MLE for variance
sum1 = sum(Dat(p:p+T-1,:),1);%sum_1^T y_(t-1)
sum2 = sum(Dat(p+1:p+T,:),1);%sum_1^T y_(t)
sum3 = sum(Dat(p:p+T-1,:).^2,1);%sum_1^T y_(t-1)^2
sum4 = sum(Dat(p+1:p+T,:).^2,1);%sum_1^T y_(t)^2
sum5 = sum(Dat(p+1:p+T,:).*Dat(p:p+T-1,:),1);%sum_1^T y_(t)y_(t-1)
sigma2_hat = zeros(1,M);
for m = 1:M
    c_phi_hat = [T,sum1(m);sum1(m),sum3(m)]\[sum2(m);sum5(m)];% c_phi_hat = [c_hat,phi_hat]
    %temp_sum = 0;
    sigma2_hat(m) = 1.0/(T)*(sum4(m)+c_phi_hat(2)^2*sum3(m)+2*c_phi_hat(1)*c_phi_hat(2)*sum1(m)...
                    -2*c_phi_hat(1)*sum2(m)-2*c_phi_hat(2)*sum5(m))+c_phi_hat(1)^2;
end
omega_d = zeros(1,K);
omega_d(1) = [vc/lam^2];
for j = 1:p
    omega_d(1+(j-1)*M+1:1+j*M) = (sigma2_hat*j^2).^(-1);
end
omega_d = lam^2 * omega_d;
end