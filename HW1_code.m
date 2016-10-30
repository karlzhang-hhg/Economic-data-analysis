clear();
% Load data into the working space
load('dataVARmedium','-mat');
p = 5;
T = 64-5;
N = size(y,1);
M = 7;
K = M*p+1;
h = 4;

%% (1) OLS solution for prediction
%% Moving window prediction
tic
z_comp = zeros(N-(T+p+h)+1,8); % the first four columns are predicted z (average growth rate of GDP) and the second four columns is the real z
for i = 1:(N-(T+p+h)+1)
    % design matrix
    X = zeros(T,K);% a T*K matrix
    for j = 1:T
        X(j,:) = [1,reshape(y((i+j+p-2:-1:i+j-1),:)',1,K-1)];% have to transpose because the reshape function operate in column
    end
    % OLS solution for coefficients
    B_hat = (X'*X)\X'*y((i+p:i+p+T-1),:);% a K*M matrix
    % for h = 1 & 4
    pred_y = zeros(h,M);% a h*M matrix
    y_lag = [1,reshape(y(i+p+T-1:-1:i+T,:)',1,K-1)];% lag of y for the prediction
    for k = 1:h
        % forecast
        pred_y(k,:) = y_lag*B_hat;
        if k<h
            % reform the lag for y of prediction
            y_lag = [1,reshape(pred_y(k:-1:1,:)',1,k*M),reshape(y(i+p+T-1:-1:i+T+k,:)',1,K-1-k*M)];
        end
    end
    z_comp(i,1:2) = pred_y(1,1:2)-y(i+p+T-1,1:2);%predicted average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,3:4) = (pred_y(h,1:2)-y(i+p+T-1,1:2))/h;%predicted average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,5:6) = y(i+p+T,1:2)-y(i+p+T-1,1:2);%average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,7:8) = (y(i+p+T+h-1,1:2)-y(i+p+T-1,1:2))/h;%average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
end
% Calculate the average squared forecast error
MSFE1 = sum((z_comp(:,1:4)-z_comp(:,5:8)).^2,1)/size(z_comp,1)
csvwrite('mse_ex1-1-moving.csv',MSFE1)
toc
%% Expanding window prediction
tic
z_comp = zeros(N-(T+p+h)+1,8); % the first four columns are predicted z (average growth rate of GDP) and the second four columns is the real z
for i = 1:(N-(T+p+h)+1)
    if i ==1
        % design matrix
        X = zeros(T,K);% a T*K matrix
        for j = 1:T
            X(j,:) = [1,reshape(y((j+p-1:-1:j),:)',1,K-1)];% have to transpose because the reshape function operate in column
        end
    else
        new_row = [1,reshape(y((i+T+p-2:-1:i+T-1),:)',1,K-1)];
        X = [X;new_row]; % a (T+i-1)*K matrix
    end
    % OLS solution for coefficients
    B_hat = (X'*X)\X'*y((p+1:p+T+i-1),:);% a K*M matrix
    % for h = 1 & 4
    pred_y = zeros(h,M);% a h*M matrix
    y_lag = [1,reshape(y(T+i+p-1:-1:T+i,:)',1,K-1)];% lag of y for the prediction
    for k = 1:h
        % forecast
        pred_y(k,:) = y_lag*B_hat;
        if k<h
            % reform the lag for y of prediction
            y_lag = [1,reshape(pred_y(k:-1:1,:)',1,k*M),reshape(y(T+i+p-1:-1:T+i+k,:)',1,K-1-k*M)];
        end
    end
    z_comp(i,1:2) = pred_y(1,1:2)-y(i+p+T-1,1:2);%predicted average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,3:4) = (pred_y(h,1:2)-y(i+p+T-1,1:2))/h;%predicted average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,5:6) = y(i+p+T,1:2)-y(i+p+T-1,1:2);%average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,7:8) = (y(i+p+T+h-1,1:2)-y(i+p+T-1,1:2))/h;%average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
end
% Calculate the average squared forecast error
MSFE2 = sum((z_comp(:,1:4)-z_comp(:,5:8)).^2,1)/size(z_comp,1)
csvwrite('mse_ex1-1-extending.csv',MSFE2)
toc
%% (2)Minnesota prior with a given lambda
%% Expanding window prediction
tic
z_comp = zeros(N-(T+p+h)+1,8); % the first four columns are predicted z (average growth rate of GDP) and the second four columns is the real z
% Given lambda
lam = 0.2;
% Estimate b
b_hat = [zeros(M,1),eye(M,M*p)];
for i = 1:(N-(T+p+h)+1)
    if i ==1
        % design matrix
        X = zeros(T,K);% a T*K matrix
        for j = 1:T
            X(j,:) = [1,reshape(y((j+p-1:-1:j),:)',1,K-1)];% have to transpose because the reshape function operate in column
        end
    else
        new_row = [1,reshape(y((i+T+p-2:-1:i+T-1),:)',1,K-1)];
        X = [X;new_row]; % a (T+i-1)*K matrix
    end
    % Estimate elements of Omega
    % The first observation is given (regarded as index 0) and the rest of observation up to the
    % end of current window (index 1 to T) are for prediction as stated in the problem
    % statement
    % Follow Hamilton 1994 to get the MLE for variance
    sum1 = sum(y(1:p+T+i-2,:),1);%sum_1^T y_(t-1)
    sum2 = sum(y(2:p+T+i-1,:),1);%sum_1^T y_(t)
    sum3 = sum(y(1:p+T+i-2,:).^2,1);%sum_1^T y_(t-1)^2
    sum4 = sum(y(2:p+T+i-1,:).^2,1);%sum_1^T y_(t)^2
    sum5 = sum(y(2:p+T+i-1,:).*y(1:p+T+i-2,:),1);%sum_1^T y_(t)y_(t-1)
    sigma2_hat = zeros(1,M);
    for m = 1:M
        c_phi_hat = [p+T+i-2,sum1(m);sum1(m),sum3(m)]\[sum2(m);sum5(m)];% c_phi_hat = [c_hat,phi_hat]
        sigma2_hat(m) = 1.0/(p+T+i-2)*(sum4(m)+c_phi_hat(2)^2*sum3(m)+2*c_phi_hat(1)*c_phi_hat(2)*sum1(m)...
                        -2*c_phi_hat(1)*sum2(m)-2*c_phi_hat(2)*sum5(m))+c_phi_hat(1)^2;
    end
    omega_d = [10^6];
    for j = 1:p
        omega_d = [omega_d,(sigma2_hat*j^2).^(-1)];
    end
    % Times the lambda^2
    omega_d = lam^2 * omega_d;
    % OLS solution for coefficients
    B_hat = (X'*X+diag(omega_d.^(-1)))\(X'*y((p+1:p+T+i-1),:)+diag(omega_d.^(-1))*b_hat');% a K*M matrix
    % for h = 1 & 4
    pred_y = zeros(h,M);% a h*M matrix
    y_lag = [1,reshape(y(T+i+p-1:-1:T+i,:)',1,K-1)];% lag of y for the prediction
    for k = 1:h
        % forecast
        pred_y(k,:) = y_lag*B_hat;
        if k<h
            % reform the lag for y of prediction
            y_lag = [1,reshape(pred_y(k:-1:1,:)',1,k*M),reshape(y(T+i+p-1:-1:T+i+k,:)',1,K-1-k*M)];
        end
    end
    z_comp(i,1:2) = pred_y(1,1:2)-y(i+p+T-1,1:2);%predicted average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,3:4) = (pred_y(h,1:2)-y(i+p+T-1,1:2))/h;%predicted average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,5:6) = y(i+p+T,1:2)-y(i+p+T-1,1:2);%average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,7:8) = (y(i+p+T+h-1,1:2)-y(i+p+T-1,1:2))/h;%average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
end
% Calculate the average squared forecast error
MSFE3 = sum((z_comp(:,1:4)-z_comp(:,5:8)).^2,1)/size(z_comp,1)
csvwrite('mse_ex1-1-lambda_given.csv',MSFE3)
toc
%% (3)Minnesota prior without a given lambda (optimize marginal likelihood to get lambda each time)
%% Expanding window prediction
tic
%Cannot set lambda = 0, because when take log of lambda, 0 is not
%meaningful
lam_bd = [0.0001,1];
N_try = 1000;
d = M + 2;
z_comp = zeros(N-(T+p+h)+1,8); % the first four columns are predicted z (average growth rate of GDP) and the second four columns is the real z
% Record optimal value of lambda
rec_lam_opt = zeros(N-(T+p+h)+1,1);
% Estimate b
b_hat = [zeros(M,1),eye(M,M*p)];
for i = 1:(N-(T+p+h)+1)
    if i ==1
        % design matrix
        X = zeros(T,K);% a T*K matrix
        for j = 1:T
            X(j,:) = [1,reshape(y((j+p-1:-1:j),:)',1,K-1)];% have to transpose because the reshape function operate in column
        end
    else
        new_row = [1,reshape(y((i+T+p-2:-1:i+T-1),:)',1,K-1)];
        X = [X;new_row]; % a (T+i-1)*K matrix
    end
    % Estimate elements of Omega
    % The first observation is given and the rest of observation up to the
    % end of current window are for prediction as stated in the problem
    % statement
    % Follow Hamilton 1994 to get the MLE for variance
    sum1 = sum(y(1:p+T+i-2,:),1);%sum_1^T y_(t-1)
    sum2 = sum(y(2:p+T+i-1,:),1);%sum_1^T y_(t)
    sum3 = sum(y(1:p+T+i-2,:).^2,1);%sum_1^T y_(t-1)^2
    sum4 = sum(y(2:p+T+i-1,:).^2,1);%sum_1^T y_(t)^2
    sum5 = sum(y(2:p+T+i-1,:).*y(1:p+T+i-2,:),1);%sum_1^T y_(t)y_(t-1)
    sigma2_hat = zeros(1,M);
    for m = 1:M
        c_phi_hat = [p+T+i-2,sum1(m);sum1(m),sum3(m)]\[sum2(m);sum5(m)];% c_phi_hat = [c_hat,phi_hat]
        temp_sum = 0;
        sigma2_hat(m) = 1.0/(p+T+i-2)*(sum4(m)+c_phi_hat(2)^2*sum3(m)+2*c_phi_hat(1)*c_phi_hat(2)*sum1(m)...
                        -2*c_phi_hat(1)*sum2(m)-2*c_phi_hat(2)*sum5(m))+c_phi_hat(1)^2;
    end
    omega_d = [10^6];
    for j = 1:p
        omega_d = [omega_d,(sigma2_hat*j^2).^(-1)];
    end
    % Optimize in terms of lambda
    % Use grid-search method
    %options_optim = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','quasi-newton','SpecifyObjectiveGradient',true,'MaxIterations',1000);
    rec_lam_opt(i) = lam_mle_grid(@logmlikelih,lam_bd,N_try,X,y((p+1:p+T+i-1),:),sigma2_hat,omega_d,b_hat',M,T+i-1,d);
    % Times the lambda^2
    omega_d = rec_lam_opt(i)^2 * omega_d;
    % OLS solution for coefficients
    B_hat = (X'*X+diag(omega_d.^(-1)))\(X'*y((p+1:p+T+i-1),:)+diag(omega_d.^(-1))*b_hat');% a K*M matrix
    % for h = 1 & 4
    pred_y = zeros(h,M);% a h*M matrix
    y_lag = [1,reshape(y(T+i+p-1:-1:T+i,:)',1,K-1)];% lag of y for the prediction
    for k = 1:h
        % forecast
        pred_y(k,:) = y_lag*B_hat;
        if k<h
            % reform the lag for y of prediction
            y_lag = [1,reshape(pred_y(k:-1:1,:)',1,k*M),reshape(y(T+i+p-1:-1:T+i+k,:)',1,K-1-k*M)];
        end
    end
    z_comp(i,1:2) = pred_y(1,1:2)-y(i+p+T-1,1:2);%predicted average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,3:4) = (pred_y(h,1:2)-y(i+p+T-1,1:2))/h;%predicted average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,5:6) = y(i+p+T,1:2)-y(i+p+T-1,1:2);%average growth rates for one quarter of (i) log-real GDP and (ii) log-GDP delector
    z_comp(i,7:8) = (y(i+p+T+h-1,1:2)-y(i+p+T-1,1:2))/h;%average growth rates for h quarters of (i) log-real GDP and (ii) log-GDP delector
end
% Calculate the average squared forecast error
MSFE4 = sum((z_comp(:,1:4)-z_comp(:,5:8)).^2,1)/size(z_comp,1)
csvwrite('mse_ex1-1-lambda_opt.csv',MSFE4)
csvwrite('lambda_opt.csv',rec_lam_opt)
figure();
hold on;
plot(rec_lam_opt);
int_str = ['[',num2str(lam_bd(1)),',',num2str(lam_bd(2)),']'];
title(['Grid search interval:',int_str]);
cwd = '/Users/kungangzhang/Documents/OneDrive/Northwestern/Study/Courses/ECON-482/HW1/';
saveas(gca,[cwd,['lambda_opt_1']],'fig');
saveas(gca,[cwd,['lambda_opt_1']],'jpg');
toc