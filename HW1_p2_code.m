clear();
% Load data into the working space
load('dataVARmedium','-mat');

lgGDP = y(:,1);%log-real GDP
qpriinf = diff(y(:,2));%quarterly price inflation
fedfr = y(:,3);%federal funds rate
qnomwainf = diff(y(:,7)+y(:,2));%quarterly nominal wage inflation
lglabsh = y(:,7)+y(:,6)-y(:,1);%log-labor share
lgcomra = y(:,4)-y(:,1);%log-consumption ratio

p = 5;
N = size(y,1);
T = N-p;%Initial time horizon for forecasting
M = 7;
K = M*p+1;
h = 0;%The maximum horizon of out-of-sample forecast

%% (1) OLS solution for prediction
tic
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
    % Do projection based on MLE of coefficients
    yproj1 = y(1:p,:);
    for i = p+1:N
        yproj1(i,:) = [1,reshape(yproj1(i-1:-1:i-p,:)',1,K-1)]*B_hat; 
    end
    lgGDP1 = yproj1(:,1);
    qpriinf1 = diff(yproj1(:,2));
    fedfr1 = yproj1(:,3);
    qnomwainf1 = diff(yproj1(:,7)+yproj1(:,2));
    lglabsh1 = yproj1(:,7)+yproj1(:,6)-yproj1(:,1);
    lgcomra1 = yproj1(:,4)-yproj1(:,1);
end
xt = ((1959+0.25): 0.25: (2008+1))';
figure();
hold on;
box on;
%title('Prejected v.s. real data');

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Prejected v.s. real data for different priors', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

legendinfo{1} = ['Observed data'];
handle(1) = plot_all([xt,lgGDP],[xt(2:end),qpriinf],[xt,fedfr],[xt(2:end),qnomwainf],[xt,lglabsh],[xt,lgcomra],'-b','Observed data');
xt1 = xt(1:end);
legendinfo{2} = ['Flat prior projection'];
handle(2) = plot_all([xt1,lgGDP1],[xt1(2:end),qpriinf1],[xt1,fedfr1],[xt1(2:end),qnomwainf1],[xt1,lglabsh1],[xt1,lgcomra1],'-.r','Flat prior projection');
toc

%% (2)Minnesota prior with a given lambda
tic
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
    % MAP solution for coefficients
    B_hat = (X'*X+diag(omega_d.^(-1)))\(X'*y((p+1:p+T+i-1),:)+diag(omega_d.^(-1))*b_hat');% a K*M matrix
    % Do projection based on MLE of coefficients
    yproj2 = y(1:p,:);
    for i = p+1:N
        yproj2(i,:) = [1,reshape(yproj2(i-1:-1:i-p,:)',1,K-1)]*B_hat; 
    end
    lgGDP2 = yproj2(:,1);
    qpriinf2 = diff(yproj2(:,2));
    fedfr2 = yproj2(:,3);
    qnomwainf2 = diff(yproj2(:,7)+yproj2(:,2));
    lglabsh2 = yproj2(:,7)+yproj2(:,6)-yproj2(:,1);
    lgcomra2 = yproj2(:,4)-yproj2(:,1);
end
xt2 = xt(1:end);
legendinfo{3} = ['Minnesota prior projection (\lambda = 0.2)'];
handle(3) = plot_all([xt2,lgGDP2],[xt2(2:end),qpriinf2],[xt2,fedfr2],[xt2(2:end),qnomwainf2],[xt2,lglabsh2],[xt2,lgcomra2],'-.g','Minnesota prior projection (\lambda = 0.2)');

toc
%% (3)Minnesota prior without a given lambda (optimize marginal likelihood to get lambda each time)
tic
% Given lambda
lam = 0.2;
mu = 1;
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
    
    %stack the artificial observations at the bottom to get soc prior
    Xplus = [zeros(M,1),repmat(diag(mu*mean(y(1:p,:),1)),[1 p])];    
    X= [X;Xplus]; % stack x at the bottom
    yplus = diag(mu*mean(y(1:p,:),1));
    ynew = [y((p+1:p+T+i-1),:);yplus];
    
    % MAP solution for coefficients
    B_hat = (X'*X+diag(omega_d.^(-1)))\(X'*ynew+diag(omega_d.^(-1))*b_hat');% a K*M matrix
    % Do projection based on MLE of coefficients
    yproj3 = y(1:p,:);
    for j = p+1:N
        yproj3(j,:) = [1,reshape(yproj3(j-1:-1:j-p,:)',1,K-1)]*B_hat; 
    end
    lgGDP3 = yproj3(:,1);
    qpriinf3 = diff(yproj3(:,2));
    fedfr3 = yproj3(:,3);
    qnomwainf3 = diff(yproj3(:,7)+yproj3(:,2));
    lglabsh3 = yproj3(:,7)+yproj3(:,6)-yproj3(:,1);
    lgcomra3 = yproj3(:,4)-yproj3(:,1);
end
xt3 = xt(1:end);
legendinfo{4} = ['SOC prior (\mu = 1)'];
handle(4) = plot_all([xt3,lgGDP3],[xt3(2:end),qpriinf3],[xt3,fedfr3],[xt3(2:end),qnomwainf3],[xt3,lglabsh3],[xt3,lgcomra3],'-.m','SOC prior (\mu = 1)');

toc

tic
% Given lambda
lam = 0.2;
mu = 5;
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
    
    %stack the artificial observations at the bottom to get soc prior
    Xplus = [zeros(M,1),repmat(diag(mu*mean(y(1:p,:),1)),[1 p])];    
    X= [X;Xplus]; % stack x at the bottom
    yplus = diag(mu*mean(y(1:p,:),1));
    ynew = [y((p+1:p+T+i-1),:);yplus];
    
    % MAP solution for coefficients
    B_hat = (X'*X+diag(omega_d.^(-1)))\(X'*ynew+diag(omega_d.^(-1))*b_hat');% a K*M matrix
    % Do projection based on MLE of coefficients
    yproj4 = y(1:p,:);
    for j = p+1:N
        yproj4(j,:) = [1,reshape(yproj4(j-1:-1:j-p,:)',1,K-1)]*B_hat; 
    end
    lgGDP4 = yproj4(:,1);
    qpriinf4 = diff(yproj4(:,2));
    fedfr4 = yproj4(:,3);
    qnomwainf4 = diff(yproj4(:,7)+yproj4(:,2));
    lglabsh4 = yproj4(:,7)+yproj4(:,6)-yproj4(:,1);
    lgcomra4 = yproj4(:,4)-yproj4(:,1);
end
xt4 = xt(1:end);
legendinfo{5} = ['SOC prior (\mu = 5)'];
handle(5) = plot_all([xt4,lgGDP4],[xt4(2:end),qpriinf4],[xt4,fedfr4],[xt4(2:end),qnomwainf4],[xt4,lglabsh4],[xt4,lgcomra4],'-.c','SOC prior (\mu = 5)');

toc

var1=[var(lgGDP),var(lgGDP1),var(lgGDP2),var(lgGDP3),var(lgGDP4)];
var2=[var(qpriinf),var(qpriinf1),var(qpriinf2),var(qpriinf3),var(qpriinf4)];
var3=[var(fedfr),var(fedfr1),var(fedfr2),var(fedfr3),var(fedfr4)];
var4=[var(qnomwainf),var(qnomwainf1),var(qnomwainf2),var(qnomwainf3),var(qnomwainf4)];
var5=[var(lglabsh),var(lglabsh1),var(lglabsh2),var(lglabsh3),var(lglabsh4)];
var6=[var(lgcomra),var(lgcomra1),var(lgcomra2),var(lgcomra3),var(lgcomra4)];
var_all = [var1;var2;var3;var4;var5;var6];
legend(handle,legendinfo,'Location','southeast');