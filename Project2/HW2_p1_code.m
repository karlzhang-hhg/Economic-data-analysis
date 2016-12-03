clear();
%Load the data set
Dat = xlsread('SZdata.xlsx');
Dat = Dat(:,2:7);
%Set constants
p = 13;
lam = 0.2;
M = 6;
T = size(Dat,1)-p;
K = M*p+1;
vc = 10^6;
mu = 1;
Num = 500000;%Number of MCMC simulation
C = 0.5;%Scaling parameter for the proposal distribution variation
ir_t = 48;
mp_shock = zeros(M,1);
mp_shock(6) = 1;
multistarts = 2500;
t_proj = 36;
%Create matrix of X and Y combining observations and dummy observations
[X,Y] = createXY(Dat,p,T,M,K,mu);

%Estimate elements of Omega
omega_d = createOmega(Dat,p,T,M,K,vc,lam);

%Estimate the MAP for B and constant for b_hat in prior
b_hat = [zeros(M,1),eye(M,M*p)]';
B_hat = ((X'*X)+diag(omega_d.^(-1)))\(X'*Y + diag(omega_d.^(-1))*b_hat);

%% Optimize posterior for A_0
func = @(a0)-logPostA0(a0,X,Y,omega_d/lam^2,b_hat,B_hat,M,T,lam);
solver = 'fminunc';
options = optimoptions(@fminunc,'Algorithm','quasi-newton','Display', 'off');
%Solve once and obtain hessian matrix for the proposal distribution in
%Metropolis alg
a0_0 = initPoint(X,Y,B_hat,b_hat,omega_d);%a0_0 is a column vector
[a0_init,post_val_init,exitflag,output_init,grad_init,Hessian]=...
    fminunc(func,a0_0,options);

%Use global optimization toolbox to find "global" optimial for the
%posterior of A_0
tic;
filename = 'p1-1.csv';
[a0_opt,manymin] = findGlobalOpt(func,solver,options,multistarts,a0_0,filename);
A0_opt = vec2mat(a0_opt,M);
toc;

%% Construct Impulse Response Function
% Monetary policy shock
time_range = [1,ir_t];
tic;
%Sign constraint, so that the one positive monetary shock will result in
%immediate decrease of the federal funds rate (R)
imd_resp = A0_opt\mp_shock;
sign_con = sign((imd_resp(M)));
imf = sign_con*IMF(ir_t,M,p,A0_opt,B_hat,mp_shock,time_range);
toc;
h_errbd = figure(1);
hold on;
box on;
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Error bands of impulse response functions', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'FontSize',20);
legendinfo{1}=['Impulse response to Monetary Policy Shock'];
handle(1) = plot_all_TS(imf,'-k',[time_range,-inf,inf]);
%Check the results and they are close

%% Metropolis algorithm and error bands
%acceptance rate
eps_hat = Y - X*B_hat;
S_dhat = eps_hat'*eps_hat+(B_hat-b_hat)'/diag(omega_d)*(B_hat-b_hat);
tic;
[ar,trace_dat] = Metropolis_alg(Num,a0_0,C,Hessian,S_dhat,T,M);%Draw A0 from its marginal posterior distribution
toc;
imf_MCMC = zeros(ir_t,M,Num);%Store the IMF for each MCMC draw 
tic;
for i = 1:Num
    A0_temp = vec2mat(trace_dat(i,:)',M);
    inv_A0_temp = inv(A0_temp);
    imd_resp = inv_A0_temp*mp_shock;%Sign constraint as before
    sign_con = sign(imd_resp(M));
    V_temp = kron(inv_A0_temp*inv_A0_temp',inv(X'*X+diag(omega_d.^(-1))));
    B_new = reshape(mvnrnd(reshape(B_hat,K*M,1),V_temp),K,M);%Draw new B based on each draw of A0
    imf_temp = sign_con*IMF(ir_t,M,p,A0_temp,B_new,mp_shock,time_range);%Calculate IMF for this sample of A0 and B
    imf_MCMC(:,:,i) = imf_temp(:,2:end);
end
err_bd = quantile(imf_MCMC,[0.05,0.95,0.16,0.84],3);%Empirical quantile for each IMF
legendinfo{2} = ['%5 error bound'];
handle(2) = plot_all_TS([(1:ir_t)',err_bd(:,:,1)],'b-.',[time_range,-inf,inf]);
legendinfo{3} = ['%95 error bound'];
handle(3) = plot_all_TS([(1:ir_t)',err_bd(:,:,2)],'r-.',[time_range,-inf,inf]);
legendinfo{4} = ['%16 error bound'];
handle(2) = plot_all_TS([(1:ir_t)',err_bd(:,:,3)],'b.',[time_range,-inf,inf]);
legendinfo{5} = ['%84 error bound'];
handle(3) = plot_all_TS([(1:ir_t)',err_bd(:,:,4)],'r.',[time_range,-inf,inf]);
toc;
saveas(h_errbd,'Error_Bands','jpg');
saveas(h_errbd,'Error_Bands','fig');

%Trace plot
tic;
h_trace_plot = figure(2);
hold on;
box on;
annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Trace plots for free variables', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center',...
    'FontSize',20);
handle1 = plot_trace(trace_dat,6,3,'b-');
toc;
saveas(h_trace_plot,'Trace_plots','jpg');
saveas(h_trace_plot,'Trace_plots','fig');

%% Variance Decomposition
var_imf_mp = imf_var(M,p,A0_opt,B_hat,mp_shock,t_proj);
var_imf_all = imf_var(M,p,A0_opt,B_hat,ones(M,1),t_proj);
display(['The portion of variance of GDP due to monetary policy shock is:',10]);
display(var_imf_mp(1,1)/var_imf_all(1,1));