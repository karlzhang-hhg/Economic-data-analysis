function lam_opt = lam_mle(objFunc,lam_bd,N_try,X,Y,d_psi,d_phi,b_hat,M,T,d,options_optim)
% Use fminunc to find the optimal solution for lam with several starting points
% Input:
% objFunc: objective function which return kernel of marginal likelihood
% (or log of it)
%===============================
    % Calculate the marginal likelihood
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
    % val: the marginal likelihood up to a constant factor
%===============================
% lam_bd: lambda value upper and lower bound for initial value exploring
%===============================
    % [lam_lower_bound,lam_upper_bound]
%===============================
% N_try: how many initial tries to do inside the lam_bd
% X: matrix X in time series (lag);
% Y: matrix Y which are to be forecasted;
% d_psi: diagnoal elements of Psi matrix (which is diagonal);
% d_ome: diagonal elements of Omega matrix (which is diagonal);
% M: the dimension of y
% T: time horizon to be forecasted
% d: some constant specified in prior
% options_optim: options when calling the fminunc function
% Output:
% opt_lam: optimal lambda that maximizing the marginal likelihood

%% code
%options_optim = optimoptions(@fminunc,'Display','iter-detailed','Algorithm','trust-region','SpecifyObjectiveGradient',true);

fun = @(lam)objFunc(lam,X,Y,d_psi,d_phi,b_hat,M,T,d);


lam_init = lam_bd(1):(lam_bd(2)-lam_bd(1))/N_try:lam_bd(2);
rec_lam = zeros(length(lam_init),3);
for i = 1:length(lam_init)
    %rec_lam(i,:)=fminunc(fun,lam_init(i),options_optim);
    fminunc(fun,lam_init(i),options_optim)
end
[~,opt_ind] = max(rec_lam(:,2));
lam_opt = rec_lam(opt_ind,1);
end