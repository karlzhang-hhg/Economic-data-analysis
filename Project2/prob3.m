clear();
%Load the data set
Dat_all = xlsread('SZdataExtended.xlsx');
Dat = Dat_all(:,2:7);
%Set constants
p = 13;%Maximum lag
lam = 0.2;%Hyperparameter
M = 6;%Dimension of vector
T = size(Dat,1)-p;%Time for forcasting
K = M*p+1;
vc = 10^6;%First element in Minnesota prior
mu = 1;
Num = 500000;%Number of MCMC simulation
C = 0.5;%Scaling parameter for the proposal distribution variation
ir_t = 48;%Impulse response function horizon
mp_shock = zeros(M,1);%Monetary policy shock
mp_shock(M) = 1;
multistarts = 5000;%Number of starts
t_proj = 36;%Variance decomposition
text = 'Prob 3-2-';
% Use the entire data set
IMF_Errbd(Dat,p,lam,M,T,K,vc,mu,Num,C,ir_t,mp_shock,multistarts,t_proj,text);
% Excluding ZLB period
Dat = Dat_all(1:600,2:7);
T = size(Dat,1)-p;%Time for forcasting
text = 'Prob 3-1-';
IMF_Errbd(Dat,p,lam,M,T,K,vc,mu,Num,C,ir_t,mp_shock,multistarts,t_proj,text);
