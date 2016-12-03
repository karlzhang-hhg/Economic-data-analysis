clear();
%Load the data set
Dat = xlsread('SZdata.xlsx');
Dat = Dat(:,2:7);
%Set constants
p = 13;%Maximum lag
lam = 0.2;%Hyperparameter
M = 6;%Dimension of vector
T = size(Dat,1)-p;%Time for forcasting
K = M*p+1;
vc = 10^6;%First element in Minnesota prior
mu = 1;
Num = 1000000;%Number of MCMC simulation
C = 0.5;%Scaling parameter for the proposal distribution variation
ir_t = 48;%Impulse response function horizon
mp_shock = zeros(M,1);%Monetary policy shock
mp_shock(M) = 1;
multistarts = 2500;%Number of starts
t_proj = 36;%Variance decomposition
text = 'Prob 1-';

IMF_Errbd(Dat,p,lam,M,T,K,vc,mu,Num,C,ir_t,mp_shock,multistarts,t_proj,text);