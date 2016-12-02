%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This program estimates the Flat Prior VAR for Part 1 of Excersise 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This uses the technique developed in class, as opposed to the technique
%outlined in Hamilton.  All numbers checked manually
clear
tic
%Read the Macro Data
%(1) GDP (2) GDP Deflator (3) Fed Funds (4) Consumption 
%(5) Investment (6) Hours (7) Wages
M=csvread('macrodata.csv');
%Set number of lags and number of variables, calcluate x vector size
p=5;[~,n]=size(M);xsize=n*p+1;
%Set degrees of freedom, and tuning parameter, alpha=1 so ignore
d=n+2;
%Set number of simulations, where to start the windows
simnum=1000;winstart=64;
optlams=zeros(196-winstart+1,1);
lamstart=0.05;lamend=0.25;
lamstep=0.0001;steps=1+(lamend-lamstart)/lamstep;
%Expanding window, starting at 1974Q4, going to 2007 Q4
for win=winstart:196
    %Extract the window of interest
    Msmall=M(1:win,:);T=win-p;
    %Run the AR(1) regression using the data up to this point
    varvec=zeros(n,1);
    for i=1:n
        y=M(2:win,i); x=M(1:win-1,i); v1=ones(win-1,1);x=[v1 x];
        beta=(x'*x)\x'*y; %All betas correct
        resid=y-x*beta;
        varvec(i)=var(resid); %Checked this number it works
    end
    varmat=diag(varvec);
    varvecinv=1./varvec;
    bigvarvec=(10^6);%Pick some huge number to be flat on the constants
    for i =1:p
       bigvarvec=[bigvarvec; i^(-2).*varvecinv];
    end
    
    %Omega(1)=0; %Thought this was messing up the result but it's not
    %Create the unit root prior vector b
    b=zeros(n*(n*p+1),1);k=1;index=1;
    %Loop over equations
    for j=1:n
       b(k)=0;k=k+1; %Constant
       %Loop over terms in equation
       for i=1:(n*p)
           if i==index
               b(k)=1;k=k+1;
           else
               b(k)=0;k=k+1;
           end
       end
       index=index+1;
    end    
    %Turn b into a matrix
    btilde=reshape(b,[n*p+1,n]);
    %Calculate Bhat, note T=59 for first set, so x is TxK, y is Txn
    %Break up the data, vector of ones, and the different parts of x
    v1=ones(T,1);x1=Msmall(p:win-1,:);x2=Msmall(p-1:win-2,:);
    x3=Msmall(p-2:win-3,:);x4=Msmall(p-3:win-4,:);x5=Msmall(p-4:win-5,:);
    %Put the data together in the "x" matrix
    x=[v1 x1 x2 x3 x4 x5];
    %Create the y matrix
    y=Msmall(p+1:win,:);
    
    ll=zeros(steps,2);
    ll_index=1;
    for lambda=lamstart:lamstep:lamend
        %Calculate omega with new lambda
        Omega=lambda^2*diag(bigvarvec);
        %Calculate the bhat
        bhat=(x'*x+inv(Omega))\(x'*y+Omega\btilde);
        %Calculate ehat
        ehat=y-x*bhat;
        %Calculate Shat and Shathat. Y/X means multiply by the inverse of X as
        %opposed to Y\X which means multiply by the inverse of Y
        shat=ehat'*ehat+(bhat-btilde)'/Omega*(bhat-btilde);
        shathat=shat+varmat; %varmat=psi
        %Get T
        [T,~]=size(y);
        %Get variance
        varhat=inv(x'*x+inv(Omega));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Want to calculate P(Y|hyperparemeters), to optimally select lambda
        %Cholesky factor Omega and Psi
        %Romega=chol(Omega); Rpsi=chol(inv(varmat));
        %NOTE: Have the multivariate gamma Y = gampdf(X,A,B), this is just a 
        %constant so it won't change the mode, leave out for now
        %((1/pi)^(n*T/2))*(gamma((T+d)/2)/gamma(d/2))*det(varmat)^(-T/2)
        %Numerical issues, check to see if INF or 0
        %e1=eig(Romega'*(x'*x)*Romega); e1=1+e1; e1p=prod(e1);
        %e2=eig(Rpsi'*(ehat'*ehat+(bhat-btilde)'/Omega*(bhat-btilde))*Rpsi); 
        %e2=1+e2; e2p=prod(e2);
        ll(ll_index,1)=lambda;
        %e1p^(-n/2)si
        %e2p^(-(T+d)/2)
        
        %Numerical issues with computing the likelihood, so just took
        %logML from the code provided online
        omega=diag(Omega);psi=varvec;epshat=ehat;betahat=bhat;b=btilde;
        
        aaa=diag(sqrt(omega))*(x'*x)*diag(sqrt(omega));
        bbb=diag(1./sqrt(psi))*(epshat'*epshat + (betahat-b)'*diag(1./omega)*(betahat-b))*diag(1./sqrt(psi));

        eigaaa=real(eig(aaa)); 
        eigaaa(eigaaa<1e-12)=0; 
        eigaaa=eigaaa+1;
        eigbbb=real(eig(bbb)); 
        eigbbb(eigbbb<1e-12)=0; 
        eigbbb=eigbbb+1;

        logML = - n*T*log(pi)/2 + sum(gammaln((T+d)/2)-gammaln(d/2)) +...
            - T*sum(log(psi))/2 - n*sum(log(eigaaa))/2 - (T+d)*sum(log(eigbbb))/2;
        ll(ll_index,2)=logML;
        %ll(ll_index,2)=(-n/2)*log(e1p)-((T+d)/2)*log(e2p);
        %Direct computation not stable, get all infinity
        %ll(ll_index,2)=log(det(Omega)^(-n/2))+log((det((x'*x)+inv(Omega))^(-n/2)))+log(det(shathat)^(-(T+d)/2));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        ll_index=ll_index+1;
    end
    %Just getting a bunch of zeros
    ll;
    [~,I]=max(ll);lambda=ll(I(2),1);optlams(win-winstart+1)=lambda;
    %Calculate omega with new lambda
    Omega=lambda^2*diag(bigvarvec);
    %Calculate the bhat
    bhat=(x'*x+inv(Omega))\(x'*y+Omega\btilde);
    %Calculate ehat
    ehat=y-x*bhat;
    %Calculate Shat and Shathat. Y/X means multiply by the inverse of X as
    %opposed to Y\X which means multiply by the inverse of Y
    shat=ehat'*ehat+(bhat-btilde)'/Omega*(bhat-btilde);
    shathat=shat+varmat; %varmat=psi
    %Get T
    [T,~]=size(y);
    %Get variance
    varhat=inv(x'*x+inv(Omega));    
    
    %Need to create a new vector for forecasting
    yfore=zeros(xsize,1);yfore(1)=1;k=2;
    for i=1:p
        for j=1:n
            yfore(k)=Msmall(win-i+1,j);k=k+1;
        end
    end
    %Forecast the next set of y's
    fore=yfore'*bhat;
    %Want to modify yfore, remove the last rows, and add forecast to top
    for j=2:4
        yfore_temp=yfore(2:end-(j-1)*n);
        yfore_temp2=[1;fore';yfore_temp];
        forecast=yfore_temp2'*bhat;
        %Add new forecast to old
        fore=[forecast fore];        
    end
    %Make each row a forecast for each time period
    f=reshape(fore,[n,4]);
    %Slice data from master, first obs is T, 5th is T+4
    y_sample=M(win:win+4,1:2);
    %Calculate implied growth rate, only do this for first 2 variables 
    %The first column of y_sample has most recent dat
    g=zeros(1,4);
    %GDP Growth  - Checked against data
    g(1)=f(1,4)-y_sample(1,1);g(2)=0.25*(f(1,1)-y_sample(1,1));
    %GDP Deflator growth  - Checked against data
    g(3)=f(2,4)-y_sample(1,2);g(4)=0.25*(f(2,1)-y_sample(1,2));
    %Get true growth rates
    tg=zeros(1,4);
    %GDP growth - Checked against data
    tg(1)=y_sample(2,1)-y_sample(1,1);tg(2)=0.25*(y_sample(5,1)-y_sample(1,1));
    %GDP Deflator Growth - Checked against data
    tg(3)=y_sample(2,2)-y_sample(1,2);tg(4)=0.25*(y_sample(5,2)-y_sample(1,2));
    %Calculate error
    e=zeros(1,4);
    for i=1:4
        e(i)=g(i)-tg(i);
    end
    if win==winstart
       emat=e; 
    else
        emat=[emat; e];
    end
end
emat2=emat.^2;
[s1,s2]=size(emat2); mse=zeros(1,4);
%Compute mean squared error
mse=(1/s1)*sum(emat2)
csvwrite('mse_opt.csv',mse)
csvwrite('lambdas.csv',optlams)
toc
plot(optlams)