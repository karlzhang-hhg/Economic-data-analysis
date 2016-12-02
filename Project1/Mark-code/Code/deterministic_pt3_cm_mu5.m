clear
tic
%Read the Macro Data
%(1) GDP (2) GDP Deflator (3) Fed Funds (4) Consumption 
%(5) Investment (6) Hours (7) Wages
format long g
load('dataVARmedium');
M=y;
%Set mu
mu=5;
%Set number of lags and number of variables, calcluate x vector size
p=5;[~,n]=size(M);xsize=n*p+1;
%Set degrees of freedom, and tuning parameter, alpha=1 so ignore
d=n+2;lambda=0.2;
%Set number of simulations, where to start the windows
%START AT 196 SO WE ONLY DO ONCE
simnum=1000;winstart=200;
%Expanding window, starting at 1974Q4, going to 2007 Q4
for win=winstart:200
    
    %Extract the window of interest
    Msmall=M(1:win,:);T=win-p;
    %Run the AR(1) regression using the data up to this point
    varvec=zeros(n,1);
    for i=1:n
        y=M(p+2:win,i); x=M(p+1:win-1,i); v1=ones(win-p-1,1);x=[v1 x];
        beta=(x'*x)\(x'*y); %All betas correct
        resid=y-x*beta;
        varvec(i)=var(resid); %Checked this number it works
    end
    varmat=diag(varvec);
    varvecinv=1./varvec;
    bigvarvec=(10^6)/(lambda^2);
    %Pick some huge number to be flat on the constants
    for i =1:p
       bigvarvec=[bigvarvec; i^(-2).*varvecinv];
    end
    Omega=lambda^2*diag(bigvarvec);
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
    ytemp=Msmall(p+1:win,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %ADD OUR SYNTHETIC OBSERVATIONS TO THE BOTTOM%
    y0bar=zeros(n,1);
    for i=1:n
        runsum=0;
        for j=1:p
            runsum=runsum+Msmall(1,i);
        end
        y0bar(i)=runsum/p;
    end
    y0bar=y0bar*mu;
    ystar=diag(y0bar);ytemp=[ytemp;ystar];v1zero=zeros(n,1);y=vec(ytemp);
     xstar=[v1zero ystar ystar ystar ystar ystar];x=[x;xstar];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    bhat=(x'*x+inv(Omega))\(x'*ytemp+Omega\btilde);
    pi=bhat;
    %Need to create a new vector for forecasting
    yfore=zeros(xsize,1);yfore(1)=1;k=2;
    yfore1=Msmall(5,:);yfore2=vec(yfore1);yfore=[1;yfore2];
    for i=1:p
        for j=1:n
            yfore(k)=Msmall(p-i+1,j);k=k+1;
        end
    end
    %Forecast the next set of y's
    fore=yfore'*pi;  
    ynew_temp=yfore(2:xsize-7);
    ynew=[fore';ynew_temp];
    fmat=fore;
    %Want to modify yfore, remove the last rows, and add forecast to top.
    %We've done one step ahead, want to do T-p(I call it T) step ahead
    for j=2:T
        %Step 1, add constant
        ypred=[1;ynew]; 
        %Step 2, forecast
        fore=ypred'*pi;
        %Step 3, update y
        ynew_temp=ypred(2:xsize-7);
        ynew=[fore'; ynew_temp];
        %Save forecast path, add to bottom so lines up with Msmall
        fmat=[fmat;fore];  
    end
    
    %Look at test paths
    plot1=Msmall(p+1:T+p,:);plot2=fmat(:,:);
    %Export this and plot in stata
    exportmat=[plot1 plot2];
    add_actual=M(1:5,:);
    add_actual2=[add_actual add_actual];
    exportmat=[add_actual2;exportmat];
    %Put the data together in the "x" matrix
    %y=exportmat(:,1);x=exportmat(:,8); one_vec=ones(200,1); x=[one_vec x];
    %r=ols1(y,x)
    csvwrite('q1pt3_mu5.csv',exportmat);
end

toc