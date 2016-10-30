clear
tic
format long g
load('dataVARmedium');
M=y;
p=5;[~,n]=size(M);xsize=n*p+1;winstart=200;
%Only estimate using whole sample with winstart=200
for win=winstart:200
    Msmall=M(1:win,:);T=win-p;
    %Break up the data, vector of ones, and the different parts of x
    v1=ones(T,1);x1=Msmall(p:win-1,:);x2=Msmall(p-1:win-2,:);
    x3=Msmall(p-2:win-3,:);x4=Msmall(p-3:win-4,:);x5=Msmall(p-4:win-5,:);
    %Put the data together in the "x" matrix
    x=[v1 x1 x2 x3 x4 x5];
    %Create the y matrix and vectorize
    ytemp=Msmall(p+1:win,:);y=vec(ytemp);
    %Compute the matrix of coefficents, and vectorize
    pi=inv(x'*x)*(x'*ytemp);
    Pi=vec(pi);
    %Calculate predicted values
    I=eye(n);X = kron(I,x);yhat=X*Pi;
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
    %y=exportmat(:,1);x=exportmat(:,8); 
    %[time,vars]=size(x);
    %one_vec=ones(time,1); x=[one_vec x];
    %r=ols1(y,x);
    csvwrite('q1pt1.csv',exportmat);
end
toc