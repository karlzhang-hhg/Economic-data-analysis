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
p=5;[~,n]=size(M);xsize=n*p+1;winstart=64;
%Expanding window, starting at 1974Q4, going to 2007 Q4
for win=winstart:196
    %Extract the window of interest
    Msmall=M(1:win,:);T=win-p;
    %Break up the data, vector of ones, and the different parts of x
    v1=ones(T,1);x1=Msmall(p:win-1,:);x2=Msmall(p-1:win-2,:);
    x3=Msmall(p-2:win-3,:);x4=Msmall(p-3:win-4,:);x5=Msmall(p-4:win-5,:);
    %Put the data together in the "x" matrix
    x=[v1 x1 x2 x3 x4 x5];
    %det(x'*x)
    %Create the y matrix and vectorize
    ytemp=Msmall(p+1:win,:);y=vec(ytemp);
    %Compute the matrix of coefficents, and vectorize
    pi=x'*x\x'*ytemp;Pi=vec(pi);
    %Calculate predicted values
    I=eye(n);X = kron(I,x);yhat=X*Pi;
    %Calculate errors
    epsilon=y-yhat;e = reshape(epsilon,[T,n]);SSR=e'*e;
    %Need to create a new vector for forecasting
    yfore=zeros(xsize,1);yfore(1)=1;k=2;
    for i=1:p
        for j=1:n
            yfore(k)=Msmall(win-i+1,j);k=k+1;
        end
    end
    %Forecast the next set of y's
    fore=yfore'*pi;
    %Want to modify yfore, remove the last rows, and add forecast to top
    for j=2:4
        yfore_temp=yfore(2:end-(j-1)*n);
        yfore_temp2=[1;fore';yfore_temp];
        forecast=yfore_temp2'*pi;
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
csvwrite('mse_flat.csv',mse)
toc