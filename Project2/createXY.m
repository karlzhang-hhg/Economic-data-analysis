function [X,Y] = createXY(Dat,p,T,M,K,mu)
%Assemble the ordinary Y
Y = Dat((p+1:p+T),:);
%Assemble the ordinary X
X = zeros(T,K);% a T*K matrix
for j = 1:T
    X(j,:) = [1,reshape(Dat((j+p-1:-1:j),:)',1,K-1)];% have to transpose because the reshape function operate in column
end

%Add sum-of-coefficient
Y_bar = diag(mu*mean(Dat(1:p,:),1));
Y = [Y;Y_bar];
Xplus = [zeros(M,1),repmat(diag(mu*mean(Dat(1:p,:),1)),[1 p])];    
X= [X;Xplus]; % stack x at the bottom
end