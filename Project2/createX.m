function X = createX(Dat,p,T,M,K)
%Assemble the ordinary X
X = zeros(T,K);% a T*K matrix
for j = 1:T
    X(j,:) = [1,reshape(Dat((j+p-1:-1:j),:)',1,K-1)];% have to transpose because the reshape function operate in column
end
%Add sum-of-coefficient

end