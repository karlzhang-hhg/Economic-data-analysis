function a = mat2vec(A)
[n,m] = size(A);
A = reshape(A,n*m,1);
a = A([1,2,3,4,5,8,9,10,11,15,16,22,28,29,30,34,35,36]);
end