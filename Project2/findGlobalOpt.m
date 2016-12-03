function [a0_opt, manymin] = findGlobalOpt(objFunc,solver,options,multistarts,a0_0,filename)
%% Use Global Optimization Toolbox in Matlab to find a "global" optimal of
%the objective function.
%Input: 
    %objFunc: The objFunc has to have defined the variables and parameters
%========================================
    %Input:
        %lam: lambda;
        %A0: The A_0 matrix (M-by-M);
        %X: X matrix ((T+M)-by-K);
        %Y: Y matrix ((T+M)-by-M);
        %d_phi: Diagonal entries for square root of Phi matrix (K entries)
        %b_hat: Matrix form for the mean of prior for beta (which is a vector form of B);
        %M: Dimension of vector variables;
        %T: Time periods for forcasting;
        %p: The maximum lag;
    %Output:
        %logfval: The log posterior function value;
%========================================
    %solver: The specified solver. Here use 'fminunc';
    %options: Options for the solver (optimoptions(@fminunc,'Algorithm','quasi-newton','Display', 'off')).
    %          Notice that here we are necessarily requiring gradient of objFunc
    %multistarts: Number of multistart;
    %a0_0: The initial guess for a_0;
    %filename: The file name for writing the resutls
%Output:
    %a0_opt: The optimal a0 with the optimal function value.
%% ==========================================================================
%Create an optimization problem
prob = createOptimProblem(solver,'objective',objFunc,'x0',a0_0,'options',options);
%Number of multistart
MS = MultiStart;
%Run multistarts times of this problem
[x,f,~,~,manymin] = run(MS,prob,multistarts);
%Write to a file and get the optimal value for free variables
csvwrite(filename,x);
a0_opt = x;
end