function [ar,trace_dat] = Metropolis_alg(Num,a0_0,C,Hessian,S_dhat,T,M,text)

if (issymmetric(Hessian)~=1)
    Hessian = (Hessian+Hessian')/2;
end
delta = min(eig(Hessian));
if (delta < 0)
    Hessian = Hessian -2*delta*eye(size(Hessian,1));
end

V = C^2*inv(Hessian);
trace_dat = zeros(Num,length(a0_0));
accept = 0;
for i = 1:Num
    A0_0 = vec2mat(a0_0,M);
    a0_new = mvnrnd(a0_0,V);
    A0_new = vec2mat(a0_new,M);
    lgp_0 = (T+M)*log(det(A0_0))-0.5*trace(S_dhat*(A0_0'*A0_0));
    lgp_new = (T+M)*log(det(A0_new))-0.5*trace(S_dhat*(A0_new'*A0_new));
    rate = exp(lgp_new-lgp_0);
    if rate >= 1
        a0_0 = a0_new;
        accept = accept + 1;
    else
        u = rand;
        if u <= rate
            a0_0 = a0_new;
            accept = accept + 1;
        end
    end
    trace_dat(i,:) = a0_0';
end
ar = accept / Num;
csvwrite([text,'trace_plot.csv'],trace_dat);
end