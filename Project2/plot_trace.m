function h=plot_trace(trace_dat,n,m,patten)

[Num,~] = size(trace_dat);

for i = 1:(n*m)

subplot(n,m,i);
hold on;
box on;
h=plot(1:Num,trace_dat(:,i),patten);
% legend('Observed Data','DC Prediction')
%legend(h,text)
title(['Free variable ',num2str(i)]);

end

end