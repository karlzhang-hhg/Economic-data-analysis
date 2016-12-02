function rt = hw1ex2plot(data,p)

T = size(data,1);
n = size(data,2);

  rt1 = hw1ex2aforecast(data,200,p);
%   rt1 = hw1ex2cforecast(data,200,p,1);
  rt2 = hw1ex2bforecast(data,200,p);
  rt3 = hw1ex2cforecast(data,200,p,1);
  rt4 = hw1ex2cforecast(data,200,p,5);
  
  xt = ((1959+0.25): 0.25: (2008+1))';
%   rt2 = 
%     h = Tulim - p;
%     realpara = zeros(1,6);
%     realpara(1,1) = data(Tulim+h,1); % log real GDP
%     realpara(1,2) = (exp(data(Tulim+h,2))-exp(data(Tulim+h-1,2)))/...
%         exp(data(Tulim+h-1,2));   % price inflation
%     realpara(1,3) = data(Tulim+h,3); % federal funds rate
% 
%     % ngdp = exp(col2)*exp(col1)/100
%     realpara(1,4) = ((exp(data(Tulim+h,2))*exp(data(Tulim+h,1)))-...
%         (exp(data(Tulim+h-1,2))*exp(data(Tulim+h-1,1)))) /...
%         (exp(data(Tulim+h,2))*exp(data(Tulim+h,1))); % nominal wage inflation
% 
%     realpara(1,5) = log((exp(data(Tulim+h,2))*exp(data(Tulim+h,1))/100)*...
%         exp(data(Tulim+h,6))/(exp(data(Tulim+h,2))*exp(data(Tulim+h,1)))); 
%         % log labor share
% 
%     realpara(1,6) = log(exp(data(Tulim+h,4))/exp(data(Tulim+h,1)));
%     realmat(Tulim-59,:) = realpara;



    
figure 
title('Predicted vs Realized: T-p ahead')

subplot(3,2,1);
plot(xt,rt1.realgdp,'-k',xt,rt1.predgdp,'-g',xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
axis([1959,2008,-inf,inf])
% legend('Observed Data','DC Prediction')
title('Log Real GDP')

subplot(3,2,2);
plot(xt(2:end),rt1.realp,'-k',xt(2:end),rt1.predp,'-g',xt(2:end),rt2.predp,'-.m',xt(2:end),rt3.predp,'--c',xt(2:end),rt4.predp,':r')
axis([1959,2008,-inf,inf])
% legend('Observed Data','DC Prediction')
title('Price Inflation')

subplot(3,2,3);
plot(xt,rt1.realff,'-k',xt,rt1.predff,'-g',xt,rt2.predff,'-.m',xt,rt3.predff,'--c',xt,rt4.predff,':r')
axis([1959,2008,-inf,inf])
legend('Observed Data','Flat DC','Minnesota DC','SOC DC mu = 1','SOC DC mu = 5')
title('Federal Funds Rate')

subplot(3,2,4);
plot(xt(2:end),rt1.realnwi,'-k',xt(2:end),rt1.prednwi,'-g',xt(2:end),rt2.prednwi,'-.m',xt(2:end),rt3.prednwi,'--c',xt(2:end),rt4.prednwi,':r')
axis([1959,2008,-inf,inf])
% legend('Observed Data','DC Prediction')
title('Nominal Wage Inflation')

subplot(3,2,5);
plot(xt,rt1.realls,'-k',xt,rt1.predls,'-g',xt,rt2.predls,'-.m',xt,rt3.predls,'--c',xt,rt4.predls,':r')
axis([1959,2008,-inf,inf])
% legend('Observed Data','DC Prediction')
title('Log Labor Share')

subplot(3,2,6);
plot(xt,rt1.reallcr,'-k',xt,rt1.predlcr,'-g',xt,rt2.predlcr,'-.m',xt,rt3.predlcr,'--c',xt,rt4.predlcr,':r')
axis([1959,2008,-inf,inf])
% legend('Observed Data','DC Prediction')
title('Log Consumption Ratio')

annotation('textbox', [0 0.9 1 0.1], ...
    'String', 'Observed vs DC with Different Priors', ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center')

display(rt1.rsq);
display(rt2.rsq);
display(rt3.rsq);
display(rt4.rsq);

rt = 0;




