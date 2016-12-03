function h=plot_all_TS(TSdat,pattern,text)
    subplot(3,2,1);
    hold on;
    box on;
    h=plot(dat1(:,1),dat1(:,2),pattern);
    axis([1959,2008,-inf,inf])
    % legend('Observed Data','DC Prediction')
    %legend(h,text)
    title('log-real GDP')

    subplot(3,2,2);
    hold on;
    box on;
    plot(dat2(:,1),dat2(:,2),pattern);
    axis([1959,2008,-inf,inf])
    % legend('Observed Data','DC Prediction')
    title('quarterly price inflation')
    
    subplot(3,2,3);
    hold on;
    box on;
    plot(dat3(:,1),dat3(:,2),pattern);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis([1959,2008,-inf,inf])
    % legend('Observed Data','DC Prediction')
    title('the federal funds rate')
    
    subplot(3,2,4);
    hold on;
    box on;
    plot(dat4(:,1),dat4(:,2),pattern);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis([1959,2008,-inf,inf])
    % legend('Observed Data','DC Prediction')
    title('quarterly nominal wage inflation')
    
    subplot(3,2,5);
    hold on;
    box on;
    plot(dat5(:,1),dat5(:,2),pattern);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis([1959,2008,-inf,inf])
    % legend('Observed Data','DC Prediction')
    title('log-labor share')
    
    subplot(3,2,6);
    hold on;
    box on;
    plot(dat6(:,1),dat6(:,2),pattern);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis([1959,2008,-inf,inf])
    % legend('Observed Data','DC Prediction')
    title('log-consumption ratio')
end