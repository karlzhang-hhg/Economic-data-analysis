function h=plot_all_TS(TSdat,patten,axlmt)
    subplot(3,2,4);
    hold on;
    box on;
    grid on;
    plot(abs(TSdat(:,1)),TSdat(:,2),patten);
    axis(axlmt)
    % legend('Observed Data','DC Prediction')
    %legend(h,text)
    title('Y')

    subplot(3,2,5);
    hold on;
    box on;
    grid on;
    plot(abs(TSdat(:,1)),TSdat(:,3),patten);
    axis(axlmt)
    % legend('Observed Data','DC Prediction')
    title('P')
    
    subplot(3,2,6);
    hold on;
    box on;
    grid on;
    plot(abs(TSdat(:,1)),TSdat(:,4),patten);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis(axlmt)
    % legend('Observed Data','DC Prediction')
    title('U')
    
    subplot(3,2,1);
    hold on;
    box on;
    grid on;
    plot(abs(TSdat(:,1)),TSdat(:,5),patten);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis(axlmt)
    % legend('Observed Data','DC Prediction')
    title('Pcom')
    
    subplot(3,2,2);
    hold on;
    box on;
    grid on;
    plot(abs(TSdat(:,1)),TSdat(:,6),patten);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis(axlmt)
    % legend('Observed Data','DC Prediction')
    title('M2')
    
    subplot(3,2,3);
    hold on;
    box on;
    grid on;
    h=plot(abs(TSdat(:,1)),TSdat(:,7),patten);%,xt,rt2.predgdp,'-.m',xt,rt3.predgdp,'--c',xt,rt4.predgdp,':r')
    axis(axlmt)
    % legend('Observed Data','DC Prediction')
    title('R')
end