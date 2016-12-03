%Calculates impluse response
%Inputs: beta, shock, size of little x, number of variables and ir time
%Outputs: Matrix with impluse response path for all n variables
function [ fmat ] = Impulse_Response_Macro(bhat,shock,xsize,n,ir_time)
    %Forecast the next set of y's
    yfore=zeros(xsize,1);
    fore=yfore'*bhat+shock;  
    %Add new observations
    ynew_temp=yfore(2:xsize-n);
    ynew=[fore';ynew_temp];
    fmat=fore;
    %Want to modify yfore, remove the last rows, and add forecast to top.
    %We've done one step ahead, want to do T-p(I call it T) step ahead
    for j=2:ir_time
        %Step 1, add constant->Changed to 0 to get IR rather than determ.
        ypred=[0;ynew]; 
        %Step 2, forecast
        fore=ypred'*bhat;
        %Step 3, update y
        ynew_temp=ypred(2:xsize-n);
        ynew=[fore'; ynew_temp];
        %Save forecast path, add to bottom so lines up with Msmall
        fmat=[fmat;fore];  
    end
end

