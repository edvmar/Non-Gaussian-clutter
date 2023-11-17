%Matthew Ritchie
%08/04/09

%Fitting function

function [KdistShape,total_sum_sqr_diff,Kdist_PFAdB_fit] = RDV_fit_05(X,pfadB,Y,nu_array,RDV_overlay_plots,RDV_pfa_seperate_Cmp)

size_pfa=size(pfadB);
pfadB=pfadB(1:size_pfa(1)-1,:);
size_pfa = size(pfadB);
X = X((1:length(X)-1));
Y = Y((1:length(Y)-1));

% figure
% plot(X ,Y);

sum_sqr_diff = zeros(1,size_pfa(2));
for count = 1:size_pfa(2)
    sqr_diff = (Y' - pfadB(:,count)).^2;
    sum_sqr_diff(count) =sum(sqr_diff);
end
[fitdif fitnu] = min(sum_sqr_diff);
total_sum_sqr_diff = sum(sum_sqr_diff);
KdistShape = nu_array(fitnu);
Kdist_PFAdB_fit = pfadB(:,fitnu)';


if RDV_overlay_plots ==1;
    figure
    %Plot all the pfa curves in green
    %     plot(X,pfadB(:,1:size_pfa(2))','g-');
    %Plot the raw data pfa curve
    plot(X,Y,'b-');
    hold on
    %Plot the fitted pfa curve in red.
    
    plot(X,pfadB(:,fitnu)','r-');
    
    axis([0 35 -5 0])
    if RDV_pfa_seperate_Cmp==0
        legend('Raw Data',sprintf('Fitted Shape Parameter %0.5g',nu_array(fitnu)))
    end        
    xlabel('Threshold (dB)')
    ylabel('Log10(Pfa))')    
end

end



