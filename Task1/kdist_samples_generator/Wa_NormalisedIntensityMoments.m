function w_outfNormalisedIntensityMoments = Wa_NormalisedIntensityMoments(varargin )
%returns the intensity moments up to order w_pnOrder, 
%claculate intensity moments. 
%2nd arg is order
% if the last arugmement is false the data is not plotted

%% calculating arguments w_pfDataIn,  w_pnOrder, plot or not

w_pbPlot= true;

w_pnOrder = 4;
switch nargin
    case {1} 
        
        w_pfDataIn = varargin{1};
    case {2}
        w_pfDataIn = varargin{1};  
        w_pnOrder = varargin{2};
   case{3}
        w_pfDataIn = varargin{1};    
        w_pnOrder = varargin{2};
        w_pbPlot = varargin{3};
        
    otherwise 
        error(' Number of args between 1 and 3');
end 



%% calculate output
w_pfNumberSamples=length(w_pfDataIn);
w_outfNormalisedIntensityMoments = zeros(1, w_pnOrder);
w_outfNormalisedIntensityMoments(1) = sum(w_pfDataIn)/w_pfNumberSamples;

for w_nidx=2:w_pnOrder    
    w_outfNormalisedIntensityMoments(w_nidx) = sum(w_pfDataIn.^w_nidx) / w_pfNumberSamples/ w_outfNormalisedIntensityMoments(1)^w_nidx;  
end

% w_outfNormalisedIntensityMoments(1)=1;

% %% plot the results
% if (w_pbPlot~=false)
%     
% w_nrange=1:w_pnOrder;
% 
% plot(factorial(w_nrange),'r')
% hold on
% plot(w_outfNormalisedIntensityMoments);
% set(gca,'XTick', w_nrange)
% hold off
% grid on
% 
% legend('Normaliosed intensity moments of Gaussian noise', 'Normalised moments of the data')
% 
% end

    



