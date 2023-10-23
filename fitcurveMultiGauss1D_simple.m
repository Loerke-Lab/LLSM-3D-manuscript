function estimates = fitcurveMultiGauss1D_simple(x, data, guess, fix)
% INPUT:    x       = x data
%
%           data    = y data (should have same length as x)
%
%           guess   = guess vector for parameter values; this input is
%           required, as it defines the number of Gaussian populations
%           The vector should have the length 3*n where n is the number of
%           populations (peaks); for each individual population, the three
%           parameters are (offset is assumed as zero for all):
%           [An xn sn] where An = amplitude, xn = max position, sn = sigma
%           Thus, for two Gaussians, the guess vector should be
%           [A1 x1 s1 A2 x2 s2]
%
%           fixvector (optional) = fixing vector for fitting; this vector
%           needs to have the same length as guess, and its values should
%           be either 0/1 to specify whether the parameter at this position
%           is fixed (in which case it is reset to the guess value after
%           each fitting iteration) or free (in which case it's fit freely)  
%
% OUTPUT:   estimates = fit results, following the same format as guess
%
% Call fminsearch with a specified starting point.
%
% NOTE:     In its current implementation, this function uses the amplitude
%           A to scale each Gaussian; the individual Gaussian peaks are not
%           normalized (by diving by sqrt(2*pi*sigma^2) first. As a result,
%           the amplitude values represent the absolute *height* of the
%           Gaussian peaks, not the *area* under the peaks. 

start_point = guess;

global fixvec
fixvec = 0*guess;

if nargin>3
    fixvec = fix;
end


% order is: offset amplitude centerpoint width
estimates = fminsearch(@Gauss1Dfun, start_point);


% expfun accepts curve parameters as inputs and outputs sse,
% the sum of squares error for A * exp(-lambda * t) - Data.
    function sse = Gauss1Dfun(params)
        
        %global fixvec
        
        if max(fixvec)>0
            fixpos = find(fixvec==1);
            params(fixpos) = guess(fixpos);
        end
        
        np = length(params);
        nump = (np)/3;
        amp    = params(1:3:np);
        x0     = params(2:3:np);
        sig    = abs(params(3:3:np));
        
        for k=1:nump
            IndGauss(k,:) = amp(k) .* exp( - ((x-x0(k)).^2) / (2*sig(k)^2) );
        end
        
        FittedCurve = sum(IndGauss,1);
        
        if max(fixvec)>0
            fixpos = find(fixvec==1);
            params(fixpos) = guess(fixpos);
        end
        
        %hold off
        %plot(x,FittedCurve)
        %pause(0.01)
        
        ErrorVector = FittedCurve - data;
        sse = sum(ErrorVector .^ 2);      
        
    end

% return estimates to fixed guess positions where necessary
if max(fixvec)>0
    fixpos = find(fixvec==1);
    estimates(fixpos) = guess(fixpos);
end

% plot fitting results: mixed population in red
plot(x,data,'bo','MarkerSize',4); hold on;
%plot(x,FittedCurve,'r-');

% plot fitting results: individual subpopulations in green
np = length(estimates);
nump = (np)/3;
amp    = estimates(1:3:np);
x0     = estimates(2:3:np);
sig    = abs(estimates(3:3:np));
%sigAvg = mean(sig); % added 5/3/23 KB
    
SumPopulations = zeros(size(x));
for k=1:nump
    IndivPopulation = amp(k) .* exp( - ((x-x0(k)).^2) / (2*sig(k)^2) );
    %IndivPopulation = amp(k) .* exp( - ((x-x0(k)).^2) / (2*sigAvg^2) ); % added 5/3/23 KB for equal sigmas %
    plot(x,IndivPopulation,'LineWidth',2);
    
    SumPopulations = SumPopulations + IndivPopulation;
end
plot(x,SumPopulations,'LineWidth',2,'Color',[0.4940 0.1840 0.5560]);
%legend('Data','Fit 1','Fit 2','Sum Fits')
yticks([0 0.02 0.04 0.06 0.08 0.1 0.12 0.14])
yticks('manual')
grid on
legend('Data','Fit 1','Fit 2','Sum Fits','Location','northoutside')
ylim([0 0.15])

end % of function