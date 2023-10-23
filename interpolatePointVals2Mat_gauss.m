function [ result ] = interpolatePointVals2Mat_gauss( valuematrix, binsX, binsY, sigmaX, sigmaY, displayVariable )
% construct a value matrix using values at unevenly spaced locations; the 
% parameter vectors is interpolated using gaussian weighting of distances; 
% in addition it is possible to weight contributions using an external
% weight vector
%
% INPUT:    valuematrix     n x 3 matrix containing:
%               first column:   features' x-coordinates (e.g. distance)
%               second column:  features' y-coordinates (e.g. angle)
%               third column:  features' value (e.g. correlation magnitude)
%               fourth column: (optional) weight
%                           
%
%           binsX           vector of the desired value bins of the
%                           x-vector to which the value will be
%                           interpolated
%           binsY           vector of the desired value bins of the
%                           y-vector to which the value will be
%                           interpolated
%
%           sigmaX          sigma of the Gaussian weighting function in x
%           sigmaY          sigma of the Gaussian weighting function in y
%                           
% OUPUT:   result          n x m matrix containing:
%           rows:       interpolated values for x bins (for given y)
%           columns:    interpolated values for y bins (for given x)
%
% last modified 10/07/2014 Dinah Loerke

display = 0;
if nargin>5
    if displayVariable==1
        display = 1;
    end
end

%   extract x,y coordinate and parameter vectors from input
xvec = valuematrix(:,1);
yvec = valuematrix(:,2);
zvec = valuematrix(:,3);

weightvec_sig = 1+0*xvec;
if size(valuematrix,2)>3
    weightvec_sig = valuematrix(:,4);
end


% loop over x-bins
for i=1:length(binsX)
    
    % display progress
    % display(['processing mask section ',num2str(i),' of ', num2str(length(stepvector))]);
    
    % distance vector 1 (from current x bin) for each data point
    distvec1 = abs(xvec-binsX(i));
    
    % weighting vector for each data point is Gaussian function with 
    % distance from bin 
    weightvec_dist1 = exp(-(distvec1.^2)/(2*sigmaX^2));
    
    % loop over y-bins
    for k=1:length(binsY)
        
        % distance vector (from current y-bin) for each data point
        distvec2 = abs(yvec-binsY(k));

        % weighting vector for each data point is Gaussian function with 
        % distance from bin (multiplying factors for both bins)
        weightvec_dist2 = exp(-(distvec2.^2)/(2*sigmaY^2));
        
        weightvec_distTotal = weightvec_dist1 .* weightvec_dist2;

        nanpos = isnan(zvec);
        weightvec_distTotal(nanpos) = nan;

        % multiply this weight vector with the 'external' weights, which may
        % represent e.g. statistical significance of the point
        weightvec = weightvec_distTotal.*weightvec_sig;

        % weighted average of the data parameter (yvec) for this bin
        valuvec = nansum((zvec.*weightvec),1)./nansum(weightvec,1);

        % write results (averaged parameter values for the current pixels) into
        % the results at the specified positions 
        result(i,k) = valuvec;
        
    end
    
end % of for i-loop

% display results
if display == 1
    figure; 
    subplot(1,2,1);
    plot(xvec,zvec,'b.');
    subplot(1,2,2);
    plot(yvec,zvec,'r.');
    
    figure;
    imshow(result,[]); colormap jet; set(gca,'Position',[0 0 1 1]);
end

end

