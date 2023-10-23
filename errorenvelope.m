function errorenvelope(x,y,err,faceColor,alphaVal)
% errorenvelope generates a symmetric error envelope, y-err to y+err, over
% the user-input x-vector.
%
% INPUTS:               x: vector of x-values
%                       y: vector of y-values which is the mean of the
%                          envelope
%                     err: vector of error values which will be half the
%                          width of the envelope ( y +/- err).
%               faceColor: optional input of the envelope color, dfault is
%                          blue
%                alphaVal: transparency value where 1 is opaque and zero is
%                          transparent, the default is 0.2;
%
%


% if no facecolor is given as input, set to blue
if nargin < 5
    faceColor = [0 0 1];
end
if nargin < 6
    alphaVal = 0.2;
end

% we use an error envelope that is symmetric about the mean
ymin = y - err;
ymax = y + err;

% if NaN's exist in the y-values, take just the finite values
finitePts = isfinite(ymin);
xfinite = x(finitePts);
yminfinite = ymin(finitePts);
ymaxfinite = ymax(finitePts);

x = xfinite;
ymin = yminfinite;
ymax = ymaxfinite;

X=[x,fliplr(x)];                %#create continuous x value array for plotting
Y=[ymin,fliplr(ymax)];              %#create y values for out and then back
h = fill(X,Y,'b');   
set(h,'FaceColor',faceColor,'EdgeColor','none')


% make the filled area transparent, this can be adjusted by changing the
% value: 1 is opaque, and 0 is completely transparent
alpha(alphaVal)

end
