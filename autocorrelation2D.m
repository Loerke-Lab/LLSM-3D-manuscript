function [ acmat ] = autocorrelation2D( intmat, range_row, range_column )
%2-dimensional autocorrelation function
% INPUT:    intmat:     intensity matrix; this function assumes that this
%                       input matrix is already normalized appropriately
%                       and that inappropriate zeros are converted to nans
%           range_row:  value range for the dimension corresponding to the
%                       rows in the input matrix intmat, e.g. 2 (for -2:2)
%           range_column: value range for the dimension corresponding to
%                       the columns in intmat, use e.g. 5 (for -5:5)
%
% OUTPUT:   acmat:      2-dimensional autocorrelation matrix
%                       dimensions [2*range_row+1 2*range_column+1]
%
% last modified DL Dec 10,2103

% determine size of input matrix of image
[size_row, size_column] = size(intmat);

% initialize output matrix
acmat = zeros(2*range_row+1,2*range_column+1);

% loop over all desired row shifts 
for shift_row = -range_row:range_row
    
    % loop over all desired column shifts
    for shift_column = -range_column:range_column
        
        % for each specific shift in row and column direction, two partial
        % matrices are cropped out of the original matrix, which are of
        % equal size and where the second matrix is shifted by shift_row in
        % the row direction and by shift_column in column direction.
        % For example: for a shift of shift_row=-1 and shift_column=2, the 
        % second matrix is cropped such that the intensity at 
        % intensity(i,j) will be correlated to the intensity of 
        % intensity(i-1,j+2)
        
        % DETERMINE APPROPRIATE CROP MARGINS FOR BOTH PARTIAL MATRICES
        % in max and min expressions: first argument reflects the case if
        % shift_row or shift_column are negative; the second argument
        % reflects the case if the shifts are positive
        rowstart1       = max(1-shift_row,1);
        rowend1         = min(size_row,size_row-shift_row);
        rowstart2       = max(1,1+shift_row);
        rowend2         = min(size_row+shift_row,size_row);
        
        columnstart1   = max(1-shift_column,1);
        columnend1     = min(size_column,size_column-shift_column);
        columnstart2   = max(1,1+shift_column);
        columnend2     = min(size_column+shift_column,size_column);
        
        % cropped matrix 1  
        area1 = intmat(rowstart1:rowend1,columnstart1:columnend1);
             
        % cropped matrix 2 
        area2 = intmat(rowstart2:rowend2,columnstart2:columnend2);
        
        % product of shifted*unshifted
        intProduct = area1.*area2;
        % enter result into autocorrelation matrix
        acmat(range_row+shift_row+1,range_column+shift_column+1) = nanmean(intProduct(:));
        % acmat(range_row-shift_row+1,range_column-shift_column+1) = nanmean(intProduct(:));
    end % of for
end % of for

% normalize matrix by 0-0 value
acmat = acmat/acmat(range_row+1,range_column+1);

end % of function