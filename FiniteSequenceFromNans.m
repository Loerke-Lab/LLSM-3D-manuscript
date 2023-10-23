function [VectorOut, idx_out] = FiniteSequenceFromNans(VectorIn)
% This function takes a vector, which may contain NaN values throughout,
% and outputs the longest continuous sequence which contains no NaNs. It 
% also outputs the indices of that sequence. For example, the input 
% [ NaN NaN 1 4 NaN 7 9 12 NaN] gives the out put [7 9 12] with indices
% [6 7 8].

if all(isnan(VectorIn(:)))  % These 5 lines were added to return empty outputs
    VectorOut = [];         % when the input vector is all nans
    idx_out = [];
    return;
end


finiteIdx = find(isfinite(VectorIn));
x = [0 cumsum(diff(finiteIdx)~=1)];
idx_out = finiteIdx(x==mode(x));
VectorOut = VectorIn(idx_out);
    
end % end FiniteSequenceFromNans