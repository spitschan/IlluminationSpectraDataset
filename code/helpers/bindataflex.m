function [xx, binnedVals, minVal, maxVal, meanVal, medianVal, sdVal] = bindataflex(x, y, xx)
% [xx, binnedVals, minVal, maxVal, meanVal, medianVal, sdVal] = bindataflex(x, y, xx)
% 
% Bins data.
%
% 3/26/2016 	spitschan	 Added to repository.

[~,ind] = histc(x, xx);
uniqueInds = unique(ind);
for m = 1:length(xx)
    idx = find(xx(ind) == xx(m));
    if isempty(idx)
        minVal(m) = NaN;
        maxVal(m) = NaN;
        meanVal(m) = NaN;
        medianVal(m) = NaN;
        sdVal(m) = NaN;
        binnedVals{m} = [];
    else
        minVal(m) = min(y(idx));
        maxVal(m) = max(y(idx));
        meanVal(m) = nanmean(y(idx));
        medianVal(m) = nanmedian(y(idx));
        sdVal(m) = nanstd(y(idx));
        binnedVals{m} = y(idx);
    end
end