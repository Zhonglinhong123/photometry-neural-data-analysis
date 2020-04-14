%Makes one row of data (one trial) for eventual PSTH array

function [PsthRow] = processPhotDataRow_normDat(normDat, thisIndex, nTsPrev, nTsPost)

if nTsPrev < thisIndex && size(normDat,1) > (thisIndex + nTsPost);
    PsthRow = (normDat((thisIndex - nTsPrev):(thisIndex + nTsPost)))';
elseif nTsPrev >= thisIndex && size(normDat,1) > (thisIndex + nTsPost);
    mismatch = (nTsPrev - thisIndex)+1;
    PsthRow1 = NaN (1, mismatch);
    PsthRow2 = (normDat((1):(thisIndex + nTsPost)))';
    PsthRow = [PsthRow1, PsthRow2];
else
    mismatch = ((thisIndex + nTsPost) - size(normDat,1));
    PsthRow1 = (normDat((thisIndex - nTsPrev):(size(normDat,1))))';
    PsthRow2 = NaN (1, mismatch);
    PsthRow = [PsthRow1, PsthRow2];
end
