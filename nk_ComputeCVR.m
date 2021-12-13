function [CVR, SE] = nk_ComputeCVR(I, SEL)

SUM2 = nm_nansum(I.^2,2);
SUM = nm_nansum(I,2);
SQ = sqrt(SUM2 ./ SEL - (SUM ./ SEL).^2);
SE = SQ./sqrt(SEL)*1.96;

% Mean relevance/weight vector across partitions
MN = nm_nansum(SUM, 2) ./ SEL;

% Compute CVR
CVR = MN./SE;