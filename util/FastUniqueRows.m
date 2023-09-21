function idx = FastUniqueRows(A)
% Fast Unique Rows
% ----------------------------------------------------------------
% This MATLAB function provides a faster version of MATLAB's unique
% rows method (i.e. 'unique(points,''rows'')').
%
% Sample Usage:
% points=[1 1; 2 3; 0 1; 0 0; 0 0; 2 3; 3 5; 1 1; 0 0; 2 3; 3 6; 3 5]
% [IA2]=myFastUniqueRows(points); %equivalent is '[~,IA1]=unique(points,'rows');'
% indices=sort(IA2);
% points_Unique=points(indices,:);
%
% Input:
%   A: Rows of n-Dimensional Vectors
% Output:
%   idx: Row indices containing unique vectors (not sorted)
%
% Please see execFUR.m and execFUR_PT.m scripts for comparison, unit tests
% and automated performance evaluations.
%
% Author: Caglar Arslan
% Date  : May 2020
% License : 'No License'
% ----------------------------------------------------------------
[A_sorted, idx1] = sortrows(A);
k    = find([true; any(diff(A_sorted, 1, 1), 2); true]);
idx2 = k(diff(k) >= 1);
idx=idx1(idx2);
% B    = A(idx1(idx2), :); %This line is left for experimental and educational purposes.
end
%% Caglar Arslan, 2020, No License