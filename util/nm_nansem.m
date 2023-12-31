function y = nm_nansem(x,dim)
% FORMAT: Y = NANSEM(X,DIM)
% 
%    Standard error of the mean ignoring NaNs
%
%    NANSTD(X,DIM) calculates the standard error of the mean along any
%    dimension of the N-D array X ignoring NaNs.  
%
%    If DIM is omitted NANSTD calculates the standard deviation along first
%    non-singleton dimension of X.
%
%    Similar functions exist: NANMEAN, NANSTD, NANMEDIAN, NANMIN, NANMAX, and
%    NANSUM which are all part of the NaN-suite.

% -------------------------------------------------------------------------
%    author:      Jan Gl�scher
%    affiliation: Neuroimage Nord, University of Hamburg, Germany
%    email:       glaescher@uke.uni-hamburg.de
%    
%    $Revision: 1.1 $ $Date: 2004/07/22 09:02:27 $

if isempty(x)
	y = NaN;
	return
end

if nargin < 2
	dim = min(find(size(x)~=1));
	if isempty(dim)
		dim = 1; 
	end	  
end

% Find NaNs in x and nm_nanmean(x)
nans = isnan(x);
fullnans = sum(nans,dim) == size(x,2);

count = size(x,dim) - sum(nans,dim);

% Protect against a  all NaNs in one dimension
i = find(count==0);
count(i) = 1;

y = nm_nanstd(x,dim)./sqrt(count);

y(i) = i + NaN;
y(fullnans) = nan;

% $Id: nm_nansem.m,v 1.1 2004/07/22 09:02:27 glaescher Exp glaescher $
