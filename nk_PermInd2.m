% =========================================================================
% FORMAT permmat = nk_PermInd(nPerms,Labels,Constraint)
% =========================================================================
% 
% Generate 'nPerms' within-group permutations of class membership indices
% 'Labels' and guarantee that replicated permutations do not occur
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Neurominer, (c) Nikolaos Koutsouleris 1/2012
function permmat = nk_PermInd2(nPerms, Labels, Constraint)

uLabels = unique(Labels);

if any(~isfinite(uLabels)) 
    uLabels(~isfinite(uLabels))=[];
    uLabels(end+1)=NaN;
end
mL = numel(uLabels);

ll = size(Labels,1);
trylim = 500;
permmat = zeros(nPerms,ll);

if ~exist('Constraint','var') || isempty(Constraint)
    mC = 1; uConstraint = []; cfl = false;
else
    uConstraint = unique(Constraint);
    mC = numel(uConstraint);  cfl = true;
end

for i=1:mL % Loop through classes
    
   for h=1:mC % Eventually loop through constraint vector
        
       if cfl % Permute within class and within grouping 
           indcl = find(Labels == uLabels(i) & Constraint == uConstraint(h));
       else % Permute within class
           if ~isfinite(uLabels(i))
               indcl = find(~isfinite(Labels));
           else
               indcl = find(Labels == uLabels(i));
           end
       end    
  
       if isempty(indcl) 
           indcl = (1:numel(Labels))'; 
       end
       
       lindcl = length(indcl);
       j=1; 

       while j <= nPerms % Generate permutations
        
           rclass = randperm(lindcl);
           rindcl = indcl(rclass)';
           permmat(j,indcl) = rindcl;

           % Check whether permutation occured before
           if j > 1
               if ~sum(ismember(permmat(1:j-1,indcl),rindcl,'rows'))
                   j = j + 1;    
               else
                   fprintf('\nReplicated permutation indices detected. Skip current permutation.')
                   trycnt = trycnt + 1;
                   if trycnt == trylim
                       if cfl
                           error(sprintf(['\nTry limit reached. Permutation of class membership data failed!', ...
                                          '\n\tConstrain index = %g', ...
                                          '\n\tNo. of observations in current class and constrain index = %g\n'], ...
                                          uConstraint(h), numel(rindcl)));
                       else
                           error('Try limit reached. Permutation of class membership data failed.')
                       end
                   end
               end
           else
               trycnt = 1;
               j = j +1;
           end
       end
   end
end