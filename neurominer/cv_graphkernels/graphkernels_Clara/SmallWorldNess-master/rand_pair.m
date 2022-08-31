function y = rand_pair(a,b,n,varargin)

% RAND_PAIR Uniformly-sampled random pair combinations
%
%   RAND_PAIR(A,B,N) is a 2-column matrix of length N, where each row is a unique pair of numbers
%   chosen at random from the arrays A and B.
%
%   RAND_PAIR(A,B,N,S) starts the pseudo-random number generator using seed S. 
%
%   RAND_PAIR(A,B,N,S,'bi') ensures that each pair is unique in both directions 
%   (i.e. if [2,3] already chosen then [3,2] rejected). If random seed not required then
%   use the empty matrix [] instead of S.
%
%   RAND_PAIR(....,'diff') ensures that no pair is chosen which consists of the same numbers (i.e. [2,2] rejected);
%   set S and 'bi' flag to [] if not required. 
%
%   Mark Humphries 15/10/2004

if nargin > 3 & ~isempty(varargin{1})
    rand('state',varargin{1});
end

if n > length(a) * length(b)
    error('Cannot produce that many permutations of pairs')
end

if nargin == 6 & strcmp('diff',varargin{3})
    diff = 1;
else
    diff = 0;
end
    
y = zeros(n,2);

if nargin >= 5 & strcmp('bi',varargin{2})
    % bidirectional uniqueness - generate pair first then check
    for loop = 1:n
        isunique = 0;
        while ~isunique
            first_idx = ceil(rand * length(a));           
            second_idx = ceil(rand * length(b));
            
            % check for uniqueness if either (a) not checking for sameness OR (b) are checking, and not the same
            if ~diff | (diff & a(first_idx) ~= b(second_idx))
                % test both combinations
                n1 = find(y(:,1) == a(first_idx));
                exists = find(b(second_idx) == y(n1,2));    % is second number in second column at same row number?
                
                if exists
                    % this pair exists
                    isunique = 0;
                else
                    % otherwise test reversed pair
                    n2 = find(y(:,2) == a(first_idx));
                    exists = find(b(second_idx) == y(n2,1));    % is second number in *first* column at same row number?
                    if exists
                        isunique = 0;
                    else
                        % pair does not exist in either configuration
                        isunique = 1;
                    end
                end                    
            end
               
        end
        y(loop,:) = [a(first_idx) b(second_idx)];
	end    
elseif diff
    % sameness checking without bidirection checking...
    
    % check for fast version if pairs drawn from same distribution and
    % there are many links to find
    if length(a)==length(b) & a == b % & n > 8000
        % all possible pair numbers
        possible_pair_nos = 1:length(a)^2;

        % pair numbers which are duplicates
        duplicates = 1:length(a)+1:length(a)^2;
        
        % remove duplicates
        possible_pair_nos(duplicates) = [];
        
        % create random permutation sequence
        seq = randperm(length(possible_pair_nos));
        
        % select first n of that sequence to create pairs
        selected_pairs = possible_pair_nos(seq(1:n));
        
        first_idx = ceil(selected_pairs ./ length(b));
        second_idx = rem(selected_pairs,length(b));
        second_idx(second_idx==0) = length(b);
        y = [a(first_idx)',b(second_idx)']; 
        
    else
        old_pairs = zeros(n,1);
	
        % number of unique permutations
        perms = length(a) * length(b);
        
        % fastish version for no reverse duplicate checking
		for loop = 1:n
            isunique = 0;
            while ~isunique
                pair_no = ceil(rand * perms);              % assumes pairs generated thus [a1,b1; a1,b2; a1,b3; a2,b1; a2,b2....]
                if all(old_pairs ~= pair_no)
                    isunique = 1;
                       
                    first_idx = ceil(pair_no / length(b));      % index for first column
                    second_idx = rem(pair_no,length(b)); 
                    if second_idx == 0 second_idx = length(b); end
                
                    first_num = a(first_idx);
                    second_num = b(second_idx);    % number for second column
                    
                    if first_num == second_num      % if same then no longer unique
                        is_unique = 0;
                    end
                end
            end
            old_pairs(loop) = pair_no;
            y(loop,:) = [first_num second_num];
        end
    end % version if..then...
else
    old_pairs = zeros(n,1);

    % number of unique permutations
    perms = length(a) * length(b);
    
    % fast version for no reverse duplicate or sameness checking
	for loop = 1:n
        isunique = 0;
        while ~isunique
            pair_no = ceil(rand * perms);              % assumes pairs generated thus [a1,b1; a1,b2; a1,b3; a2,b1; a2,b2....]
            if all(old_pairs ~= pair_no)
                isunique = 1;
            end            
        end
        old_pairs(loop) = pair_no;
        
        first_idx = ceil(pair_no / length(b));      % number for first column
        second_idx = rem(pair_no,length(b)); 
        if second_idx == 0 second_idx = length(b); end
        
        first_num = a(first_idx);
        second_num = b(second_idx);    % number for second column
        y(loop,:) = [first_num second_num];
	end
end



