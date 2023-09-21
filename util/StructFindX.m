function FoundFields = StructFind(search_struct,search_object,varargin)

% Function to search for entries within a structure. This is done
% recursively, running through all elements of array-structures.
%
% Input:
%
%       search_struct:      structure to be searched, could also be array of
%                           struct
%
%       search_object:      string, integer, cell, array or other thing to
%                           be searched for
%
%    [optional]
%
%       structure_name:     name of the structure to be searched. Used for
%                           full output of the structure content
%
% Output:
%
%       FoundFields:        Cell array of fields where the search object
%                           is found.
%
% 
% Alexander Mering
% Version 1.0; 09. February 2012
%

if ~isstruct(search_struct)
   error('Input should be a structure!')
   return
end

FoundFields = Really_StructFind(search_struct,search_object);

% Prepend structure name to the output from the search function
if nargin > 2 && length(varargin{1}) > 0
  FoundFields = strcat(repmat({varargin{1}},length(FoundFields),1),FoundFields');
end


%% Search function
% used for recursive search through the structure
function FoundFieldsList = Really_StructFind(search_struct,search_object)

% initialize output
FoundFieldsList = cell('');

% get fieldnames
struct_fields  = fieldnames(search_struct);

% outer loop is array of struct
for n = 1:length(search_struct)
    
    % inner loop runs through the fields of the struct
    for m= 1: length(struct_fields)
        
        % get field to be worked with
        working_field = search_struct(n).(struct_fields{m});
       
        if isstruct(working_field)
            % run search routine recursively
            InsideFound = Really_StructFind(working_field,search_object);
            
            % append search results to current output
            for k = 1:length(InsideFound)
                FoundFieldsList{end+1} = strcat('(',num2str(n),').',struct_fields{m},InsideFound{k});
            end
            
        elseif isequal(working_field,search_object)   % HIT!!
            % append found fields to the result cell
            FoundFieldsList{end+1} = strcat('(',num2str(n),').',struct_fields{m});
        end
        
    end
end
