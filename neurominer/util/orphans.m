function[list] = orphans(varargin)
% Find m-files whose name is not found in other m-files within a directory tree
% Outputs a cell array of strings containing file names
% Notes : Calls MGREP and RDIR, available from File Exchange
% Example : list = orphans('C:\Matlab');  char(list)
% Dimitri Shvorob, dimitri.shvorob@gmail.com, 7/23/08, with improvements by Ryan Hamilton
if nargin == 1
    cpath = varargin{1};
else
    cpath = cd;                                  %#ok
end
if exist(cpath,'dir') ~=7
   error('??? Path %s not found.',cpath)
end
f = rdir([cpath filesep '**' filesep '*.m']);
n = length(f);
c = cell(n,1);
file = cell(n,1);
for i = 1:n
    s = f(i).name;
    j = max(strfind(s,filesep))+ 1;
    c{i} = s(j:end-2);
    file{i} = s;
end
o = {};
for i = 1:n
    s = mgrep(c{i},cpath,'recurse','on','showline','off');
    s(strcmp(f(i).name,s)) = [];
    if isempty(s)
       o{end+1} = file{i};                       %#ok
    end
end
list = o(:);