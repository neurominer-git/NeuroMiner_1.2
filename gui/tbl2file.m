function [ERR, STATUS, fil, typ] = tbl2file(tbl, filename, sheetname )
% Write table data to either a Excel or a text-based file
ERR=[]; STATUS = 0;
[pth,filename,ext] = fileparts(filename); 
filename = regexprep(filename,' ','_');
fil = fullfile(pth, [filename ext]);
try
    if ispc 
        typ='xls';  
        if isempty(ext)
            fil = sprintf('%s.%s',filename,typ);
        end
        s_rownames = size(tbl.rownames);
        if s_rownames(1)==1
            tbl.rownames = tbl.rownames';
        end
        T = [tbl.rownames array2table(tbl.array)];
        T.Properties.VariableNames = tbl.colnames;
        if numel(sheetname)>31
            sheetname = inputdlg(sprintf('Sheet name %s is too long. Provide alternative sheet name:',sheetname),'Sheetname problem',1,{sheetname},'on');
        end
        writetable(T, fil, 'Sheet', char(sheetname));                
        STATUS = 1;
    else
        typ='csv';
        if isempty(ext)
            fil = sprintf('%s_%s.%s', filename, sheetname, typ);
        end
        fid = fopen(fil,'w');
        %Print header row
        for i=1:size(tbl.colnames,2)
            fprintf(fid,'%s',tbl.colnames{i});
            if i<size(tbl.colnames,2),fprintf(fid,'\t'); end
        end
        %Print table rows
        for i=1:numel(tbl.rownames)
            fprintf(fid,'\n%s', tbl.rownames{i});
            if iscell(tbl.array)
                for j=1:size(tbl.array,2)
                    fprintf(fid,'\t%g', tbl.array{i,j});
                end
            else
                for j=1:size(tbl.array,2)
                    fprintf(fid,'\t%g', tbl.array(i,j));
                end
            end
        end
        fclose(fid);
        STATUS = 1;
    end
    fprintf('\n%s successfully written to disk.\n', fil );
catch ERR
    errordlg(ERR.message);
end
   