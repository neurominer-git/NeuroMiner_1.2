% create mex files all 


% %% list all c files 
% cd 'trainpredict'
% trainpred_c = dir('**/*.c');
% trainpred_cpp = dir('**/*.cpp');
% 
% cd 'preproc'
% preproc_c = dir('**/*.c');
% preproc_cpp = dir('**/*.cpp');

% 
% %% trainpredict c
% error_count_trainpredict_c = 0;
% for i= 1:size(trainpred_c,1)
%     cur_c = trainpred_c(i);
%     %fullpath = sprintf('%s/%s', cur_c.folder, cur_c.name);
%     %disp(fullpath)
%     
%     cd(cur_c.folder); 
%     try 
%         eval(sprintf('mex CFLAGS=''$CFLAGS -std=c99'' -largeArrayDims %s', cur_c.name));
%     catch
%         disp(cur_c.name)
%         error_count_trainpredict_c = error_count_trainpredict_c + 1;
%     end
% 
% end
% 
% %% trainpredict cpp
% error_count_trainpredict_cpp = 0;
% for i= 1:size(trainpred_cpp,1)
%     cur_cpp = trainpred_cpp(i);
%     %fullpath = sprintf('%s/%s', cur_cpp.folder, cur_cpp.name);
%     %disp(fullpath)
%     
%     cd(cur_cpp.folder); 
%     try 
%         eval(sprintf('mex CFLAGS=''$CFLAGS -std=c99'' -largeArrayDims %s', cur_cpp.name));
%     catch
%         disp(cur_cpp.name)
%         error_count_trainpredict_cpp = error_count_trainpredict_cpp + 1;
%     end
% end
% %% preproc c
% error_count_preproc_c = 0;
% for i= 1:size(preproc_c,1)
%     cur_c = preproc_c(i);
%     %fullpath = sprintf('%s/%s', cur_c.folder, cur_c.name);
%     %disp(fullpath)
%     
%     cd(cur_c.folder); 
%     try 
%         eval(sprintf('mex CFLAGS=''$CFLAGS -std=c99'' -largeArrayDims %s', cur_c.name));
%     catch
%         disp(cur_c.name)
%         error_count_preproc_c = error_count_preproc_c + 1;
%     end
% 
% end
% 
% %% preproc cpp
% error_count_preproc_cpp = 0;
% for i= 1:size(preproc_cpp,1)
%     cur_cpp = preproc_cpp(i);
%     %fullpath = sprintf('%s/%s', cur_cpp.folder, cur_cpp.name);
%     %disp(fullpath)
%     
%     cd(cur_cpp.folder); 
%     try 
%         eval(sprintf('mex CFLAGS=''$CFLAGS -std=c99'' -largeArrayDims %s', cur_cpp.name));
%     catch
%         disp(cur_cpp.name)
%         error_count_preproc_cpp = error_count_preproc_cpp + 1;
%     end
% 
% end

%% list all c-files and cpp-files in NM folder
% cd NeuroMiner_Current
nm_cfiles = dir('**/*.c')
nm_cppfiles = dir('**/*.cpp')

%%
error_count_nm_cfiles = 0;
errors_nm_cfiles = {};
for i= 1:size(nm_cfiles,1)
    cur_c = nm_cfiles(i);
    fullpath = sprintf('%s/%s', cur_c.folder, cur_c.name);
    %disp(fullpath)
    
    cd(cur_c.folder); 
    try 
        eval(sprintf('mex CFLAGS=''$CFLAGS -std=c99'' -largeArrayDims %s', cur_c.name));
    catch
        disp(cur_c.name)
        error_count_nm_cfiles = error_count_nm_cfiles + 1;
        errors_nm_cfiles{error_count_nm_cfiles} = fullpath;
    end

end
%%
error_count_nm_cppfiles = 0;
errors_nm_cppfiles = {};
for i= 1:size(nm_cppfiles,1)
    cur_cpp = nm_cppfiles(i);
    fullpath = sprintf('%s/%s', cur_cpp.folder, cur_cpp.name);
    %disp(fullpath)
    
    cd(cur_cpp.folder); 
    try 
        eval(sprintf('mex CFLAGS=''$CFLAGS -std=c99'' -largeArrayDims %s', cur_cpp.name));
    catch
        disp(cur_cpp.name)
        error_count_nm_cppfiles = error_count_nm_cppfiles + 1;
        errors_nm_cppfiles{error_count_nm_cppfiles} = fullpath;
    end

end
