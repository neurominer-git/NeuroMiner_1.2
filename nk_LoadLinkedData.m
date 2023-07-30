function datacontainer = nk_LoadLinkedData(datacontainer, varind)
global VERBOSE 
if iscell(datacontainer.Y)
    for i=1:numel(varind)
        if ischar(datacontainer.Y{varind(i)}) 
            if exist(datacontainer.Y{i},"file")
                fprintf('\nLoading data file %s into modality #%g', datacontainer.Y{varind(i)}, varind(i))
                load(datacontainer.Y{varind(i)});
                datacontainer.Y{varind(i)} = Y;
            else
                error('\nData file %s not found. Please update link using NM Modality Manager!', datacontainer.Y{varind(i)});
            end
        else
            if VERBOSE, fprintf('\nModality #%g contains already matrix data', varind(i)); end
        end
    end
% For compatibility with NM 1.1
elseif ischar(datacontainer.Y) 
   if ~exist(datacontainer.Y,'file')
       error('Linked data container not found in expected path: %s',datacontainer.Y);
   end
   fprintf('\nLoading linked data container: %s', datacontainer.Y);
   load(datacontainer.Y);
   if exist("OOCV","var")
       datacontainer = OOCV;
   else
       error('\nOOCV structure not found.')
   end
end