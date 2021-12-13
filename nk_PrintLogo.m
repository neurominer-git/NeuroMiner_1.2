function nk_PrintLogo(maindlg)
global NMinfo EXPERT SPMAVAIL

clc
if ~SPMAVAIL
    %cl = [0.1,0.5,0]; 
    mode = 'non-imaging mode';
else
    %cl = NMinfo.cllogo; 
    mode = [];
end
if exist('maindlg','var') && maindlg
    fprintf('\n\t*******************************************')
    fprintf('\n\t****\\                                 /****')
    fprintf('\n\t*****\\     ');fprintf('~~~~~~~~~~~~~~~~~~~~~ '); fprintf('    /*****');
    fprintf('\n\t******\\     ');fprintf('N E U R O M I N E R ');  fprintf('    /******');
    fprintf('\n\t*******\\   ');fprintf('~~~~~~~~~~~~~~~~~~~~~ '); fprintf('  /*******');
    fprintf('\n\t*******/                           \\*******')
    fprintf('\n\t******/     ');fprintf('pattern recognition ');   fprintf('    \\******');
    fprintf('\n\t*****/      ');fprintf('for neurodiagnostic ');   fprintf('     \\*****');
    fprintf('\n\t****/           ');fprintf('applications ');      fprintf('         \\****');
    fprintf('\n\t***/                                   \\***')
    fprintf('\n\t*******************************************')
    if ~isempty(mode),fprintf('\n\t%s',mode);end
else
    fprintf('\t~~~~~~~~~~~~~~~~~~~~~~~ \n');
    fprintf('\t  N E U R O M I N E R \n');
    fprintf('\t~~~~~~~~~~~~~~~~~~~~~~~ ');
end
if EXPERT
    fprintf('\n')
    fprintf('\t>>> EXPERT MODE <<< ')
end
fprintf('\n\t%s \n', NMinfo.info.ver); fprintf('\n')
fprintf('(c) %s, %s ', NMinfo.info.author, NMinfo.info.datever)
%disp(['<a href="matlab: sendmail(''' NMinfo.info.email ''',''NeuroMiner'')">' NMinfo.info.email '</a>'])
fprintf('\n    nm@pronia.eu \n')