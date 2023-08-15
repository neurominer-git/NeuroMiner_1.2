% =========================================================================
% =                  MULTI-GROUP CLASSIFICATION PLOTS                     =
% =========================================================================
function handles = display_multiclassplot(handles)

axes(handles.axes1);
%uistack(handles.axes1,'top')
cla;
hold on

h_onevsall_val  = get(handles.selOneVsAll_Info,'Value');
probabilities   = handles.MultiClass.probabilities;
errors          = logical(handles.MultiClass.errors');
[m,n]           = size(probabilities);
subjvec         = (1:m)';
errx            = find(errors == 1);
erry            = zeros(size(errx));
[~,maxI]        = max(probabilities,[],2);
for i=1:length(errx)
    erry(i) = probabilities(errx(i),maxI(errx(i)));
end
ind0 = false(m,n); ind0(sub2ind([m n], 1:m, maxI')) = true;
if m>1000
    MrkSizeP = 6;
    MrkSizeL = 2;
elseif m>600
    MrkSizeP = 8;
    MrkSizeL = 3;
elseif m>300
    MrkSizeP = 10;
    MrkSizeL = 4;
else
    MrkSizeP = 12;
    MrkSizeL = 6;
end
offy=0;
lg_groups = [];
lg_handles = [];

for h=1:n
    
    % Index to current group
    indh = maxI == h;
    
    if ~isempty(indh)
        
        % Check if one-vs-all or one-vs-one
        % if former use grey color for the 'all' group
        if h_onevsall_val > 1 && h ~= h_onevsall_val-1
            hcolptin = [0.7 0.7 0.7];
        else
            % Otherwise use group 2 color
            hcolptin = handles.colptin(h,:);
        end
        
        % Paint background in light color of current group
        lxL = find(handles.MultiClass.labels==h);
        v = [ [ lxL(1) offy ]; [ lxL(1) 1 ]; [ lxL(end) 1 ]; [ lxL(end) offy ] ];             
        patch('Faces',[1 2 3 4], 'Vertices', v, 'FaceColor', hcolptin, 'EdgeColor', 'none', 'FaceAlpha', 0.1)

        handles.MultiClass.bx(h)   = plot(subjvec,probabilities(:,h), ...
                                        handles.colpt, ...
                                        'MarkerSize', MrkSizeL, ...
                                        'MarkerFaceColor',hcolptin,...
                                        'Color',hcolptin, ...
                                        'LineWidth',1, ...
                                        'LineStyle','none');
        lg_groups = [lg_groups {sprintf('%s: probabilities',handles.NM.groupnames{h})} ];
        lg_handles = [lg_handles handles.MultiClass.bx(h)];
        indhx =  indh & ~errors; indh0 = indh & errors;

        if sum(indhx)>0
            handles.MultiClass.b(h)    = plot(subjvec(indhx),probabilities(indhx,h), ...
                                        handles.colpt, ...
                                        'MarkerSize',MrkSizeP,...
                                        'MarkerFaceColor',hcolptin,...
                                        'LineWidth',handles.DataMarkerWidth,...
                                        'Color',hcolptin, ...
                                        'LineStyle','none');
            lg_groups = [lg_groups {sprintf('%s: correct classifications',handles.NM.groupnames{h})} ];
            lg_handles = [lg_handles handles.MultiClass.b(h)];
        else
            handles.MultiClass.b(h) = plot(0,0);
        end
        if sum(indh0)>0
            handles.MultiClass.be(h) = plot(subjvec( indh0 ),probabilities(indh0,h), ...
                                        handles.colpx, ...
                                        'MarkerSize',MrkSizeP,...
                                        'MarkerFaceColor',hcolptin,...
                                        'Color',hcolptin, ...
                                        'LineWidth',handles.DataMissMarkerWidth,...
                                        'LineStyle','none');
            lg_groups = [lg_groups {sprintf('%s: misclassifications',handles.NM.groupnames{h})} ];
            lg_handles = [lg_handles handles.MultiClass.be(h)];
        else
            handles.MultiClass.be(h) = plot(0,0);
        end
    else
        handles.MultiClass.b(h)     = plot (0,0);
        handles.MultiClass.bx(h)    = plot (0,0);
    end
end

xlim([min(subjvec) max(subjvec)]);
ylim([0 1]);

% Define textbox info data 
pss = cell(1,m); subjects = handles.MultiClass.cases;
figdata.y = zeros(m,1);
uL = unique(handles.MultiClass.labels);
for i=1:numel(pss)
    iI = uL == handles.MultiClass.labels(i);
    expgroupi   = handles.NM.groupnames{iI}; 
    predgroupi  = handles.NM.groupnames{maxI(i)}; 
    pss{i} = sprintf(['Subject ID: %s' ...
                    '\nExpected group: %s' ...
                    '\nPredicted Group: %s'], subjects{i}, expgroupi, predgroupi);
    for j=1:handles.ngroups
        pss{i} = sprintf('%s\nProbability %s: %1.2f', pss{i}, handles.NM.groupnames{j}, probabilities(i,j));
    end
    figdata.y(i) = probabilities(i,maxI(i));
end
hText = uicontrol('Style','text','String', pss{numel(pss)},'FontSize',11, 'Units','normalized', 'Parent', gcf,'Visible','off');
figdata.cases       = handles.MultiClass.cases;
figdata.x           = subjvec;
figdata.patterntext = pss;
figdata.parentui    = handles.pnBinary;
figdata.pnpos       = handles.pnBinary.Position;
figdata.figpos      = handles.figure1.Position;
figdata.hPanel      = uipanel('Units','norm', 'Position',hText.Extent, 'BorderType','etchedout', 'BackgroundColor', [.6 .7 .6], 'Visible','off');
figdata.textHdl     = annotation(figdata.hPanel, 'textbox', 'String','', ...
                            'Interpreter','none', ... %'VerticalAlign', 'Top', ...
                            'Color', 'black', ...
                            'BackgroundColor',[.6 .7 .6], ...
                            'Position', [0 0 0.99 0.99], ...
                            'EdgeColor', [.6 .7 .6], ...
                            'LineWidth', 0.1, ...
                            'Margin', 5, ...
                            'FitBoxToText','on', ...
                            'Visible','off');
set(handles.axes1,'UserData',figdata);
ylabel('Multi-class probabilities'); 

xLimitsVec = linspace(1,numel(handles.MultiClass.labels), numel(handles.axes1.XAxis.TickValues)-1);
zeroline = ones(1,numel(xLimitsVec))*0.5;

plot(xLimitsVec,zeroline,'LineWidth',handles.ZeroLineWidth,'Color',rgb('Grey'))

if h_onevsall_val == 1 
    legend(lg_handles, lg_groups, 'Location','Best', 'FontSize',handles.LegendFontSize);
else
    legend('off')
end