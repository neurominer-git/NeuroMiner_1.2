function hi = plotshaded(x, y, fstr, mrksize, linsize, rotateflag)
% x: x coordinates
% y: either just one y vector, or 2xN or 3xN matrix of y-data
% fstr: format ('r' or 'b--' etc)
%
% example
% x=[-10:.1:10];plotshaded(x,[sin(x.*1.1)+1;sin(x*.9)-1],'r');
 
if ~exist('mrksize', 'var') || isempty(mrksize), mrksize = 10; end
if ~exist('linsize', 'var') || isempty(linsize), linsize = 1; end
if ~exist('rotateflag', 'var') || isempty(rotateflag), rotateflag = 0; end
if ~mrksize
    MarkerType = 'none';
    mrksize=0.1;
else
    MarkerType = 'o';
end

if size(y,1)>size(y,2)
    y=y';
end
 
if size(y,1)==1 % just plot one line
    plot(x,y,fstr);
 
elseif size(y,1)==2 %plot shaded area
    px=[x,fliplr(x)]; % make closed patch
    py=[y(1,:), fliplr(y(2,:))];
    hx = patch(px,py,1,'FaceColor',fstr,'EdgeColor','none','FaceAlpha',0.25);
    hi=hx;
 
elseif size(y,1)==3 % also draw mean
    px=[x,fliplr(x)];
    py=[y(1,:), fliplr(y(3,:))];
    if rotateflag
        hi = plot(y(2,:),x,'Color', fstr, 'Marker', MarkerType, 'MarkerSize', mrksize, 'MarkerEdgeColor','w','MarkerFaceColor', fstr,'LineWidth', linsize);
        hx = patch(py,px,1,'FaceColor',fstr,'EdgeColor','none','FaceAlpha',0.15);
    else
        hi = plot(x, y(2,:),'Color', fstr, 'Marker', MarkerType, 'MarkerSize', mrksize, 'MarkerEdgeColor','w','MarkerFaceColor', fstr,'LineWidth', linsize);
        hx = patch(px, py,1,'FaceColor',fstr,'EdgeColor','none','FaceAlpha',0.15);
   end
end
