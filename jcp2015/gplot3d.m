function [Xout,Yout]=gplot3d(A,xy,lc)
%GPLOT Plot graph, as in "graph theory".
%   GPLOT(A,xy) plots the graph specified by A and xy. A graph, G, is
%   a set of nodes numbered from 1 to n, and a set of connections, or
%   edges, between them.  
%
%   In order to plot G, two matrices are needed. The adjacency matrix,
%   A, has a(i,j) nonzero if and only if node i is connected to node
%   j.  The coordinates array, xy, is an n-by-2 matrix with the
%   position for node i in the i-th row, xy(i,:) = [x(i) y(i)].
%   
%   GPLOT(A,xy,LineSpec) uses line type and color specified in the
%   string LineSpec. See PLOT for possibilities.
%
%   [X,Y] = GPLOT(A,xy) returns the NaN-punctuated vectors
%   X and Y without actually generating a plot. These vectors
%   can be used to generate the plot at a later time if desired.  As a
%   result, the two argument output case is only valid when xy is of type
%   single or double.
%   
%   See also SPY, TREEPLOT.

%   John Gilbert
%   Copyright 1984-2013 The MathWorks, Inc. 

[i,j] = find(A);
[~, p] = sort(max(i,j));
i = i(p);
j = j(p);

X = [ xy(i,1) xy(j,1)]';
Y = [ xy(i,2) xy(j,2)]';
Z = [ xy(i,3) xy(j,3)]';

if isfloat(xy) || nargout ~= 0
    X = [X; NaN(size(i))'];
    Y = [Y; NaN(size(i))'];
    Z = [Z; NaN(size(i))'];
end

if nargout == 0
    if ~isfloat(xy)
        if nargin < 3
            lc = '';
        end
        [lsty, csty, msty] = gplotGetRightLineStyle(gca,lc);    
        plot3(X,Z,Y,'LineStyle',lsty,'Color',csty,'Marker',msty);
    else
        if nargin < 3
            plot3(X(:),Z(:),Y(:));
        else
            plot3(X(:),Z(:),Y(:),lc);
        end
    end
else
    Xout = X(:);
    Yout = Y(:);
    Zout = Z(:);
end

function [lsty, csty, msty] = gplotGetRightLineStyle(ax, lc)
%  gplotGetRightLineStyle
%    Helper function which correctly sets the color, line style, and marker
%    style when plotting the data above.  This style makes sure that the
%    plot is as conformant as possible to gplot from previous versions of
%    MATLAB, even when the coordinates array is not a floating point type.
co = get(ax,'ColorOrder');
lo = get(ax,'LineStyleOrder');
holdstyle = getappdata(gca,'PlotHoldStyle');
if isempty(holdstyle)
    holdstyle = 0;
end
lind = getappdata(gca,'PlotLineStyleIndex');
if isempty(lind) || holdstyle ~= 1
    lind = 1;
end
cind = getappdata(gca,'PlotColorIndex');
if isempty(cind) || holdstyle ~= 1
    cind = 1;
end
nlsty = lo(lind);
ncsty = co(cind,:);
nmsty = 'none';
%  Get the linespec requested by the user.
[lsty,csty,msty] = colstyle(lc);
if isempty(lsty)
    lsty = nlsty;
end
if isempty(csty)
    csty = ncsty;
end
if isempty(msty)
    msty = nmsty;
end

