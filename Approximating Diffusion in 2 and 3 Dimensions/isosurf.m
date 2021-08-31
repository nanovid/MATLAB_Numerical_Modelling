% isosurf( U, n )
%
% Plot isosurfaces for solutions to various equations
% in three dimensions.
%
%  Author:  Douglas Wilhelm Harder
%  Copyright (c) 2010 by Douglas Wilhelm Harder.  All rights reserved.
%
% The input is an m1 x m2 x m3 array A which holds
% values within the region and plots n + 1 isosurfaces
% interpolated from the smallest to the largest value.
%
function [] = isosurf( U, n )
    % Find the minimum and maximum values in the plot
    lo = min( min( min( U ) ) );
    hi = max( max( max( U ) ) );
    vals = linspace( lo, hi, n );

    newplot;

    for i = 1:n
        Usurface = patch(                                             ...
            isosurface( U, vals(i) ),                                 ...
            'FaceColor', [(i - 1)/(n - 1), 0, (1 - (i - 1)/(n - 1))], ...
            'EdgeColor', 'none'                                       ...
        );

        patch(                         ...
            isocaps( U, vals(i) ), ...
            'FaceColor', 'interp',     ...
            'EdgeColor', 'none'        ...
        );

        isonormals( U, Usurface );
    end

    view(3);
    camlight('headlight');
    lighting phong;
    alpha(0.1);
    axis equal;
end
