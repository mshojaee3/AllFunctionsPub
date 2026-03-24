function [Rx,Ry] = computeReactionFromNodalCSV(file, edgeName)
% computeReactionFromNodalCSV
% Approximate boundary reaction force from nodal Cauchy stress data.
%
% INPUT
%   file     : CSV file path
%   edgeName : 'right', 'left', 'top', 'bottom'
%
% OUTPUT
%   Rx, Ry   : reaction force components
%
% REQUIRED COLUMNS
%   X, Y, U_U1, U_U2, S11, S22, S12
%
% METHOD
%   1) Extract nodes on requested boundary in reference coordinates
%   2) Build deformed coordinates x = X+U_U1, y = Y+U_U2
%   3) Sort nodes along the edge
%   4) For each segment, compute current tangent and outward normal
%   5) Evaluate traction t = sigma*n at segment midpoint
%   6) Integrate force = sum( t * ds )

    T = readtable(file);

    requiredCols = {'X','Y','U_U1','U_U2','S11','S22','S12'};
    for k = 1:numel(requiredCols)
        if ~ismember(requiredCols{k}, T.Properties.VariableNames)
            error('Column "%s" not found in file: %s', requiredCols{k}, file);
        end
    end

    X   = T.X;
    Y   = T.Y;
    U1  = T.U_U1;
    U2  = T.U_U2;
    S11 = T.S11;
    S22 = T.S22;
    S12 = T.S12;

    x = X + U1;
    y = Y + U2;

    tol = 1e-8 * max(1, max(abs([X;Y])));

    switch lower(edgeName)
        case 'right'
            edgeVal = max(X);
            idx = abs(X - edgeVal) < tol;
            sortCoord = Y(idx);
            outwardHint = [1;0];

        case 'left'
            edgeVal = min(X);
            idx = abs(X - edgeVal) < tol;
            sortCoord = Y(idx);
            outwardHint = [-1;0];

        case 'top'
            edgeVal = max(Y);
            idx = abs(Y - edgeVal) < tol;
            sortCoord = X(idx);
            outwardHint = [0;1];

        case 'bottom'
            edgeVal = min(Y);
            idx = abs(Y - edgeVal) < tol;
            sortCoord = X(idx);
            outwardHint = [0;-1];

        otherwise
            error('Unknown edgeName "%s". Use right, left, top, or bottom.', edgeName);
    end

    if ~any(idx)
        error('No nodes found on edge "%s".', edgeName);
    end

    xb   = x(idx);
    yb   = y(idx);
    S11b = S11(idx);
    S22b = S22(idx);
    S12b = S12(idx);

    [~, order] = sort(sortCoord);
    xb   = xb(order);
    yb   = yb(order);
    S11b = S11b(order);
    S22b = S22b(order);
    S12b = S12b(order);

    % Remove duplicate points if present
    pts = [xb(:), yb(:)];
    [~, ia] = unique(round(pts*1e12)/1e12, 'rows', 'stable');
    xb   = xb(ia);
    yb   = yb(ia);
    S11b = S11b(ia);
    S22b = S22b(ia);
    S12b = S12b(ia);

    nSeg = numel(xb) - 1;
    if nSeg < 1
        error('Not enough boundary points found on edge "%s".', edgeName);
    end

    Rx = 0;
    Ry = 0;

    for i = 1:nSeg
        dx = xb(i+1) - xb(i);
        dy = yb(i+1) - yb(i);
        ds = hypot(dx, dy);

        if ds < 1e-14
            continue;
        end

        % tangent
        txg = dx / ds;
        tyg = dy / ds;

        % candidate normal (90 deg rotation)
        n = [tyg; -txg];

        % choose sign so it points outward
        if dot(n, outwardHint) < 0
            n = -n;
        end

        nx = n(1);
        ny = n(2);

        % midpoint stress (average of endpoints)
        s11m = 0.5 * (S11b(i) + S11b(i+1));
        s22m = 0.5 * (S22b(i) + S22b(i+1));
        s12m = 0.5 * (S12b(i) + S12b(i+1));

        % traction t = sigma * n
        tvec_x = s11m * nx + s12m * ny;
        tvec_y = s12m * nx + s22m * ny;

        % segment force
        Rx = Rx + tvec_x * ds;
        Ry = Ry + tvec_y * ds;
    end
end