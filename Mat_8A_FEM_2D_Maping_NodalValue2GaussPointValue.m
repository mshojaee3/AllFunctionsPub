function GP_Value = Mat_8A_FEM_2D_Maping_NodalValue2GaussPointValue(X, Y, Field, Connectivity, Node_Labels)
% Mat_8A_FEM_2D_Maping_NodalValue2GaussPointValue
% Maps a nodal scalar value to Gauss-point values for 2D Abaqus elements.
%
% Example call:
%   GP_Value = Mat_8A_FEM_2D_Maping_NodalValue2GaussPointValue(X, Y, U1, Conn, Node_Labels);
%
% INPUTS:
%  X, Y         : [N x 1] Nodal coordinates
%  Field        : [N x 1] Nodal scalar field (e.g., U1, U2, Temperature, SDV at nodes, etc.)
%  Connectivity : [E x M] Element matrix [ElemID, TypeID, N1, N2...]
%  Node_Labels  : [N x 1] User node labels (Abaqus node IDs)
%
% OUTPUT:
%  GP_Value     : [Total_IPs x 5]
%                 [ElemID, IP_ID, X_gp, Y_gp, Val_gp]

    % 1) Map node labels -> indices in X/Y/Field
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2) Gauss integration point locations (same ordering as your Mat_2B)
    GaussRules = struct();

    % T3: 1 point centroid
    GaussRules.T3.xi  = 1/3;
    GaussRules.T3.eta = 1/3;

    % T6: 3-point (Hammer)
    GaussRules.T6.xi  = [1/6, 2/3, 1/6];
    GaussRules.T6.eta = [1/6, 1/6, 2/3];

    % Q4: 2x2 tensor product order: (-,-), (+,-), (-,+), (+,+)
    pt = 1/sqrt(3);
    GaussRules.Q4.xi  = [-pt,  pt, -pt,  pt];
    GaussRules.Q4.eta = [-pt, -pt,  pt,  pt];

    % Q8: 3x3 full integration, row-major (bottom row -> top row)
    r = sqrt(3/5);
    % Bottom row eta=-r
    r1_xi  = [-r, 0, r]; r1_eta = [-r, -r, -r];
    % Mid row eta=0
    r2_xi  = [-r, 0, r]; r2_eta = [ 0,  0,  0];
    % Top row eta=+r
    r3_xi  = [-r, 0, r]; r3_eta = [ r,  r,  r];

    GaussRules.Q8.xi  = [r1_xi, r2_xi, r3_xi];
    GaussRules.Q8.eta = [r1_eta, r2_eta, r3_eta];

    % 3) Preallocate output (worst-case 9 IPs per element)
    nElem = size(Connectivity, 1);
    GP_Value = zeros(nElem * 9, 5);
    rowCount = 0;

    [~, totalCols] = size(Connectivity);

    % 4) Loop elements
    for e = 1:nElem
        if totalCols < 2, continue; end

        ElemID = Connectivity(e, 1);
        TypeID = Connectivity(e, 2);

        switch TypeID
            case 3,  mType = 'T3'; numNodes = 3;
            case 6,  mType = 'T6'; numNodes = 6;
            case 4,  mType = 'Q4'; numNodes = 4;
            case 8,  mType = 'Q8'; numNodes = 8;
            otherwise
                continue; % unsupported
        end

        if totalCols < (2 + numNodes), continue; end

        nodes = Connectivity(e, 3:2+numNodes);
        if any(isnan(nodes)), continue; end

        idx = Map(round(nodes));
        if any(idx == 0), continue; end

        xe = X(idx);
        ye = Y(idx);
        fe = Field(idx);

        rule = GaussRules.(mType);
        numIPs = length(rule.xi);

        for ip = 1:numIPs
            xi  = rule.xi(ip);
            eta = rule.eta(ip);

            % Shape functions
            N = get_shape_funcs_only(xi, eta, mType);

            % Interpolate to Gauss point
            X_gp   = N * xe;
            Y_gp   = N * ye;
            Val_gp = N * fe;

            rowCount = rowCount + 1;
            GP_Value(rowCount, :) = [ElemID, ip, X_gp, Y_gp, Val_gp];
        end
    end

    GP_Value = GP_Value(1:rowCount, :);
end

%% ---- Helper: Shape functions only (matches Mat_2B definitions) ----
function N = get_shape_funcs_only(xi, eta, type)
    switch upper(type)
        case 'T3'
            N = [1 - xi - eta, xi, eta];

        case 'T6'
            L1 = 1 - xi - eta; L2 = xi; L3 = eta;
            N = [ ...
                L1*(2*L1-1), ...
                L2*(2*L2-1), ...
                L3*(2*L3-1), ...
                4*L1*L2, ...
                4*L2*L3, ...
                4*L3*L1 ];

        case 'Q4'
            % Nodes: (-1,-1), (1,-1), (1,1), (-1,1)
            N = 0.25 * [ ...
                (1-xi)*(1-eta), ...
                (1+xi)*(1-eta), ...
                (1+xi)*(1+eta), ...
                (1-xi)*(1+eta) ];

        case 'Q8'
            % 1:(-1,-1), 2:(1,-1), 3:(1,1), 4:(-1,1)
            % 5:(0,-1),  6:(1,0),  7:(0,1), 8:(-1,0)
            N1 = 0.25*(1-xi)*(1-eta)*(-xi-eta-1);
            N2 = 0.25*(1+xi)*(1-eta)*( xi-eta-1);
            N3 = 0.25*(1+xi)*(1+eta)*( xi+eta-1);
            N4 = 0.25*(1-xi)*(1+eta)*(-xi+eta-1);
            N5 = 0.5*(1-xi^2)*(1-eta);
            N6 = 0.5*(1+xi)*(1-eta^2);
            N7 = 0.5*(1-xi^2)*(1+eta);
            N8 = 0.5*(1-xi)*(1-eta^2);
            N = [N1, N2, N3, N4, N5, N6, N7, N8];

        otherwise
            error('Unsupported element type: %s', type);
    end
end
