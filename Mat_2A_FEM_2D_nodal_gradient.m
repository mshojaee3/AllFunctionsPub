function [gradX, gradY] = Mat_2A_FEM_2D_nodal_gradient(X, Y, Field, Connectivity, Node_Labels)
% Mat_2A_FEM_2D_nodal_gradient: Mixed-mesh compatible gradient calculation.
%
% INPUTS:
%  X            : [N x 1] Vector of Nodal X-coordinates
%  Y            : [N x 1] Vector of Nodal Y-coordinates
%  Field        : [N x 1] Vector of Nodal field values (e.g., Displacement)
%  Connectivity : [E x M] Matrix from CSV [ElementID, TypeID, N1, N2, ...]
%                 TypeID mapping: 3=T3, 6=T6, 4=Q4, 8=Q8
%  Node_Labels  : [N x 1] Vector of Nodal IDs
%
% OUTPUTS:
%  gradX, gradY : [N x 1] Vectors of nodal gradients

    % Create a label-to-index map to handle non-sequential Abaqus IDs
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);
    
    nN = length(X);
    AccX = zeros(nN, 1); 
    AccY = zeros(nN, 1); 
    Cnt  = zeros(nN, 1);
    
    % Define local coordinates for each supported type once for performance
    localCoords = struct();
    localCoords.T3 = struct('x', [0, 1, 0], 'y', [0, 0, 1]);
    localCoords.T6 = struct('x', [0, 1, 0, 0.5, 0.5, 0], 'y', [0, 0, 1, 0, 0.5, 0.5]);
    localCoords.Q4 = struct('x', [-1, 1, 1, -1], 'y', [-1, -1, 1, 1]);
    localCoords.Q8 = struct('x', [-1, 1, 1, -1, 0, 1, 0, -1], 'y', [-1, -1, 1, 1, -1, 0, 1, 0]);

    % Determine matrix dimensions
    [~, totalCols] = size(Connectivity);

    for e = 1:size(Connectivity, 1)
        % Ensure the row has at least an ID and a TypeID column
        if totalCols < 2
            continue;
        end
        
        typeID = Connectivity(e, 2); % Read TypeID from 2nd column
        
        % Determine mesh parameters based on TypeID
        switch typeID
            case 3
                mType = 'T3'; numNodes = 3;
            case 6
                mType = 'T6'; numNodes = 6;
            case 4
                mType = 'Q4'; numNodes = 4;
            case 8
                mType = 'Q8'; numNodes = 8;
            otherwise
                continue; % Skip unsupported or empty rows
        end
        
        % Safety Check: Ensure the matrix has enough columns for this specific element
        % We start nodes at column 3, so we need at least (2 + numNodes) columns.
        if totalCols < (2 + numNodes)
            continue; 
        end
        
        % Extract only the relevant nodes for this element
        nodes = Connectivity(e, 3:2+numNodes);
        
        % Map node labels to internal array indices
        % Round is used in case readmatrix imported labels as floating point
        idx = Map(round(nodes));
        
        % Verify all nodes are valid (mapped correctly) and no NaNs exist in the connectivity slice
        if any(idx == 0) || any(isnan(nodes))
            continue; 
        end
        
        xe = X(idx); 
        ye = Y(idx); 
        fe = Field(idx);
        xN = localCoords.(mType).x; 
        yN = localCoords.(mType).y;
        
        % Evaluate gradient at each node of the element using shape functions
        for n = 1:numNodes
            [dN_dxi, dN_deta] = get_shape_derivatives(xN(n), yN(n), mType);
            
            % Compute Jacobian Matrix
            J = [dN_dxi * xe, dN_dxi * ye; 
                 dN_deta * xe, dN_deta * ye];
            
            % Compute Global Derivatives of Shape Functions
            % g_N = [dN/dx; dN/dy] = inv(J) * [dN/dxi; dN/deta]
            g_N = J \ [dN_dxi; dN_deta];
            
            % Accumulate gradient contribution at the node
            p = idx(n);
            AccX(p) = AccX(p) + g_N(1, :) * fe;
            AccY(p) = AccY(p) + g_N(2, :) * fe;
            Cnt(p)  = Cnt(p) + 1;
        end
    end
    
    % Average accumulated gradients by the number of elements sharing each node
    gradX = AccX ./ max(Cnt, 1);
    gradY = AccY ./ max(Cnt, 1);
end

%% Helper: Shape Function Derivatives
function [dN_dxi, dN_deta] = get_shape_derivatives(xi, eta, type)
    switch upper(type)
        case 'T3' % Linear Triangle
            dN_dxi  = [1, 0, -1];
            dN_deta = [0, 1, -1];
        case 'T6' % Quadratic Triangle
            L1 = xi; L2 = eta; L3 = 1 - xi - eta;
            dN_dxi  = [4*L1-1, 0, 1-4*L3, 4*L2, -4*L2, 4*(L3-L1)];
            dN_deta = [0, 4*L2-1, 1-4*L3, -4*L1, 4*L1, 4*(L3-L2)];
        case 'Q4' % Linear Quad
            dN_dxi  = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi),  -(1+xi),  (1+xi),   (1-xi)];
        case 'Q8' % Quadratic Quad
            [dN_dxi, dN_deta] = sd_Q8(xi, eta);
    end
end

%% Helper: Serendipity Quadratic Quad Derivatives
function [dN_dxi, dN_deta] = sd_Q8(xi, eta)
    % Derivatives for 8-node serendipity quad
    % Corners: 1, 2, 3, 4 | Midpoints: 5, 6, 7, 8
    dN_dxi = 0.25 * [
        (1-eta)*(2*xi+eta), (1-eta)*(2*xi-eta), (1+eta)*(2*xi+eta), (1+eta)*(2*xi-eta), ...
        -4*xi*(1-eta), 2*(1-eta^2), -4*xi*(1+eta), -2*(1-eta^2) ];
    dN_deta = 0.25 * [
        (1-xi)*(2*eta+xi), (1+xi)*(2*eta-xi), (1+xi)*(2*eta+xi), (1-xi)*(2*eta-xi), ...
        -2*(1-xi^2), -4*eta*(1+xi), 2*(1-xi^2), -4*eta*(1-xi) ];
end

