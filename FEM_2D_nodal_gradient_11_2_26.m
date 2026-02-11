function [gradX, gradY] = FEM_2D_nodal_gradient_11_2_26(X, Y, Field, Connectivity, Node_Labels, meshType)
% GET_NODAL_GRADIENT: Calculates nodal-wise gradients using FEM shape functions.
% Method: Nodal averaging of element-level gradients (Direct FEM method).
%
% ========================== INPUT TYPES ==========================
%  X            : [N x 1] Column Vector. Nodal X-coordinates (from data.X)
%  Y            : [N x 1] Column Vector. Nodal Y-coordinates (from data.Y)
%  Field        : [N x 1] Column Vector. The parameter to differentiate 
%                 (e.g., data.U_U1, data.E22, or temperature).
%  Connectivity : [E x M] Matrix. Element connectivity matrix. 
%                 Each row is an element; columns are node labels.
%  Node_Labels  : [N x 1] Column Vector. Labels/IDs for each node 
%                 (from data.Node_Label). Used for mapping IDs to indices.
%  meshType     : [String] 'T3' (3-node Tri), 'T6' (6-node Tri), 
%                 'Q4' (4-node Quad), 'Q8' (8-node Quad).
% =================================================================

    % 1. Create a label-to-index map (Handles Abaqus non-sequential IDs)
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);
    
    nN = length(X);
    AccX = zeros(nN, 1); % Accumulator for dF/dx
    AccY = zeros(nN, 1); % Accumulator for dF/dy
    Cnt  = zeros(nN, 1); % Count of elements sharing each node
    
    % 2. Setup element-specific local coordinates
    switch upper(meshType)
        case 'T3'
            xN = [0, 1, 0]; yN = [0, 0, 1]; numNodes = 3;
        case 'T6'
            xN = [0, 1, 0, 0.5, 0.5, 0]; yN = [0, 0, 1, 0, 0.5, 0.5]; numNodes = 6;
        case 'Q4'
            xN = [-1, 1, 1, -1]; yN = [-1, -1, 1, 1]; numNodes = 4;
        case 'Q8'
            xN = [-1, 1, 1, -1, 0, 1, 0, -1]; yN = [-1, -1, 1, 1, -1, 0, 1, 0]; numNodes = 8;
        otherwise
            error('Mesh type %s not supported.', meshType);
    end
    
    % 3. Loop through elements
    for e = 1:size(Connectivity, 1)
        nodes = Connectivity(e, :);
        valid_nodes = nodes(~isnan(nodes)); % Ignore NaN padding
        idx = Map(valid_nodes);
        
        % Skip if element doesn't match the expected node count
        if any(idx == 0) || length(idx) ~= numNodes
            continue; 
        end
        
        % Local element nodal data
        xe = X(idx); ye = Y(idx); fe = Field(idx);
        
        % Evaluate gradient at each node of the element
        for n = 1:numNodes
            [dN_dxi, dN_deta] = get_shape_derivatives(xN(n), yN(n), meshType);
            
            % Jacobian Matrix: J = [dx/dxi, dy/dxi; dx/deta, dy/deta]
            J = [dN_dxi * xe, dN_dxi * ye; 
                 dN_deta * xe, dN_deta * ye];
            
            % Global derivatives of shape functions: [dN/dx; dN/dy]
            g_N = J \ [dN_dxi; dN_deta];
            
            % Calculate nodal gradient contributed by this element
            p = idx(n);
            AccX(p) = AccX(p) + g_N(1, :) * fe;
            AccY(p) = AccY(p) + g_N(2, :) * fe;
            Cnt(p)  = Cnt(p) + 1;
        end
    end
    
    % 4. Final nodal averaging
    gradX = AccX ./ max(Cnt, 1);
    gradY = AccY ./ max(Cnt, 1);
end

%% --- Helper: Shape Function Derivatives ---
function [dN_dxi, dN_deta] = get_shape_derivatives(xi, eta, type)
    switch upper(type)
        case 'T3'
            dN_dxi  = [1, 0, -1];
            dN_deta = [0, 1, -1];
        case 'T6'
            L1 = xi; L2 = eta; L3 = 1 - xi - eta;
            dN_dxi  = [4*L1-1, 0, 1-4*L3, 4*L2, -4*L2, 4*(L3-L1)];
            dN_deta = [0, 4*L2-1, 1-4*L3, -4*L1, 4*L1, 4*(L3-L2)];
        case 'Q4'
            dN_dxi  = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi),  -(1+xi),  (1+xi),   (1-xi)];
        case 'Q8'
            % Derivatives for 8-node serendipity quad
            [dN_dxi, dN_deta] = sd_Q8(xi, eta);
    end
end

function [dN_dxi, dN_deta] = sd_Q8(xi, eta)
    % Internal implementation of Quadratic Quad derivatives
    % Corners: 1, 2, 3, 4 | Midpoints: 5, 6, 7, 8
    dN_dxi = 0.25 * [
        (1-eta)*(2*xi+eta), (1-eta)*(2*xi-eta), (1+eta)*(2*xi+eta), (1+eta)*(2*xi-eta), ...
        -4*xi*(1-eta), 2*(1-eta^2), -4*xi*(1+eta), -2*(1-eta^2) ];
    dN_deta = 0.25 * [
        (1-xi)*(2*eta+xi), (1+xi)*(2*eta-xi), (1+xi)*(2*eta+xi), (1-xi)*(2*eta-xi), ...
        -2*(1-xi^2), -4*eta*(1+xi), 2*(1-xi^2), -4*eta*(1-xi) ];
end