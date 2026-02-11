function [gradX, gradY] = FEM_2D_nodal_gradient_11_2_26(X, Y, Field, Connectivity, Node_Labels)
% FEM_2D_NODAL_GRADIENT_11_2_26: Mixed-mesh compatible gradient calculation.
%
% INPUTS:
%  X, Y, Field, Node_Labels : [N x 1] Vectors
%  Connectivity : [E x 10] Matrix from CSV [ID, TypeID, N1...N8]
%                 TypeID mapping: 3=T3, 6=T6, 4=Q4, 8=Q8

    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);
    
    nN = length(X);
    AccX = zeros(nN, 1); AccY = zeros(nN, 1); Cnt  = zeros(nN, 1);
    
    % Define local coordinates for each supported type once
    localCoords = struct();
    localCoords.T3 = struct('x', [0, 1, 0], 'y', [0, 0, 1]);
    localCoords.T6 = struct('x', [0, 1, 0, 0.5, 0.5, 0], 'y', [0, 0, 1, 0, 0.5, 0.5]);
    localCoords.Q4 = struct('x', [-1, 1, 1, -1], 'y', [-1, -1, 1, 1]);
    localCoords.Q8 = struct('x', [-1, 1, 1, -1, 0, 1, 0, -1], 'y', [-1, -1, 1, 1, -1, 0, 1, 0]);

    for e = 1:size(Connectivity, 1)
        typeID = Connectivity(e, 2); % Read TypeID from 2nd column
        nodes  = Connectivity(e, 3:end);
        
        % Determine mesh parameters based on TypeID
        switch typeID
            case 3, mType = 'T3'; numNodes = 3;
            case 6, mType = 'T6'; numNodes = 6;
            case 4, mType = 'Q4'; numNodes = 4;
            case 8, mType = 'Q8'; numNodes = 8;
            otherwise, continue; % Unsupported or empty row
        end
        
        % Map nodes to indices
        valid_nodes = nodes(1:numNodes);
        idx = Map(valid_nodes);
        if any(idx == 0), continue; end
        
        xe = X(idx); ye = Y(idx); fe = Field(idx);
        xN = localCoords.(mType).x; yN = localCoords.(mType).y;
        
        for n = 1:numNodes
            [dN_dxi, dN_deta] = get_shape_derivatives(xN(n), yN(n), mType);
            
            % Jacobian
            J = [dN_dxi * xe, dN_dxi * ye; 
                 dN_deta * xe, dN_deta * ye];
            
            % Global Gradient of Shape Functions
            g_N = J \ [dN_dxi; dN_deta];
            
            p = idx(n);
            AccX(p) = AccX(p) + g_N(1, :) * fe;
            AccY(p) = AccY(p) + g_N(2, :) * fe;
            Cnt(p)  = Cnt(p) + 1;
        end
    end
    
    gradX = AccX ./ max(Cnt, 1);
    gradY = AccY ./ max(Cnt, 1);
end

% ... [get_shape_derivatives and sd_Q8 helper functions remain same] ...
