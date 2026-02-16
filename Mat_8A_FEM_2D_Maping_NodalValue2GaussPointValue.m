function GP_Data = Mat_8A_FEM_2D_Maping_NodalValue2GaussPointValue(X, Y, Field, Connectivity, Node_Labels)
% FEM_2D_GP_GRADIENT_12_2_26: Computes Gauss  Points values from nodal points.
%
% FIXED: Q4/Q8 Ordering Updated based on user feedback.
%        Q4: Switched from CCW to Tensor Product (BL, BR, TL, TR).
%        Q8: Reverted to Row-Major (Tensor Product) to match Q4 logic.
%
% INPUTS:
%  X, Y         : [N x 1] Nodal Coordinates
%  Field        : [N x 1] Nodal Field Value (e.g., U1 or U2)
%  Connectivity : [E x M] Element Matrix [ElemID, TypeID, N1, N2...]
%  Node_Labels  : [N x 1] User Node IDs
%
% OUTPUTS:
%  GP_Data      : Matrix [Total_IPs x 5]
%                 [ElemID, IP_ID, X_gp, Y_gp, Val_gp]

    % 1. Map Node Labels to Indices
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2. Define Gauss Integration Rules (Locations & Weights)
    GaussRules = struct();
    
    % T3 (Triangle, Linear): 1 Point at centroid
    GaussRules.T3.xi  = [1/3];
    GaussRules.T3.eta = [1/3];
    
    % T6 (Triangle, Quadratic): 3 Points (Hammer)
    GaussRules.T6.xi  = [1/6, 2/3, 1/6];
    GaussRules.T6.eta = [1/6, 1/6, 2/3];
    
    % Q4 (Quad, Linear): 4 Points (2x2)
    % PREVIOUS (CCW): BL, BR, TR, TL. (Failed)
    % NEW (Tensor Product): BL, BR, TL, TR.
    % Order: (-,-), (+,-), (-,+), (+,+)
    pt = 1/sqrt(3);
    GaussRules.Q4.xi  = [-pt,  pt, -pt,  pt];
    GaussRules.Q4.eta = [-pt, -pt,  pt,  pt];
    
    % Q8 (Quad, Quadratic): 9 Points (3x3 Full Integration)
    % PREVIOUS Attempt (Row-Major) failed? Re-verifying logic.
    % Since Q4 works with Tensor Product (Row-Major), Q8 MUST be Row-Major.
    % Order: Bottom Row (L,M,R), Mid Row (L,M,R), Top Row (L,M,R).
    sq35 = sqrt(3/5);
    r = sq35; 
    
    % Row 1 (Bottom, eta = -r)
    r1_xi = [-r, 0, r]; r1_eta = [-r, -r, -r];
    % Row 2 (Mid, eta = 0)
    r2_xi = [-r, 0, r]; r2_eta = [0, 0, 0];
    % Row 3 (Top, eta = r)
    r3_xi = [-r, 0, r]; r3_eta = [r, r, r];
    
    GaussRules.Q8.xi  = [r1_xi, r2_xi, r3_xi];
    GaussRules.Q8.eta = [r1_eta, r2_eta, r3_eta];
    
    % 3. Pre-allocate Output
    nElem = size(Connectivity, 1);
    estRows = nElem * 9; 
    GP_Data = zeros(estRows, 5);
    rowCount = 0;

    % 4. Loop over Elements
    [~, totalCols] = size(Connectivity);
    
    for e = 1:nElem
        if totalCols < 2, continue; end
        
        ElemID = Connectivity(e, 1);
        TypeID = Connectivity(e, 2);
        
        switch TypeID
            case 3, mType = 'T3'; numNodes = 3;
            case 6, mType = 'T6'; numNodes = 6;
            case 4, mType = 'Q4'; numNodes = 4;
            case 8, mType = 'Q8'; numNodes = 8;
            otherwise, continue;
        end
        
        if totalCols < (2 + numNodes), continue; end
        
        nodes = Connectivity(e, 3:2+numNodes);
        idx = Map(round(nodes));
        
        if any(idx == 0) || any(isnan(nodes)), continue; end
        
        xe = X(idx);
        ye = Y(idx);
        fe = Field(idx);
        
        rule = GaussRules.(mType);
        numIPs = length(rule.xi);
        
        for ip = 1:numIPs
            xi  = rule.xi(ip);
            eta = rule.eta(ip);
            
            % 1. Get Shape Functions
            [N] = get_shape_funcs_and_derivs(xi, eta, mType);
            
            % 2. Interpolate at Gauss Point
            X_gp = N * xe;
            Y_gp = N * ye;
            val_gp = N * fe;
            
            % 5. Store
            rowCount = rowCount + 1;
            GP_Data(rowCount, :) = [ElemID, ip, X_gp, Y_gp, val_gp];
        end
    end
    GP_Data = GP_Data(1:rowCount, :);
end

%% Helper: Shape Functions
function [N] = get_shape_funcs_and_derivs(xi, eta, type)
    switch upper(type)
        case 'T3' 
            N = [1-xi-eta, xi, eta];
            
            
        case 'T6'
            L1 = 1 - xi - eta; L2 = xi; L3 = eta;
            N = [L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1];
           
            
        case 'Q4' 
            % Nodes: (-1,-1), (1,-1), (1,1), (-1,1)
            N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
           
            
        case 'Q8'
            % Nodes 1-4 (Corners), 5-8 (Mids)
            % 1:(-1,-1), 2:(1,-1), 3:(1,1), 4:(-1,1)
            % 5:(0,-1), 6:(1,0), 7:(0,1), 8:(-1,0)
            
           
            
            % Shape Functions (reconstruct for N*xe)
            N1 = 0.25*(1-xi)*(1-eta)*(-xi-eta-1);
            N2 = 0.25*(1+xi)*(1-eta)*(xi-eta-1);
            N3 = 0.25*(1+xi)*(1+eta)*(xi+eta-1);
            N4 = 0.25*(1-xi)*(1+eta)*(-xi+eta-1);
            N5 = 0.5*(1-xi^2)*(1-eta);
            N6 = 0.5*(1+xi)*(1-eta^2);
            N7 = 0.5*(1-xi^2)*(1+eta);
            N8 = 0.5*(1-xi)*(1-eta^2);
            N = [N1, N2, N3, N4, N5, N6, N7, N8];
    end
end

