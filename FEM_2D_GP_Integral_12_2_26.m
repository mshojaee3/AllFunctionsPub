function [Integral_Val, Total_Area] = FEM_2D_GP_Integral_12_2_26(X, Y, Field, Connectivity, Node_Labels)
% FEM_2D_GP_INTEGRAL_12_2_26: Calculates surface integral of a field.
%
% Computes Int(Field) dA over the 2D domain using Gauss Quadrature.
% Matches the Gauss Point ordering of FEM_2D_GP_gradient_12_2_26.
%
% INPUTS:
%  X, Y         : [N x 1] Nodal Coordinates
%  Field        : [N x 1] Nodal Field Value (Scalar)
%  Connectivity : [E x M] Element Matrix [ElemID, TypeID, N1, N2...]
%  Node_Labels  : [N x 1] User Node IDs
%
% OUTPUTS:
%  Integral_Val : Scalar result of Int(Field(x,y)) dA
%  Total_Area   : Scalar result of Int(1) dA (Total Domain Area)

    % 1. Map Node Labels to Indices
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2. Define Gauss Rules (Locations & Weights)
    GaussRules = struct();
    
    % --- T3 (Triangle, Linear): 1 Point ---
    % Area of Ref Triangle = 0.5
    GaussRules.T3.xi  = [1/3];
    GaussRules.T3.eta = [1/3];
    GaussRules.T3.w   = [0.5]; 
    
    % --- T6 (Triangle, Quadratic): 3 Points (Hammer) ---
    % Sum of weights = 0.5
    GaussRules.T6.xi  = [1/6, 2/3, 1/6];
    GaussRules.T6.eta = [1/6, 1/6, 2/3];
    GaussRules.T6.w   = [1/6, 1/6, 1/6];
    
    % --- Q4 (Quad, Linear): 4 Points ---
    % Area of Ref Quad = 4.0. Weights are 1.0 each.
    % Ordering: Tensor Product (Row-Major): BL, BR, TL, TR
    pt = 1/sqrt(3);
    GaussRules.Q4.xi  = [-pt,  pt, -pt,  pt];
    GaussRules.Q4.eta = [-pt, -pt,  pt,  pt];
    GaussRules.Q4.w   = [1.0, 1.0, 1.0, 1.0];
    
    % --- Q8 (Quad, Quadratic): 9 Points ---
    % Ordering: Row-Major (Tensor Product).
    sq35 = sqrt(3/5);
    r = sq35;
    
    % Locations (Same as Gradient Function)
    % Row 1 (Bottom), Row 2 (Mid), Row 3 (Top)
    r1_xi = [-r, 0, r]; r1_eta = [-r, -r, -r];
    r2_xi = [-r, 0, r]; r2_eta = [0, 0, 0];
    r3_xi = [-r, 0, r]; r3_eta = [r, r, r];
    
    GaussRules.Q8.xi  = [r1_xi, r2_xi, r3_xi];
    GaussRules.Q8.eta = [r1_eta, r2_eta, r3_eta];
    
    % Weights (Tensor Product of 1D Gauss: 5/9, 8/9, 5/9)
    wA = 5/9; wB = 8/9;
    % Row 1 (w_eta = 5/9): [A*A, B*A, A*A]
    w_r1 = [wA*wA, wB*wA, wA*wA];
    % Row 2 (w_eta = 8/9): [A*B, B*B, A*B]
    w_r2 = [wA*wB, wB*wB, wA*wB];
    % Row 3 (w_eta = 5/9): [A*A, B*A, A*A]
    w_r3 = [wA*wA, wB*wA, wA*wA];
    
    GaussRules.Q8.w = [w_r1, w_r2, w_r3];

    % 3. Initialize Accumulators
    Integral_Val = 0.0;
    Total_Area = 0.0;

    % 4. Loop over Elements
    nElem = size(Connectivity, 1);
    [~, totalCols] = size(Connectivity);
    
    for e = 1:nElem
        if totalCols < 2, continue; end
        
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
            w   = rule.w(ip);
            
            % 1. Shape Functions & Derivatives
            [N, dN_dxi, dN_deta] = get_shape_funcs_and_derivs(xi, eta, mType);
            
            % 2. Jacobian
            J = [dN_dxi * xe,  dN_dxi * ye;
                 dN_deta * xe, dN_deta * ye];
             
            detJ = det(J);
            
            % 3. Interpolate Field
            val_gp = N * fe;
            
            % 4. Accumulate Integral
            % dA = det(J) * dxi * deta (captured by w)
            dArea = abs(detJ) * w;
            
            Integral_Val = Integral_Val + val_gp * dArea;
            Total_Area   = Total_Area + dArea;
        end
    end
end

%% Helper: Shape Functions (Identical to Gradient Function)
function [N, dN_dxi, dN_deta] = get_shape_funcs_and_derivs(xi, eta, type)
    switch upper(type)
        case 'T3' 
            N = [1-xi-eta, xi, eta];
            dN_dxi  = [-1, 1, 0];
            dN_deta = [-1, 0, 1];
            
        case 'T6'
            L1 = 1 - xi - eta; L2 = xi; L3 = eta;
            N = [L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1];
            dN_dxi = [1-4*L1, 4*L2-1, 0, 4*(L1-L2), 4*L3, -4*L3];
            dN_deta = [1-4*L1, 0, 4*L3-1, -4*L2, 4*L2, 4*(L1-L3)];
            
        case 'Q4' 
            N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
            dN_dxi  = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi),  -(1+xi),  (1+xi),   (1-xi)];
            
        case 'Q8'
            % Derivatives
            [dN_dxi, dN_deta] = sd_Q8(xi, eta);
            % Shape Functions
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

function [dN_dxi, dN_deta] = sd_Q8(xi, eta)
    dN_dxi = 0.25 * [
        (1-eta)*(2*xi+eta), (1-eta)*(2*xi-eta), (1+eta)*(2*xi+eta), (1+eta)*(2*xi-eta), ...
        -4*xi*(1-eta), 2*(1-eta^2), -4*xi*(1+eta), -2*(1-eta^2) ];
    dN_deta = 0.25 * [
        (1-xi)*(2*eta+xi), (1+xi)*(2*eta-xi), (1+xi)*(2*eta+xi), (1-xi)*(2*eta-xi), ...
        -2*(1-xi^2), -4*eta*(1+xi), 2*(1-xi^2), -4*eta*(1-xi) ];
end
