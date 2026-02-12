function GP_Data = FEM_2D_GP_gradient_12_2_26(X, Y, Field, Connectivity, Node_Labels)
% FEM_2D_GP_GRADIENT_12_2_26: Computes gradients at Gauss Integration Points.
%
% This function mimics Abaqus behavior by calculating the field value and 
% its spatial derivatives exactly at the integration points of the elements.
%
% INPUTS:
%  X, Y         : [N x 1] Nodal Coordinates
%  Field        : [N x 1] Nodal Field Value (e.g., U1 or U2)
%  Connectivity : [E x M] Element Matrix [ElemID, TypeID, N1, N2...]
%  Node_Labels  : [N x 1] User Node IDs
%
% OUTPUTS:
%  GP_Data      : Matrix [Total_IPs x 7]
%                 Col 1: Element Label
%                 Col 2: Integration Point Label (1, 2, 3...)
%                 Col 3: X coordinate at Gauss Point
%                 Col 4: Y coordinate at Gauss Point
%                 Col 5: Field Value at Gauss Point (Interpolated U)
%                 Col 6: dField/dX (Gradient X)
%                 Col 7: dField/dY (Gradient Y)

    % 1. Map Node Labels to Indices
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2. Define Gauss Integration Rules (Locations & Weights)
    % Coordinates (xi, eta) match Abaqus Full Integration defaults.
    GaussRules = struct();
    
    % T3 (Triangle, Linear): 1 Point at centroid
    GaussRules.T3.xi  = [1/3];
    GaussRules.T3.eta = [1/3];
    
    % T6 (Triangle, Quadratic): 3 Points
    % Standard Hammer rule
    GaussRules.T6.xi  = [1/6, 2/3, 1/6];
    GaussRules.T6.eta = [1/6, 1/6, 2/3];
    
    % Q4 (Quad, Linear): 4 Points (2x2)
    pt = 1/sqrt(3);
    GaussRules.Q4.xi  = [-pt,  pt,  pt, -pt];
    GaussRules.Q4.eta = [-pt, -pt,  pt,  pt];
    
    % Q8 (Quad, Quadratic): 9 Points (3x3 Full Integration)
    % Note: Some Abaqus settings use Reduced (4 pts), but Standard is 9.
    sq35 = sqrt(3/5);
    r = [-sq35, 0, sq35];
    [G_xi, G_eta] = meshgrid(r, r);
    GaussRules.Q8.xi  = G_xi(:)'; % Vectorize
    GaussRules.Q8.eta = G_eta(:)';
    
    % 3. Pre-allocate Output
    % Estimate size (assume max 9 IPs per element to be safe)
    nElem = size(Connectivity, 1);
    estRows = nElem * 9; 
    GP_Data = zeros(estRows, 7);
    rowCount = 0;

    % 4. Loop over Elements
    [~, totalCols] = size(Connectivity);
    
    for e = 1:nElem
        % Basic checks
        if totalCols < 2, continue; end
        
        ElemID = Connectivity(e, 1);
        TypeID = Connectivity(e, 2);
        
        % Determine Element Type
        switch TypeID
            case 3, mType = 'T3'; numNodes = 3;
            case 6, mType = 'T6'; numNodes = 6;
            case 4, mType = 'Q4'; numNodes = 4;
            case 8, mType = 'Q8'; numNodes = 8;
            otherwise, continue;
        end
        
        if totalCols < (2 + numNodes), continue; end
        
        % Get Nodes
        nodes = Connectivity(e, 3:2+numNodes);
        idx = Map(round(nodes));
        
        if any(idx == 0) || any(isnan(nodes)), continue; end
        
        % Element Nodal Coordinates & Field Values
        xe = X(idx);
        ye = Y(idx);
        fe = Field(idx);
        
        % Loop over Gauss Points for this element
        rule = GaussRules.(mType);
        numIPs = length(rule.xi);
        
        for ip = 1:numIPs
            xi  = rule.xi(ip);
            eta = rule.eta(ip);
            
            % 1. Get Shape Functions (N) and Derivatives (dN_dxi)
            [N, dN_dxi, dN_deta] = get_shape_funcs_and_derivs(xi, eta, mType);
            
            % 2. Interpolate Coordinate & Field at Gauss Point
            % X_gp = sum(N_i * x_i)
            X_gp = N * xe;
            Y_gp = N * ye;
            val_gp = N * fe;
            
            % 3. Jacobian Calculation
            % J = [dx/dxi  dy/dxi ]
            %     [dx/deta dy/deta]
            J = [dN_dxi * xe,  dN_dxi * ye;
                 dN_deta * xe, dN_deta * ye];
             
            % Check for singularity (Zero Area)
            detJ = det(J);
            if abs(detJ) < 1e-12
                % Handle bad element (store NaN)
                grad_glob = [NaN; NaN];
            else
                % 4. Global Derivatives
                % [dN/dx; dN/dy] = inv(J) * [dN/dxi; dN/deta]
                % Gradient of U = [dN/dx * u; dN/dy * u]
                % Efficient solve:
                g_N = J \ [dN_dxi; dN_deta]; % 2 x numNodes
                grad_glob = g_N * fe;        % 2 x 1 vector
            end
            
            % 5. Store Data
            rowCount = rowCount + 1;
            GP_Data(rowCount, :) = [ElemID, ip, X_gp, Y_gp, val_gp, grad_glob(1), grad_glob(2)];
        end
    end
    
    % Trim unused rows
    GP_Data = GP_Data(1:rowCount, :);
end

%% Helper: Shape Functions AND Derivatives
function [N, dN_dxi, dN_deta] = get_shape_funcs_and_derivs(xi, eta, type)
    % N: Row vector [1 x numNodes]
    % dN: Row vector [1 x numNodes]
    
    switch upper(type)
        case 'T3' % Linear Triangle
            % N1=xi, N2=eta, N3=1-xi-eta (Standard Area coords)
            % Be careful with node ordering. 
            % Abaqus T3: 1(0,0), 2(1,0), 3(0,1)? Usually 1-2-3 ccw.
            % Let's use standard isoparametric T3: N1=1-xi-eta, N2=xi, N3=eta
            N = [1-xi-eta, xi, eta];
            dN_dxi  = [-1, 1, 0];
            dN_deta = [-1, 0, 1];
            
        case 'T6' % Quadratic Triangle
            L1 = 1 - xi - eta; L2 = xi; L3 = eta;
            % Shape functions (Standard convention)
            N = [L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1];
            
            % Derivatives (Chain rule)
            % dL1/dxi = -1, dL2/dxi = 1, dL3/dxi = 0
            dN_dxi = [1-4*L1, 4*L2-1, 0, 4*(L1-L2), 4*L3, -4*L3];
            % dL1/deta = -1, dL2/deta = 0, dL3/deta = 1
            dN_deta = [1-4*L1, 0, 4*L3-1, -4*L2, 4*L2, 4*(L1-L3)];
            
        case 'Q4' % Linear Quad
            % N_i = 0.25*(1 +/- xi)*(1 +/- eta)
            % Order: (-1,-1), (1,-1), (1,1), (-1,1)
            N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
            
            dN_dxi  = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi),  -(1+xi),  (1+xi),   (1-xi)];
            
        case 'Q8' % Quadratic Quad (Serendipity)
            % Corners: 1-4, Mids: 5-8
            % 1:(-1,-1), 2:(1,-1), 3:(1,1), 4:(-1,1)
            % 5:(0,-1), 6:(1,0), 7:(0,1), 8:(-1,0)
            
            % Corner nodes
            N1 = 0.25*(1-xi)*(1-eta)*(-xi-eta-1);
            N2 = 0.25*(1+xi)*(1-eta)*(xi-eta-1);
            N3 = 0.25*(1+xi)*(1+eta)*(xi+eta-1);
            N4 = 0.25*(1-xi)*(1+eta)*(-xi+eta-1);
            % Midside nodes
            N5 = 0.5*(1-xi^2)*(1-eta);
            N6 = 0.5*(1+xi)*(1-eta^2);
            N7 = 0.5*(1-xi^2)*(1+eta);
            N8 = 0.5*(1-xi)*(1-eta^2);
            
            N = [N1, N2, N3, N4, N5, N6, N7, N8];
            
            % Reuse derivative logic from previous file or explicit
            [dN_dxi, dN_deta] = sd_Q8(xi, eta);
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
