function [Integral_Val, Total_Area, dArea_Vector] = FEM_2D_GP_Integral_12_2_26(X, Y, Field, Connectivity, Node_Labels)
% FEM_2D_GP_INTEGRAL_12_2_26: Calculates surface integral of a field.
%
% AUTOMATION UPDATE:
%   - Automatically detects if 'Field' is Nodal or Gauss Point data.
%   - If Nodal (Size == NumNodes): Interpolates to GP using Shape Functions.
%   - If GP (Size != NumNodes): Uses GP values directly (Sequential).
%
% INPUTS:
%  X, Y         : [N x 1] Nodal Coordinates
%  Field        : [N x 1] OR [Total_IPs x 1] Data Vector
%  Connectivity : [E x M] Element Matrix
%  Node_Labels  : [N x 1] User Node IDs
%
% OUTPUTS:
%  Integral_Val : Scalar result of Int(Field) dA
%  Total_Area   : Scalar result of Int(1) dA
%  dArea_Vector : [Total_IPs x 1] Vector of dA for each Gauss Point

    % 1. Map Node Labels to Indices
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2. Input Type Detection
    nNodes = length(X);
    nField = length(Field);
    isNodal = (nField == nNodes);
    
    if isNodal
        fprintf('   [Integral] Mode: Nodal Interpolation (%d Nodes)\n', nNodes);
    else
        fprintf('   [Integral] Mode: Direct Gauss Point Summation (%d Points)\n', nField);
    end

    % 3. Define Gauss Rules (Row-Major / Tensor Product)
    GaussRules = struct();
    
    % T3
    GaussRules.T3.xi = [1/3]; GaussRules.T3.eta = [1/3]; GaussRules.T3.w = [0.5];
    % T6
    GaussRules.T6.xi = [1/6, 2/3, 1/6]; GaussRules.T6.eta = [1/6, 1/6, 2/3]; GaussRules.T6.w = [1/6, 1/6, 1/6];
    % Q4 (Row-Major)
    pt = 1/sqrt(3);
    GaussRules.Q4.xi = [-pt, pt, -pt, pt]; GaussRules.Q4.eta = [-pt, -pt, pt, pt]; GaussRules.Q4.w = ones(1,4);
    % Q8 (Row-Major)
    sq35 = sqrt(3/5); r = sq35;
    r1_xi = [-r, 0, r]; r1_eta = [-r, -r, -r];
    r2_xi = [-r, 0, r]; r2_eta = [0, 0, 0];
    r3_xi = [-r, 0, r]; r3_eta = [r, r, r];
    GaussRules.Q8.xi = [r1_xi, r2_xi, r3_xi]; GaussRules.Q8.eta = [r1_eta, r2_eta, r3_eta];
    wA = 5/9; wB = 8/9;
    w_r1 = [wA*wA, wB*wA, wA*wA]; w_r2 = [wA*wB, wB*wB, wA*wB]; w_r3 = [wA*wA, wB*wA, wA*wA];
    GaussRules.Q8.w = [w_r1, w_r2, w_r3];

    % 4. Initialize Accumulators
    Integral_Val = 0.0;
    Total_Area = 0.0;
    
    % Pre-allocate output vector (estimate size)
    nElem = size(Connectivity, 1);
    estSize = nElem * 9; 
    dArea_Vector = zeros(estSize, 1); 
    rowCount = 0;
    
    % Global Counter for GP Mode
    gp_global_idx = 0;

    % 5. Loop over Elements
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
        
        % If Nodal Mode: Fetch nodal values for this element
        if isNodal
            fe = Field(idx);
        end
        
        rule = GaussRules.(mType);
        numIPs = length(rule.xi);
        
        for ip = 1:numIPs
            xi  = rule.xi(ip);
            eta = rule.eta(ip);
            w   = rule.w(ip);
            
            % Increment global counter
            gp_global_idx = gp_global_idx + 1;
            
            % A. Geometry (Jacobian)
            [N, dN_dxi, dN_deta] = get_shape_funcs_and_derivs(xi, eta, mType);
            J = [dN_dxi * xe,  dN_dxi * ye;
                 dN_deta * xe, dN_deta * ye];
            detJ = det(J);
            dArea = abs(detJ) * w;
            
            % B. Get Field Value
            if isNodal
                val_gp = N * fe; % Interpolate
            else
                % Safety check for GP mode
                if gp_global_idx > nField
                    error('Error: Mesh has more Gauss Points than the Input Field vector (%d vs %d).', gp_global_idx, nField);
                end
                val_gp = Field(gp_global_idx); % Direct Access
            end
            
            % C. Accumulate
            Integral_Val = Integral_Val + val_gp * dArea;
            Total_Area   = Total_Area + dArea;
            
            % Store dArea
            rowCount = rowCount + 1;
            if rowCount > length(dArea_Vector) % Expand if needed
                dArea_Vector = [dArea_Vector; zeros(nElem, 1)]; 
            end
            dArea_Vector(rowCount) = dArea;
        end
    end
    
    % Final Cleanup
    dArea_Vector = dArea_Vector(1:rowCount);
    
    % Warning for mismatched GP data size
    if ~isNodal && gp_global_idx ~= nField
        warning('Input Field has %d points, but mesh processed %d Gauss Points. Integral might be misaligned.', nField, gp_global_idx);
    end
end

%% Helper: Shape Functions
function [N, dN_dxi, dN_deta] = get_shape_funcs_and_derivs(xi, eta, type)
    switch upper(type)
        case 'T3' 
            N = [1-xi-eta, xi, eta]; dN_dxi = [-1, 1, 0]; dN_deta = [-1, 0, 1];
        case 'T6'
            L1 = 1-xi-eta; L2 = xi; L3 = eta;
            N = [L1*(2*L1-1), L2*(2*L2-1), L3*(2*L3-1), 4*L1*L2, 4*L2*L3, 4*L3*L1];
            dN_dxi = [1-4*L1, 4*L2-1, 0, 4*(L1-L2), 4*L3, -4*L3];
            dN_deta = [1-4*L1, 0, 4*L3-1, -4*L2, 4*L2, 4*(L1-L3)];
        case 'Q4' 
            N = 0.25 * [(1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta)];
            dN_dxi = 0.25*[-(1-eta),(1-eta),(1+eta),-(1+eta)]; dN_deta = 0.25*[-(1-xi),-(1+xi),(1+xi),(1-xi)];
        case 'Q8'
            [dN_dxi, dN_deta] = sd_Q8(xi, eta);
            N1=0.25*(1-xi)*(1-eta)*(-xi-eta-1); N2=0.25*(1+xi)*(1-eta)*(xi-eta-1);
            N3=0.25*(1+xi)*(1+eta)*(xi+eta-1); N4=0.25*(1-xi)*(1+eta)*(-xi+eta-1);
            N5=0.5*(1-xi^2)*(1-eta); N6=0.5*(1+xi)*(1-eta^2);
            N7=0.5*(1-xi^2)*(1+eta); N8=0.5*(1-xi)*(1-eta^2);
            N = [N1, N2, N3, N4, N5, N6, N7, N8];
    end
end
function [dN_dxi, dN_deta] = sd_Q8(xi, eta)
    dN_dxi = 0.25 * [(1-eta)*(2*xi+eta), (1-eta)*(2*xi-eta), (1+eta)*(2*xi+eta), (1+eta)*(2*xi-eta), -4*xi*(1-eta), 2*(1-eta^2), -4*xi*(1+eta), -2*(1-eta^2)];
    dN_deta = 0.25 * [(1-xi)*(2*eta+xi), (1+xi)*(2*eta-xi), (1+xi)*(2*eta+xi), (1-xi)*(2*eta-xi), -2*(1-xi^2), -4*eta*(1+xi), 2*(1-xi^2), -4*eta*(1-xi)];
end
