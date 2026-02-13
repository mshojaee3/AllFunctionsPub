function [Integral_Val, Total_Area] = Mat_3A_FEM_2D_GaussPoint_Integral(X, Y, Field_GP, Connectivity, Node_Labels)
% FEM_2D_GP_INTEGRAL_12_2_26: Integrates a Gauss Point field over the domain.
%
% INPUTS:
%  X, Y         : [N x 1] Nodal Coordinates (for calculating Jacobian/Area)
%  Field_GP     : [N_IP x 1] Field values exactly at Integration Points
%                 (Must match Abaqus output order: Elem 1 IPs, Elem 2 IPs...)
%  Connectivity : [E x M] Element Matrix
%  Node_Labels  : [N x 1] User Node IDs
%
% OUTPUTS:
%  Integral_Val : Sum(Field_i * dA_i)
%  Total_Area   : Sum(dA_i)

    % 1. Map Node Labels to Indices
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2. Check Input Size
    nGP_Input = length(Field_GP);

    % 3. Define Gauss Rules (Row-Major / Tensor Product to match Abaqus)
    GaussRules = struct();
    
    % T3
    GaussRules.T3.xi = [1/3]; GaussRules.T3.eta = [1/3]; GaussRules.T3.w = [0.5];
    % T6
    GaussRules.T6.xi = [1/6, 2/3, 1/6]; GaussRules.T6.eta = [1/6, 1/6, 2/3]; GaussRules.T6.w = [1/6, 1/6, 1/6];
    % Q4 (Row-Major: BL, BR, TL, TR)
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
    
    global_ip_idx = 0; % Tracks position in Field_GP vector

    % 5. Loop over Elements
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
        
        rule = GaussRules.(mType);
        numIPs = length(rule.xi);
        
        for ip = 1:numIPs
            global_ip_idx = global_ip_idx + 1;
            
            % Safety Break
            if global_ip_idx > nGP_Input
                error('Mismatch: Mesh has more Integration Points than the Input Field (%d vs %d). Check export.', global_ip_idx, nGP_Input);
            end

            xi  = rule.xi(ip);
            eta = rule.eta(ip);
            w   = rule.w(ip);
            
            % 1. Calculate Area Jacobian (dA)
            [~, dN_dxi, dN_deta] = get_shape_funcs_and_derivs(xi, eta, mType);
            J = [dN_dxi * xe,  dN_dxi * ye;
                 dN_deta * xe, dN_deta * ye];
            dArea = abs(det(J)) * w;
            
            % 2. Retrieve Field Value Directly
            val = Field_GP(global_ip_idx);
            
            % 3. Integrate
            Integral_Val = Integral_Val + val * dArea;
            Total_Area   = Total_Area + dArea;
        end
    end
    
    if global_ip_idx ~= nGP_Input
        warning('Mesh processed %d IPs, but Input Field had %d. Some data may have been ignored.', global_ip_idx, nGP_Input);
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
