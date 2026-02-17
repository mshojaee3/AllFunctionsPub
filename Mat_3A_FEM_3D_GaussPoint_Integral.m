function [Integral_Val, Total_Vol] = Mat_3A_FEM_3D_GaussPoint_Integral(X, Y, Z, Field_GP, Connectivity, Node_Labels)
% FEM_3D_GP_INTEGRAL: Integrates a Gauss Point field over a 3D domain.
%
% Fully supports Abaqus Standard, Reduced, Hybrid (H), and Modified (M) Elements:
%  - C3D8 / C3D8R / C3D8H / C3D8RH (8-node Hexahedron)
%  - C3D20 / C3D20R / C3D20H / C3D20RH (20-node Hexahedron)
%  - C3D4 / C3D4H (4-node Tetrahedron)
%  - C3D10 / C3D10M / C3D10MH (10-node Tetrahedron)
%  - C3D6 (6-node Wedge/Prism)
%
% INPUTS:
%  X, Y, Z      : [N x 1] Nodal Coordinates
%  Field_GP     : [N_IP x 1] Field values exactly at Integration Points
%  Connectivity : [E x M] Element Matrix [ElemID, NumNodes, N1, N2...]
%  Node_Labels  : [N x 1] User Node IDs

    % 1. Map Node Labels to Indices
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2. Check Input Size
    nGP_Input = length(Field_GP);
    nElem = size(Connectivity, 1);
    avg_ips = round(nGP_Input / nElem);

    % 3. Define 3D Gauss Rules
    GaussRules = build_gauss_rules();

    % 4. Initialize Accumulators
    Integral_Val = 0.0;
    Total_Vol = 0.0;
    global_ip_idx = 0; % Tracks position in Field_GP vector

    % 5. Loop over Elements
    [~, totalCols] = size(Connectivity);
    
    for e = 1:nElem
        if totalCols < 2, continue; end
        
        % TypeID is the number of nodes (8, 20, 4, 10, 6)
        TypeID = Connectivity(e, 2);
        
        switch TypeID
            case 8
                mType = 'HEX8'; numNodes = 8;
                rule = GaussRules.HEX8_8;
                if avg_ips == 1, rule = GaussRules.HEX8_1; end
            case 20
                mType = 'HEX20'; numNodes = 20;
                rule = GaussRules.HEX20_27;
                if avg_ips == 8, rule = GaussRules.HEX20_8; end
            case 4
                mType = 'TET4'; numNodes = 4;
                rule = GaussRules.TET4_1;
            case 10
                mType = 'TET10'; numNodes = 10;
                rule = GaussRules.TET10_4;
            case 6
                mType = 'WEDGE6'; numNodes = 6;
                rule = GaussRules.WEDGE6_2;
            otherwise
                continue; % Skip unsupported/1D/2D elements
        end
        
        if totalCols < (2 + numNodes), continue; end
        
        nodes = Connectivity(e, 3:2+numNodes);
        idx = Map(round(nodes));
        if any(idx == 0) || any(isnan(nodes)), continue; end
        
        xe = X(idx); ye = Y(idx); ze = Z(idx);
        numIPs = length(rule.xi);
        
        for ip = 1:numIPs
            global_ip_idx = global_ip_idx + 1;
            
            if global_ip_idx > nGP_Input
                error('Mismatch: Mesh has more Integration Points than the Input Field.');
            end

            xi = rule.xi(ip); eta = rule.eta(ip); zeta = rule.zeta(ip); w = rule.w(ip);
            
            % 1. Calculate 3D Volume Jacobian (dV)
            [~, dN_dxi, dN_deta, dN_dzeta] = get_shape_funcs_3D(xi, eta, zeta, mType);
            
            J = [dN_dxi * xe,   dN_dxi * ye,   dN_dxi * ze;
                 dN_deta * xe,  dN_deta * ye,  dN_deta * ze;
                 dN_dzeta * xe, dN_dzeta * ye, dN_dzeta * ze];
                 
            dVol = abs(det(J)) * w;
            
            % 2. Integrate
            val = Field_GP(global_ip_idx);
            Integral_Val = Integral_Val + val * dVol;
            Total_Vol    = Total_Vol + dVol;
        end
    end
    
    if global_ip_idx ~= nGP_Input
        warning('Mesh processed %d IPs, but Input Field had %d.', global_ip_idx, nGP_Input);
    end
end

%% Helper: Build Gauss Quadrature Rules
function GR = build_gauss_rules()
    GR = struct();
    
    % HEX8 / C3D8R (1 IP)
    GR.HEX8_1.xi = 0; GR.HEX8_1.eta = 0; GR.HEX8_1.zeta = 0; GR.HEX8_1.w = 8;
    
    % HEX8 / C3D8 (8 IPs)
    pt = 1/sqrt(3);
    GR.HEX8_8.xi   = [-pt, pt, -pt, pt, -pt, pt, -pt, pt];
    GR.HEX8_8.eta  = [-pt, -pt, pt, pt, -pt, -pt, pt, pt];
    GR.HEX8_8.zeta = [-pt, -pt, -pt, -pt, pt, pt, pt, pt];
    GR.HEX8_8.w    = ones(1,8);

    % HEX20 / C3D20R (8 IPs - Same locations as HEX8 Full)
    GR.HEX20_8 = GR.HEX8_8;

    % HEX20 / C3D20 (27 IPs - Full Integration)
    pts3 = [-sqrt(3/5), 0, sqrt(3/5)];
    wts3 = [5/9, 8/9, 5/9];
    [X,Y,Z] = ndgrid(pts3, pts3, pts3);
    [WX,WY,WZ] = ndgrid(wts3, wts3, wts3);
    GR.HEX20_27.xi = X(:)'; GR.HEX20_27.eta = Y(:)'; GR.HEX20_27.zeta = Z(:)';
    GR.HEX20_27.w = (WX .* WY .* WZ)';

    % TET4 / C3D4 (1 IP)
    GR.TET4_1.xi = 1/4; GR.TET4_1.eta = 1/4; GR.TET4_1.zeta = 1/4; GR.TET4_1.w = 1/6;

    % TET10 / C3D10 (4 IPs)
    a = (5+3*sqrt(5))/20; b = (5-sqrt(5))/20;
    GR.TET10_4.xi   = [b, a, b, b];
    GR.TET10_4.eta  = [b, b, a, b];
    GR.TET10_4.zeta = [b, b, b, a];
    GR.TET10_4.w    = ones(1,4) * (1/24); % V_ref = 1/6

    % WEDGE6 / C3D6 (2 IPs)
    GR.WEDGE6_2.xi   = [1/3, 1/3];
    GR.WEDGE6_2.eta  = [1/3, 1/3];
    GR.WEDGE6_2.zeta = [-pt, pt];
    GR.WEDGE6_2.w    = [0.5, 0.5];
end

%% Helper: 3D Shape Functions & Derivatives
function [N, dN_dxi, dN_deta, dN_dzeta] = get_shape_funcs_3D(xi, eta, zeta, type)
    switch upper(type)
        case 'HEX8'
            xi_i   = [-1, 1, 1,-1,-1, 1, 1,-1];
            eta_i  = [-1,-1, 1, 1,-1,-1, 1, 1];
            zeta_i = [-1,-1,-1,-1, 1, 1, 1, 1];
            N = zeros(1,8); dN_dxi = zeros(1,8); dN_deta = zeros(1,8); dN_dzeta = zeros(1,8);
            for i=1:8
                N(i)        = 1/8 * (1 + xi*xi_i(i)) * (1 + eta*eta_i(i)) * (1 + zeta*zeta_i(i));
                dN_dxi(i)   = 1/8 * xi_i(i) * (1 + eta*eta_i(i)) * (1 + zeta*zeta_i(i));
                dN_deta(i)  = 1/8 * eta_i(i) * (1 + xi*xi_i(i)) * (1 + zeta*zeta_i(i));
                dN_dzeta(i) = 1/8 * zeta_i(i) * (1 + xi*xi_i(i)) * (1 + eta*eta_i(i));
            end
            
        case 'HEX20' % Serendipity Hexahedron
            xi_i   = [-1, 1, 1,-1,-1, 1, 1,-1,  0, 1, 0,-1,  0, 1, 0,-1, -1, 1, 1,-1];
            eta_i  = [-1,-1, 1, 1,-1,-1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1,-1, 1, 1];
            zeta_i = [-1,-1,-1,-1, 1, 1, 1, 1, -1,-1,-1,-1,  1, 1, 1, 1,  0, 0, 0, 0];
            N = zeros(1,20); dN_dxi = zeros(1,20); dN_deta = zeros(1,20); dN_dzeta = zeros(1,20);
            
            for i=1:8 % Corners
                N(i) = 1/8 * (1+xi*xi_i(i))*(1+eta*eta_i(i))*(1+zeta*zeta_i(i)) * (xi*xi_i(i) + eta*eta_i(i) + zeta*zeta_i(i) - 2);
                dN_dxi(i) = 1/8 * xi_i(i)*(1+eta*eta_i(i))*(1+zeta*zeta_i(i)) * (2*xi*xi_i(i) + eta*eta_i(i) + zeta*zeta_i(i) - 1);
                dN_deta(i)= 1/8 * eta_i(i)*(1+xi*xi_i(i))*(1+zeta*zeta_i(i)) * (xi*xi_i(i) + 2*eta*eta_i(i) + zeta*zeta_i(i) - 1);
                dN_dzeta(i)=1/8 * zeta_i(i)*(1+xi*xi_i(i))*(1+eta*eta_i(i)) * (xi*xi_i(i) + eta*eta_i(i) + 2*zeta*zeta_i(i) - 1);
            end
            for i=9:20 % Midsides
                if xi_i(i) == 0
                    N(i) = 1/4 * (1-xi^2)*(1+eta*eta_i(i))*(1+zeta*zeta_i(i));
                    dN_dxi(i) = -1/2 * xi * (1+eta*eta_i(i))*(1+zeta*zeta_i(i));
                    dN_deta(i)= 1/4 * eta_i(i)*(1-xi^2)*(1+zeta*zeta_i(i));
                    dN_dzeta(i)=1/4 * zeta_i(i)*(1-xi^2)*(1+eta*eta_i(i));
                elseif eta_i(i) == 0
                    N(i) = 1/4 * (1+xi*xi_i(i))*(1-eta^2)*(1+zeta*zeta_i(i));
                    dN_dxi(i) = 1/4 * xi_i(i)*(1-eta^2)*(1+zeta*zeta_i(i));
                    dN_deta(i)= -1/2 * eta * (1+xi*xi_i(i))*(1+zeta*zeta_i(i));
                    dN_dzeta(i)=1/4 * zeta_i(i)*(1+xi*xi_i(i))*(1-eta^2);
                else
                    N(i) = 1/4 * (1+xi*xi_i(i))*(1+eta*eta_i(i))*(1-zeta^2);
                    dN_dxi(i) = 1/4 * xi_i(i)*(1+eta*eta_i(i))*(1-zeta^2);
                    dN_deta(i)= 1/4 * eta_i(i)*(1+xi*xi_i(i))*(1-zeta^2);
                    dN_dzeta(i)= -1/2 * zeta * (1+xi*xi_i(i))*(1+eta*eta_i(i));
                end
            end
            
        case 'TET4'
            N = [1-xi-eta-zeta, xi, eta, zeta];
            dN_dxi = [-1, 1, 0, 0]; dN_deta = [-1, 0, 1, 0]; dN_dzeta = [-1, 0, 0, 1];
            
        case 'TET10'
            L = [1-xi-eta-zeta, xi, eta, zeta];
            dL_dxi = [-1, 1, 0, 0]; dL_deta = [-1, 0, 1, 0]; dL_dzeta = [-1, 0, 0, 1];
            
            N = zeros(1,10); dN_dxi = zeros(1,10); dN_deta = zeros(1,10); dN_dzeta = zeros(1,10);
            for i=1:4 % Corners
                N(i) = L(i)*(2*L(i)-1);
                dN_dxi(i) = (4*L(i)-1)*dL_dxi(i); dN_deta(i) = (4*L(i)-1)*dL_deta(i); dN_dzeta(i) = (4*L(i)-1)*dL_dzeta(i);
            end
            mids = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4]; % Abaqus Tet10 Midside Pairs
            for i=1:6
                m=i+4; p1=mids(i,1); p2=mids(i,2);
                N(m) = 4*L(p1)*L(p2);
                dN_dxi(m)   = 4*(dL_dxi(p1)*L(p2)   + L(p1)*dL_dxi(p2));
                dN_deta(m)  = 4*(dL_deta(p1)*L(p2)  + L(p1)*dL_deta(p2));
                dN_dzeta(m) = 4*(dL_dzeta(p1)*L(p2) + L(p1)*dL_dzeta(p2));
            end

        case 'WEDGE6'
            L1 = 1-xi-eta; L2 = xi; L3 = eta;
            Nt = [L1, L2, L3];
            dNt_dxi = [-1, 1, 0]; dNt_deta = [-1, 0, 1];
            Nz = [0.5*(1-zeta), 0.5*(1+zeta)]; dNz_dzeta = [-0.5, 0.5];
            
            N = [Nt*Nz(1), Nt*Nz(2)];
            dN_dxi   = [dNt_dxi*Nz(1), dNt_dxi*Nz(2)];
            dN_deta  = [dNt_deta*Nz(1), dNt_deta*Nz(2)];
            dN_dzeta = [Nt*dNz_dzeta(1), Nt*dNz_dzeta(2)];
    end
end
