function [Moment_Val, Integral_Val, Total_Vol] = Mat_3B_FEM_3D_GaussPoint_Torque(X, Y, Z, Field_GP, Center, Connectivity, Node_Labels)
% Integrates the Torque (or First Moment) of a Gauss Point field over a 3D domain.
%
% INPUTS:
%  X, Y, Z      : [N x 1] Nodal Coordinates
%  Field_GP     : [N_IP x 1] (Scalar Field) OR [N_IP x 3] (Vector Field)
%  Center       : [1 x 3] Coordinates of the rotation center [Xc, Yc, Zc]
%  Connectivity : [E x M] Element Matrix
%  Node_Labels  : [N x 1] User Node IDs
%
% OUTPUTS:
%  Moment_Val   : [1 x 3] The [Mx, My, Mz] Torque or Bending Moments
%  Integral_Val : [1 x N] Standard Volume Integral of the field (0th Moment)
%  Total_Vol    : Total RVE Volume

    % 1. Map Node Labels to Indices
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    % 2. Check Input Size
    nGP_Input = size(Field_GP, 1);
    is_vector_field = (size(Field_GP, 2) == 3);
    nElem = size(Connectivity, 1);
    avg_ips = round(nGP_Input / nElem);

    % 3. Define 3D Gauss Rules
    GaussRules = build_gauss_rules();

    % 4. Initialize Accumulators
    Moment_Val = zeros(1, 3);
    Integral_Val = zeros(1, size(Field_GP, 2));
    Total_Vol = 0.0;
    global_ip_idx = 0; 

    % Center Coordinates
    Xc = Center(1); Yc = Center(2); Zc = Center(3);

    % 5. Loop over Elements
    [~, totalCols] = size(Connectivity);
    
    for e = 1:nElem
        if totalCols < 2, continue; end
        
        TypeID = Connectivity(e, 2); % Usually the number of nodes
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
                continue; % Skip 1D/2D elements
        end
        
        if totalCols < (2 + numNodes), continue; end
        
        nodes = Connectivity(e, 3:2+numNodes);
        idx = Map(round(nodes));
        if any(idx == 0) || any(isnan(nodes)), continue; end
        
        xe = X(idx); ye = Y(idx); ze = Z(idx);
        numIPs = length(rule.xi);
        
        for ip = 1:numIPs
            global_ip_idx = global_ip_idx + 1;
            
            xi = rule.xi(ip); eta = rule.eta(ip); zeta = rule.zeta(ip); w = rule.w(ip);
            
            % 1. Get Shape Functions and compute global IP coordinates
            [N, dN_dxi, dN_deta, dN_dzeta] = get_shape_funcs_3D(xi, eta, zeta, mType);
            
            x_ip = N * xe;
            y_ip = N * ye;
            z_ip = N * ze;
            
            % 2. Calculate Lever Arm vector (r)
            rx = x_ip - Xc;
            ry = y_ip - Yc;
            rz = z_ip - Zc;
            
            % 3. Calculate 3D Volume Jacobian (dV)
            J = [dN_dxi * xe,   dN_dxi * ye,   dN_dxi * ze;
                 dN_deta * xe,  dN_deta * ye,  dN_deta * ze;
                 dN_dzeta * xe, dN_dzeta * ye, dN_dzeta * ze];
            dVol = abs(det(J)) * w;
            
            % 4. Retrieve Field Value
            val = Field_GP(global_ip_idx, :);
            
            % 5. Compute Moment / Torque
            if is_vector_field
                % True Torque: cross(r, F)
                dT_x = (ry * val(3) - rz * val(2)) * dVol;
                dT_y = (rz * val(1) - rx * val(3)) * dVol;
                dT_z = (rx * val(2) - ry * val(1)) * dVol;
                Moment_Val = Moment_Val + [dT_x, dT_y, dT_z];
            else
                % Scalar First Moment: r * Field (Used for bending moments)
                dT_x = rx * val * dVol;
                dT_y = ry * val * dVol;
                dT_z = rz * val * dVol;
                Moment_Val = Moment_Val + [dT_x, dT_y, dT_z];
            end
            
            % Standard Volumetric Integration
            Integral_Val = Integral_Val + val * dVol;
            Total_Vol    = Total_Vol + dVol;
        end
    end
end

%% Helper: Build Gauss Quadrature Rules
function GR = build_gauss_rules()
    GR = struct();
    % HEX8 (1 IP)
    GR.HEX8_1.xi = 0; GR.HEX8_1.eta = 0; GR.HEX8_1.zeta = 0; GR.HEX8_1.w = 8;
    
    % HEX8 (8 IPs)
    pt = 1/sqrt(3);
    GR.HEX8_8.xi   = [-pt, pt, -pt, pt, -pt, pt, -pt, pt];
    GR.HEX8_8.eta  = [-pt, -pt, pt, pt, -pt, -pt, pt, pt];
    GR.HEX8_8.zeta = [-pt, -pt, -pt, -pt, pt, pt, pt, pt];
    GR.HEX8_8.w    = ones(1,8);

    % HEX20 (8 IPs)
    GR.HEX20_8 = GR.HEX8_8;

    % HEX20 (27 IPs)
    pts3 = [-sqrt(3/5), 0, sqrt(3/5)];
    wts3 = [5/9, 8/9, 5/9];
    [Xg,Yg,Zg] = ndgrid(pts3, pts3, pts3);
    [WX,WY,WZ] = ndgrid(wts3, wts3, wts3);
    GR.HEX20_27.xi = Xg(:)'; GR.HEX20_27.eta = Yg(:)'; GR.HEX20_27.zeta = Zg(:)';
    W_mat = WX .* WY .* WZ; GR.HEX20_27.w = W_mat(:)';

    % TET4 (1 IP)
    GR.TET4_1.xi = 1/4; GR.TET4_1.eta = 1/4; GR.TET4_1.zeta = 1/4; GR.TET4_1.w = 1/6;

    % TET10 (4 IPs)
    a = (5+3*sqrt(5))/20; b = (5-sqrt(5))/20;
    GR.TET10_4.xi   = [b, a, b, b]; GR.TET10_4.eta  = [b, b, a, b];
    GR.TET10_4.zeta = [b, b, b, a]; GR.TET10_4.w    = ones(1,4) * (1/24);

    % WEDGE6 (2 IPs)
    GR.WEDGE6_2.xi   = [1/3, 1/3]; GR.WEDGE6_2.eta  = [1/3, 1/3];
    GR.WEDGE6_2.zeta = [-pt, pt];  GR.WEDGE6_2.w    = [0.5, 0.5];
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
            
        case 'HEX20'
            xi_i   = [-1, 1, 1,-1,-1, 1, 1,-1,  0, 1, 0,-1,  0, 1, 0,-1, -1, 1, 1,-1];
            eta_i  = [-1,-1, 1, 1,-1,-1, 1, 1, -1, 0, 1, 0, -1, 0, 1, 0, -1,-1, 1, 1];
            zeta_i = [-1,-1,-1,-1, 1, 1, 1, 1, -1,-1,-1,-1,  1, 1, 1, 1,  0, 0, 0, 0];
            N = zeros(1,20); dN_dxi = zeros(1,20); dN_deta = zeros(1,20); dN_dzeta = zeros(1,20);
            for i=1:8
                N(i) = 1/8 * (1+xi*xi_i(i))*(1+eta*eta_i(i))*(1+zeta*zeta_i(i)) * (xi*xi_i(i) + eta*eta_i(i) + zeta*zeta_i(i) - 2);
                dN_dxi(i) = 1/8 * xi_i(i)*(1+eta*eta_i(i))*(1+zeta*zeta_i(i)) * (2*xi*xi_i(i) + eta*eta_i(i) + zeta*zeta_i(i) - 1);
                dN_deta(i)= 1/8 * eta_i(i)*(1+xi*xi_i(i))*(1+zeta*zeta_i(i)) * (xi*xi_i(i) + 2*eta*eta_i(i) + zeta*zeta_i(i) - 1);
                dN_dzeta(i)=1/8 * zeta_i(i)*(1+xi*xi_i(i))*(1+eta*eta_i(i)) * (xi*xi_i(i) + eta*eta_i(i) + 2*zeta*zeta_i(i) - 1);
            end
            for i=9:20
                if xi_i(i) == 0
                    N(i) = 1/4 * (1-xi^2)*(1+eta*eta_i(i))*(1+zeta*zeta_i(i));
                    dN_dxi(i) = -1/2 * xi * (1+eta*eta_i(i))*(1+zeta*zeta_i(i)); dN_deta(i)= 1/4 * eta_i(i)*(1-xi^2)*(1+zeta*zeta_i(i)); dN_dzeta(i)=1/4 * zeta_i(i)*(1-xi^2)*(1+eta*eta_i(i));
                elseif eta_i(i) == 0
                    N(i) = 1/4 * (1+xi*xi_i(i))*(1-eta^2)*(1+zeta*zeta_i(i));
                    dN_dxi(i) = 1/4 * xi_i(i)*(1-eta^2)*(1+zeta*zeta_i(i)); dN_deta(i)= -1/2 * eta * (1+xi*xi_i(i))*(1+zeta*zeta_i(i)); dN_dzeta(i)=1/4 * zeta_i(i)*(1+xi*xi_i(i))*(1-eta^2);
                else
                    N(i) = 1/4 * (1+xi*xi_i(i))*(1+eta*eta_i(i))*(1-zeta^2);
                    dN_dxi(i) = 1/4 * xi_i(i)*(1+eta*eta_i(i))*(1-zeta^2); dN_deta(i)= 1/4 * eta_i(i)*(1+xi*xi_i(i))*(1-zeta^2); dN_dzeta(i)= -1/2 * zeta * (1+xi*xi_i(i))*(1+eta*eta_i(i));
                end
            end
            
        case 'TET4'
            N = [1-xi-eta-zeta, xi, eta, zeta];
            dN_dxi = [-1, 1, 0, 0]; dN_deta = [-1, 0, 1, 0]; dN_dzeta = [-1, 0, 0, 1];
            
        case 'TET10'
            L = [1-xi-eta-zeta, xi, eta, zeta];
            dL_dxi = [-1, 1, 0, 0]; dL_deta = [-1, 0, 1, 0]; dL_dzeta = [-1, 0, 0, 1];
            N = zeros(1,10); dN_dxi = zeros(1,10); dN_deta = zeros(1,10); dN_dzeta = zeros(1,10);
            for i=1:4 
                N(i) = L(i)*(2*L(i)-1);
                dN_dxi(i) = (4*L(i)-1)*dL_dxi(i); dN_deta(i) = (4*L(i)-1)*dL_deta(i); dN_dzeta(i) = (4*L(i)-1)*dL_dzeta(i);
            end
            mids = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4];
            for i=1:6
                m=i+4; p1=mids(i,1); p2=mids(i,2);
                N(m) = 4*L(p1)*L(p2);
                dN_dxi(m)   = 4*(dL_dxi(p1)*L(p2)   + L(p1)*dL_dxi(p2));
                dN_deta(m)  = 4*(dL_deta(p1)*L(p2)  + L(p1)*dL_deta(p2));
                dN_dzeta(m) = 4*(dL_dzeta(p1)*L(p2) + L(p1)*dL_dzeta(p2));
            end

        case 'WEDGE6'
            L1 = 1-xi-eta; L2 = xi; L3 = eta;
            Nt = [L1, L2, L3]; dNt_dxi = [-1, 1, 0]; dNt_deta = [-1, 0, 1];
            Nz = [0.5*(1-zeta), 0.5*(1+zeta)]; dNz_dzeta = [-0.5, 0.5];
            N = [Nt*Nz(1), Nt*Nz(2)];
            dN_dxi   = [dNt_dxi*Nz(1), dNt_dxi*Nz(2)];
            dN_deta  = [dNt_deta*Nz(1), dNt_deta*Nz(2)];
            dN_dzeta = [Nt*dNz_dzeta(1), Nt*dNz_dzeta(2)];
    end
end
