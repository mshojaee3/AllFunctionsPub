function [gradXgp, gradYgp, gradXnp, gradYnp] = FEM_2D_nodal_gradient_12_2_26(X, Y, Field, Connectivity, Node_Labels)
% FEM_2D_NODAL_GRADIENT_11_2_26
% Mixed-mesh compatible gradient calculation at:
%   (1) Gauss/Integration points (per element)
%   (2) Nodes (averaged from element-wise nodal evaluations)
%
% INPUTS:
%  X, Y         : [N x 1] nodal coordinates
%  Field        : [N x 1] nodal scalar field values
%  Connectivity : [E x M] [ElementID, TypeID, N1, N2, ...]
%                TypeID mapping: 3=T3, 6=T6, 4=Q4, 8=Q8
%  Node_Labels  : [N x 1] Abaqus node labels (may be non-sequential)
%
% OUTPUTS:
%  gradXgp : [nIPtotal x 3] [ElementID, IP, dField/dX] at Gauss points
%  gradYgp : [nIPtotal x 3] [ElementID, IP, dField/dY] at Gauss points
%  gradXnp : [N x 1] nodal gradient dField/dX (averaged)
%  gradYnp : [N x 1] nodal gradient dField/dY (averaged)

    % ----------------------------
    % Map Abaqus node labels -> 1..N indices
    % ----------------------------
    max_lbl = max(Node_Labels);
    Map = zeros(max_lbl, 1);
    Map(Node_Labels) = 1:length(Node_Labels);

    nN = length(X);

    % ----------------------------
    % Accumulators for NODAL gradients
    % ----------------------------
    AccX = zeros(nN, 1);
    AccY = zeros(nN, 1);
    Cnt  = zeros(nN, 1);

    % ----------------------------
    % Parent-space coordinates of element NODES (for nodal gradient evaluation)
    % ----------------------------
    localCoords = struct();
    localCoords.T3 = struct('x', [0, 1, 0], 'y', [0, 0, 1]);
    localCoords.T6 = struct('x', [0, 1, 0, 0.5, 0.5, 0], 'y', [0, 0, 1, 0, 0.5, 0.5]);
    localCoords.Q4 = struct('x', [-1, 1, 1, -1], 'y', [-1, -1, 1, 1]);
    localCoords.Q8 = struct('x', [-1, 1, 1, -1, 0, 1, 0, -1], 'y', [-1, -1, 1, 1, -1, 0, 1, 0]);

    % ----------------------------
    % Gauss points in parent space for each element type
    % ----------------------------
    g = 1/sqrt(3);

    gp = struct();
    % T3: 1-point rule at centroid (xi,eta) in your triangle parent coords
    gp.T3 = struct('xi', 1/3, 'eta', 1/3);
    % T6: 3-point rule in barycentric coordinates mapped to (xi=L1, eta=L2)
    % (L1,L2,L3) = (1/6,1/6,2/3), (1/6,2/3,1/6), (2/3,1/6,1/6)
    gp.T6 = struct('xi',  [1/6, 1/6, 2/3], ...
                   'eta', [1/6, 2/3, 1/6]);
    % Q4/Q8: 2x2 Gauss points
    gp.Q4 = struct('xi',  [-g, +g, +g, -g], ...
                   'eta', [-g, -g, +g, +g]);
    gp.Q8 = gp.Q4;

    % ----------------------------
    % Storage for Gauss-point gradients
    % We'll collect rows then vertcat at the end.
    % ----------------------------
    gradXgp_rows = cell(size(Connectivity,1), 1);
    gradYgp_rows = cell(size(Connectivity,1), 1);

    [~, totalCols] = size(Connectivity);

    for e = 1:size(Connectivity, 1)

        if totalCols < 2
            continue;
        end

        elID  = Connectivity(e, 1);
        typeID = Connectivity(e, 2);

        switch typeID
            case 3
                mType = 'T3'; numNodes = 3; numGP = 1;
            case 6
                mType = 'T6'; numNodes = 6; numGP = 3;
            case 4
                mType = 'Q4'; numNodes = 4; numGP = 4;
            case 8
                mType = 'Q8'; numNodes = 8; numGP = 4;
            otherwise
                continue;
        end

        if totalCols < (2 + numNodes)
            continue;
        end

        nodes = Connectivity(e, 3:2+numNodes);

        if any(isnan(nodes))
            continue;
        end

        idx = Map(round(nodes));
        if any(idx == 0)
            continue;
        end

        xe = X(idx);
        ye = Y(idx);
        fe = Field(idx);

        % ----------------------------
        % (A) GAUSS-POINT gradients
        % ----------------------------
        xi_list  = gp.(mType).xi;
        eta_list = gp.(mType).eta;

        % Make sure scalar GP definitions expand correctly
        if numGP == 1
            xi_list = xi_list(:);
            eta_list = eta_list(:);
        end

        xrows = zeros(numGP, 3);
        yrows = zeros(numGP, 3);

        for ip = 1:numGP
            xi  = xi_list(ip);
            eta = eta_list(ip);

            [dN_dxi, dN_deta] = get_shape_derivatives(xi, eta, mType);

            % Jacobian J = [dx/dxi  dy/dxi; dx/deta  dy/deta]
            J = [dN_dxi * xe,  dN_dxi * ye; ...
                 dN_deta * xe, dN_deta * ye];

            % Global derivatives: [dN/dx; dN/dy] = inv(J)*[dN/dxi; dN/deta]
            gN = J \ [dN_dxi; dN_deta];   % 2 x numNodes

            dfdx = gN(1,:) * fe;
            dfdy = gN(2,:) * fe;

            xrows(ip,:) = [elID, ip, dfdx];
            yrows(ip,:) = [elID, ip, dfdy];
        end

        gradXgp_rows{e} = xrows;
        gradYgp_rows{e} = yrows;

        % ----------------------------
        % (B) NODAL gradients (your existing approach)
        % Evaluate gradient at each element node location in parent space
        % ----------------------------
        xN = localCoords.(mType).x;
        yN = localCoords.(mType).y;

        for n = 1:numNodes
            [dN_dxi, dN_deta] = get_shape_derivatives(xN(n), yN(n), mType);

            J = [dN_dxi * xe,  dN_dxi * ye; ...
                 dN_deta * xe, dN_deta * ye];

            gN = J \ [dN_dxi; dN_deta];

            p = idx(n);
            AccX(p) = AccX(p) + gN(1,:) * fe;
            AccY(p) = AccY(p) + gN(2,:) * fe;
            Cnt(p)  = Cnt(p)  + 1;
        end
    end

    % Concatenate Gauss results
    gradXgp = vertcat(gradXgp_rows{:});
    gradYgp = vertcat(gradYgp_rows{:});

    % Average nodal results
    gradXnp = AccX ./ max(Cnt, 1);
    gradYnp = AccY ./ max(Cnt, 1);
end


%% Helper: Shape Function Derivatives
function [dN_dxi, dN_deta] = get_shape_derivatives(xi, eta, type)
    switch upper(type)
        case 'T3' % Linear Triangle (parent coords: (0,0),(1,0),(0,1))
            % N1=xi, N2=eta, N3=1-xi-eta (or another convention)
            % Your original derivatives correspond to that convention:
            dN_dxi  = [1, 0, -1];
            dN_deta = [0, 1, -1];

        case 'T6' % Quadratic Triangle with L1=xi, L2=eta, L3=1-xi-eta
            L1 = xi; L2 = eta; L3 = 1 - xi - eta;
            dN_dxi  = [4*L1-1, 0, 1-4*L3, 4*L2, -4*L2, 4*(L3-L1)];
            dN_deta = [0, 4*L2-1, 1-4*L3, -4*L1, 4*L1, 4*(L3-L2)];

        case 'Q4' % Linear Quad
            dN_dxi  = 0.25 * [-(1-eta),  (1-eta), (1+eta), -(1+eta)];
            dN_deta = 0.25 * [-(1-xi),  -(1+xi),  (1+xi),   (1-xi)];

        case 'Q8' % Quadratic Quad (Serendipity)
            [dN_dxi, dN_deta] = sd_Q8(xi, eta);

        otherwise
            error('Unsupported element type: %s', type);
    end
end


%% Helper: Serendipity Quadratic Quad Derivatives (Q8)
function [dN_dxi, dN_deta] = sd_Q8(xi, eta)
    dN_dxi = 0.25 * [
        (1-eta)*(2*xi+eta), (1-eta)*(2*xi-eta), (1+eta)*(2*xi+eta), (1+eta)*(2*xi-eta), ...
        -4*xi*(1-eta), 2*(1-eta^2), -4*xi*(1+eta), -2*(1-eta^2) ];
    dN_deta = 0.25 * [
        (1-xi)*(2*eta+xi), (1+xi)*(2*eta-xi), (1+xi)*(2*eta+xi), (1-xi)*(2*eta-xi), ...
        -2*(1-xi^2), -4*eta*(1+xi), 2*(1-xi^2), -4*eta*(1-xi) ];
end
