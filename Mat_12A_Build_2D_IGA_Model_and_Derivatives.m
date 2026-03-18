function Out = Mat_12A_Build_2D_IGA_Model_and_Derivatives(Nu_x, p)
% =========================================================================
% BUILD_1D_IGA_DERIVATIVES
% =========================================================================
% Builds 1D IGA basis functions and their derivatives at Gauss points
% in the x-direction only.
% =========================================================================

    %% 1) Read input
    if nargin == 0
        Nu_x = 6;
        p    = 5;
    elseif nargin == 1
        if isempty(Nu_x)
            Nu_x = 6;
        end
        p = 5;
    elseif nargin == 2
        if isempty(Nu_x)
            Nu_x = 6;
        end
        if isempty(p)
            p = 5;
        end
    else
        error('Use Mat_12A_Build_2D_IGA_Model_and_Derivatives(), Mat_12A_Build_2D_IGA_Model_and_Derivatives(Nu_x), or Mat_12A_Build_2D_IGA_Model_and_Derivatives(Nu_x,p).');
    end

    epsilon = 1e-12;

    %% 2) Gauss points and weights on [-1,1]
    n_degree_gp = p + 1;

    if n_degree_gp == 1
        xgp0 = 0;
        wgp0 = 2;
    else
        beta = 0.5 ./ sqrt(1 - (2*(1:n_degree_gp-1)).^(-2));
        T = diag(beta,1) + diag(beta,-1);
        [V,D] = eig(T);
        xgp0 = diag(D);
        [xgp0, idx] = sort(xgp0);
        V = V(:,idx);
        wgp0 = 2*(V(1,:)').^2;
    end

    %% 3) Parametric element nodes in x
    xi_Node = linspace(0, 1, Nu_x + 1);

    %% 4) Map Gauss points to x-elements
    x_gauss_points   = zeros(Nu_x * n_degree_gp, 1);
    w_x_gauss_points = zeros(Nu_x * n_degree_gp, 1);

    for ix = 1:Nu_x
        id = (ix-1)*n_degree_gp + (1:n_degree_gp);

        x_min = xi_Node(ix);
        x_max = xi_Node(ix + 1);

        x_gauss_points(id)   = ((x_max - x_min)/2) * xgp0 + (x_max + x_min)/2;
        w_x_gauss_points(id) = wgp0;
    end

    desired_values = x_gauss_points;

    %% 5) Knot vector in x
    Xi = [zeros(1,p+1), (1:Nu_x-1)/Nu_x, ones(1,p+1)];
    n_xi = length(Xi) - p - 1;

    %% 6) Control-point distribution in x
    n = n_xi;

    if Nu_x < p + 3
        n = p + Nu_x;

        left_half = zeros(1, Nu_x);
        left_half(1) = 0;

        d_left = zeros(1, Nu_x - 1);
        for i = 1:(Nu_x - 1)
            d_left(i) = min(i, p);
        end
        left_half(2:Nu_x) = cumsum(d_left);

        d_right = zeros(1, p);
        for i = 1:p
            d_right(i) = min(Nu_x, p - i + 1);
        end
        right_half = left_half(end) + cumsum(d_right);

        V_hat = [left_half, right_half];
    else
        V_hat1 = zeros(1, n);
        V_hat2 = zeros(1, n);
        V_hat3 = zeros(1, n);

        V_hat1(2:n) = cumsum(1:n-1);

        V_hat2(1) = p * Nu_x;
        gh = 0;
        for i = 2:n
            gh = gh + i - 1;
            V_hat2(i) = p * Nu_x - gh;
        end
        V_hat2 = sort(V_hat2);

        k = 0;
        for i = p+2:n
            k = k + 1;
            V_hat3(i) = V_hat1(p+1) + k*p;
        end

        V_hat = [V_hat1(1:p+1), V_hat3(p+2:Nu_x-1), V_hat2(Nu_x:n)];
    end

    Ctrl_DimLess_X = V_hat / (p * Nu_x);

    %% 7) Store 1D IGA points
    IGA_Points.Gaussin.XGP = x_gauss_points;
    IGA_Points.Gaussin.WGPX = w_x_gauss_points;
    IGA_Points.Gaussin.n_degree_gp = n_degree_gp;

    IGA_Points.knot.DimLess_X = Xi;

    IGA_Points.Ctrl.DimLess_X = Ctrl_DimLess_X;
    IGA_Points.Ctrl.Number_X  = n_xi;

    IGA_Points.parametric.xi_parametric_values = linspace(0,1,n_xi);

    %% 8) Basis values in X at Gauss points
    N_values = zeros(p + 1, n_xi + 1, length(desired_values));

    for k = 1:length(desired_values)
        xi = desired_values(k);
        N = zeros(p+1, n_xi+1);   % keep n_xi+1 internally

        for i = 1:n_xi
            if i < n_xi
                if (xi >= Xi(i)) && (xi < Xi(i+1))
                    N(1, i) = 1;
                end
            else
                if (xi >= Xi(i)) && (xi <= Xi(i+1))
                    N(1, i) = 1;
                end
            end
        end

        if abs(xi - 1) < epsilon
            N(1, n_xi + 1) = 1;
        end

        for j = 2:p+1
            for i = 1:n_xi
                leftDenom  = Xi(i + j - 1) - Xi(i);
                rightDenom = Xi(i + j)     - Xi(i + 1);

                if abs(leftDenom) < epsilon
                    leftTerm = 0;
                else
                    leftTerm = ((xi - Xi(i)) / leftDenom) * N(j-1, i);
                end

                if abs(rightDenom) < epsilon
                    rightTerm = 0;
                else
                    rightTerm = ((Xi(i + j) - xi) / rightDenom) * N(j-1, i+1);
                end

                N(j, i) = leftTerm + rightTerm;
            end
        end

        N_values(:,:,k) = N;
    end

    %% 9) Extract degree-p basis
    Nlast_3D = N_values(p+1,:,:);
    Nlast = squeeze(Nlast_3D);

    if isvector(Nlast)
        Nlast = reshape(Nlast, [], 1);
    end

    N_last = Nlast(1:n_xi, :);   % final size: n_xi x nGP

    %% 10) First derivative
    d1N_values = zeros(n_xi, length(desired_values));

    if p >= 1
        for i = 1:n_xi
            denom1 = Xi(i + p)   - Xi(i);
            denom2 = Xi(i + p+1) - Xi(i + 1);

            c1 = 0;
            c2 = 0;
            if abs(denom1) > epsilon
                c1 = p / denom1;
            end
            if abs(denom2) > epsilon
                c2 = p / denom2;
            end

            for k = 1:length(desired_values)
                d1N_values(i,k) = c1 * N_values(p, i,   k) ...
                                - c2 * N_values(p, i+1, k);
            end
        end
    end

    %% 11) Second derivative
    d2N_values = zeros(n_xi, length(desired_values));

    if p >= 2
        prevd1N = zeros(n_xi+1, length(desired_values)); % internal helper
        p_minus_1 = p - 1;

        for i = 1:n_xi
            denom1 = Xi(i + p_minus_1)   - Xi(i);
            denom2 = Xi(i + p_minus_1+1) - Xi(i + 1);

            c1 = 0;
            c2 = 0;
            if abs(denom1) > epsilon
                c1 = p_minus_1 / denom1;
            end
            if abs(denom2) > epsilon
                c2 = p_minus_1 / denom2;
            end

            for k = 1:length(desired_values)
                prevd1N(i,k) = c1 * N_values(p_minus_1, i,   k) ...
                             - c2 * N_values(p_minus_1, i+1, k);
            end
        end

        for i = 1:n_xi
            denom1 = Xi(i + p)   - Xi(i);
            denom2 = Xi(i + p+1) - Xi(i + 1);

            c1 = 0;
            c2 = 0;
            if abs(denom1) > epsilon
                c1 = p / denom1;
            end
            if abs(denom2) > epsilon
                c2 = p / denom2;
            end

            for k = 1:length(desired_values)
                d2N_values(i,k) = c1 * prevd1N(i,   k) ...
                                - c2 * prevd1N(i+1, k);
            end
        end
    end

    %% 12) Third derivative
    d3N_values = zeros(n_xi, length(desired_values));

    if p >= 3
        p_minus_2 = p - 2;
        prevd1N_of_p_minus_2 = zeros(n_xi+1, length(desired_values)); % internal helper

        for i = 1:n_xi
            denom1 = Xi(i + p_minus_2)   - Xi(i);
            denom2 = Xi(i + p_minus_2+1) - Xi(i + 1);

            c1 = 0;
            c2 = 0;
            if abs(denom1) > epsilon
                c1 = p_minus_2 / denom1;
            end
            if abs(denom2) > epsilon
                c2 = p_minus_2 / denom2;
            end

            for k = 1:length(desired_values)
                prevd1N_of_p_minus_2(i,k) = ...
                    c1 * N_values(p_minus_2, i,   k) ...
                  - c2 * N_values(p_minus_2, i+1, k);
            end
        end

        p_minus_1 = p - 1;
        prevd2N_of_p_minus_1 = zeros(n_xi+1, length(desired_values)); % internal helper

        for i = 1:n_xi
            denom1 = Xi(i + p_minus_1)   - Xi(i);
            denom2 = Xi(i + p_minus_1+1) - Xi(i + 1);

            c1 = 0;
            c2 = 0;
            if abs(denom1) > epsilon
                c1 = p_minus_1 / denom1;
            end
            if abs(denom2) > epsilon
                c2 = p_minus_1 / denom2;
            end

            for k = 1:length(desired_values)
                prevd2N_of_p_minus_1(i,k) = ...
                    c1 * prevd1N_of_p_minus_2(i,   k) ...
                  - c2 * prevd1N_of_p_minus_2(i+1, k);
            end
        end

        for i = 1:n_xi
            denom1 = Xi(i + p)   - Xi(i);
            denom2 = Xi(i + p+1) - Xi(i + 1);

            c1 = 0;
            c2 = 0;
            if abs(denom1) > epsilon
                c1 = p / denom1;
            end
            if abs(denom2) > epsilon
                c2 = p / denom2;
            end

            for k = 1:length(desired_values)
                d3N_values(i,k) = ...
                    c1 * prevd2N_of_p_minus_1(i,   k) ...
                  - c2 * prevd2N_of_p_minus_1(i+1, k);
            end
        end
    end

    %% 13) Final output
    Out.IGA_Points  = IGA_Points;
    Out.IGA_parm.Nu_x = Nu_x;
    Out.IGA_parm.p    = p;
    Out.N_last      = N_last;
    Out.d1N_values  = d1N_values;
    Out.d2N_values  = d2N_values;
    Out.d3N_values  = d3N_values;
end
