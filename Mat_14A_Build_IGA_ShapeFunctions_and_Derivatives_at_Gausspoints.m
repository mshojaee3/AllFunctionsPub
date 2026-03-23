function IsoGeo = Mat_14A_Build_IGA_ShapeFunctions_and_Derivatives_at_Gausspoints(params)
% Mat_14A_Build_IGA_ShapeFunctions_and_Derivatives
% Builds 1D IGA shape functions and their derivatives at Gauss points.
%
% INPUT
%   params : struct with fields
%            params.Nu_x
%            params.p
%
% OUTPUT
%   IsoGeo : struct with fields
%            x_gp, w_gp, n_gp, Nu_x, n_xi, N, dN

Nu_x = params.Nu_x;
p    = params.p;

epsilon = 1e-12;
n_gp = p + 1;

% Gauss points on [-1,1]
if n_gp == 1
    xgp0 = 0;
    wgp0 = 2;
else
    beta = 0.5 ./ sqrt(1 - (2*(1:n_gp-1)).^(-2));
    T = diag(beta,1) + diag(beta,-1);
    [V,D] = eig(T);
    xgp0 = diag(D);
    [xgp0, idx] = sort(xgp0);
    V = V(:,idx);
    wgp0 = 2*(V(1,:)').^2;
end

% Element nodes
xi_node = linspace(0, 1, Nu_x + 1);

% Map GP to each element
x_gp = zeros(Nu_x * n_gp, 1);
w_gp = zeros(Nu_x * n_gp, 1);

for ex = 1:Nu_x
    id = (ex-1)*n_gp + (1:n_gp);
    a = xi_node(ex);
    b = xi_node(ex+1);

    x_gp(id) = ((b-a)/2)*xgp0 + (a+b)/2;
    w_gp(id) = wgp0;
end

Xi = [zeros(1,p+1), (1:Nu_x-1)/Nu_x, ones(1,p+1)];
n_xi = length(Xi) - p - 1;

N_values = zeros(p+1, n_xi+1, length(x_gp));

for k = 1:length(x_gp)
    xi = x_gp(k);
    Ntemp = zeros(p+1, n_xi+1);

    for i = 1:n_xi
        if i < n_xi
            if (xi >= Xi(i)) && (xi < Xi(i+1))
                Ntemp(1,i) = 1;
            end
        else
            if (xi >= Xi(i)) && (xi <= Xi(i+1))
                Ntemp(1,i) = 1;
            end
        end
    end

    if abs(xi - 1) < epsilon
        Ntemp(1,n_xi+1) = 1;
    end

    for j = 2:p+1
        for i = 1:n_xi
            leftDenom  = Xi(i+j-1) - Xi(i);
            rightDenom = Xi(i+j)   - Xi(i+1);

            if abs(leftDenom) < epsilon
                leftTerm = 0;
            else
                leftTerm = ((xi - Xi(i))/leftDenom) * Ntemp(j-1,i);
            end

            if abs(rightDenom) < epsilon
                rightTerm = 0;
            else
                rightTerm = ((Xi(i+j) - xi)/rightDenom) * Ntemp(j-1,i+1);
            end

            Ntemp(j,i) = leftTerm + rightTerm;
        end
    end

    N_values(:,:,k) = Ntemp;
end

Nlast = squeeze(N_values(p+1,:,:));
if isvector(Nlast)
    Nlast = reshape(Nlast, [], 1);
end
N = Nlast(1:n_xi,:);

dN = zeros(n_xi, length(x_gp));
if p >= 1
    for i = 1:n_xi
        denom1 = Xi(i+p)   - Xi(i);
        denom2 = Xi(i+p+1) - Xi(i+1);

        c1 = 0;
        c2 = 0;
        if abs(denom1) > epsilon
            c1 = p / denom1;
        end
        if abs(denom2) > epsilon
            c2 = p / denom2;
        end

        for k = 1:length(x_gp)
            dN(i,k) = c1 * N_values(p,i,k) - c2 * N_values(p,i+1,k);
        end
    end
end

IsoGeo.x_gp = x_gp;
IsoGeo.w_gp = w_gp;
IsoGeo.n_gp = n_gp;
IsoGeo.Nu_x = Nu_x;
IsoGeo.n_xi = n_xi;
IsoGeo.N    = N;
IsoGeo.dN   = dN;

end
