function [d_pgospa, decomposed_cost] = ...
    PGOSPA(x_r, x_mean, x_cov, y_r, y_mean, y_cov, p, c, alpha)
%PGOSPA Probabilistic GOSPA distance between two multi-Bernoulli densities.
%
% This implementation follows the provided reference formulation and uses
% GaussianWassersteinDistance as the base distance.

n_output_arg = nargout;
checkInput();

nx = length(x_r);
ny = length(y_r);

decomposed_cost = struct( ...
    'localisation', 0, ...
    'existence_mismatch', 0, ...
    'missed', 0, ...
    'false', 0);

if nx == 0 || ny == 0
    if nx == 0
        decomposed_cost.false = sum(y_r) * c^p / alpha;
    end
    if ny == 0
        decomposed_cost.missed = sum(x_r) * c^p / alpha;
    end
    d_pgospa = (decomposed_cost.localisation + decomposed_cost.existence_mismatch + ...
                decomposed_cost.missed + decomposed_cost.false)^(1 / p);
    return
end

loc_cost_mat = zeros(nx, ny);
existence_cost_mat = zeros(nx, ny);

if alpha == 2
    cost_mat = inf(nx, ny + nx);
    for ix = 1:nx
        for iy = 1:ny
            loc_cost_mat(ix, iy) = min(x_r(ix), y_r(iy)) * ...
                computeBaseDistance(x_mean(:, ix), x_cov(:, :, ix), ...
                                    y_mean(:, iy), y_cov(:, :, iy))^p;
            existence_cost_mat(ix, iy) = abs(x_r(ix) - y_r(iy)) * c^p / 2;
            cost_mat(ix, iy) = existence_cost_mat(ix, iy) + loc_cost_mat(ix, iy);
        end
    end
    for j = 1:ny
        cost_mat(:, j) = cost_mat(:, j) - y_r(j) * c^p / 2;
    end
    for i = 1:nx
        cost_mat(i, i + ny) = x_r(i) * c^p / 2;
    end
else
    for ix = 1:nx
        for iy = 1:ny
            loc_cost_mat(ix, iy) = min(x_r(ix), y_r(iy)) * ...
                min(computeBaseDistance(x_mean(:, ix), x_cov(:, :, ix), ...
                                        y_mean(:, iy), y_cov(:, :, iy)), c)^p;
            existence_cost_mat(ix, iy) = abs(x_r(ix) - y_r(iy)) * c^p / alpha;
        end
    end
    cost_mat = loc_cost_mat + existence_cost_mat;
end

if alpha == 2
    [x_to_y_assignment, y_to_x_assignment] = assign2D(cost_mat);
    for i = 1:nx
        if x_to_y_assignment(i) <= ny
            decomposed_cost.localisation = ...
                decomposed_cost.localisation + loc_cost_mat(i, x_to_y_assignment(i));
            decomposed_cost.existence_mismatch = ...
                decomposed_cost.existence_mismatch + existence_cost_mat(i, x_to_y_assignment(i));
        end
        if x_to_y_assignment(i) > ny
            decomposed_cost.missed = decomposed_cost.missed + x_r(i) * c^p / 2;
        end
    end
    for i = 1:ny
        if y_to_x_assignment(i) == 0
            decomposed_cost.false = decomposed_cost.false + y_r(i) * c^p / 2;
        end
    end
    d_pgospa = (decomposed_cost.localisation + decomposed_cost.existence_mismatch + ...
                decomposed_cost.missed + decomposed_cost.false)^(1 / p);
else
    dummy_cost = (c^p) / alpha;
    opt_cost = 0;
    if nx == 0
        for i = 1:ny
            opt_cost = opt_cost + y_r(i) * dummy_cost;
        end
    elseif ny == 0
        for i = 1:nx
            opt_cost = opt_cost + x_r(i) * dummy_cost;
        end
    else
        [x_to_y_assignment, y_to_x_assignment] = assign2D(cost_mat);
        for ind = 1:nx
            if x_to_y_assignment(ind) ~= 0
                opt_cost = opt_cost + cost_mat(ind, x_to_y_assignment(ind));
            else
                opt_cost = opt_cost + x_r(ind) * dummy_cost;
            end
        end
        for ind = 1:ny
            if y_to_x_assignment(ind) == 0
                opt_cost = opt_cost + y_r(ind) * dummy_cost;
            end
        end
    end
    d_pgospa = opt_cost^(1 / p);
end

    function checkInput()
        if size(x_mean, 1) ~= size(y_mean, 1)
            error('The number of rows in x_mean and y_mean should be equal.');
        end
        if ~((p >= 1) && (p < inf))
            error('The value of exponent p should be within [1, inf).');
        end
        if ~(c > 0)
            error('The value of base distance c should be larger than 0.');
        end
        if ~((alpha > 0) && (alpha <= 2))
            error('The value of alpha should be within (0, 2].');
        end
        if alpha ~= 2 && n_output_arg == 2
            warning(['decomposed_cost is not valid for alpha = ', num2str(alpha)]);
        end
    end
end

function W = computeBaseDistance(mu1, Sigma1, mu2, Sigma2)
W = GaussianWassersteinDistance(mu1, Sigma1, mu2, Sigma2);
end
