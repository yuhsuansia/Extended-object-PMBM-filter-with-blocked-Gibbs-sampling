function [gospa, gospa_decomp] = computeGOSPAOverTime(data, est, p, c, alpha)
%COMPUTEGOSPAOVERTIME Compute extended-target GOSPA for all time steps.
%
% Inputs:
%   data.truthTracks  : 4 x T x N
%   data.truthExtents : 2 x 2 x T x N
%   est               : cell array (length T), est{k} is struct array of estimates
%                       with fields .mean (4x1) and .extentMean (2x2)
%   p, c, alpha       : GOSPA parameters
%
% Outputs:
%   gospa        : T x 1
%   gospa_decomp : 1 x T struct array with fields localisation, missed, false

T = numel(est);
gospa = zeros(T,1);
gospa_decomp(1,T) = struct('localisation',0,'missed',0,'false',0);

for k = 1:T
    % -------- build x_mat (truth) --------
    valid = ~isnan(squeeze(data.truthTracks(1,k,:)));
    idx = find(valid);

    if isempty(idx)
        x_mat = struct('x', zeros(2,0), 'X', zeros(2,2,0));
    else
        N = numel(idx);
        x = zeros(2,N);
        X = zeros(2,2,N);
        for n = 1:N
            o = idx(n);
            x(:,n) = data.truthTracks(1:2,k,o);
            X(:,:,n) = data.truthExtents(:,:,k,o);
        end
        x_mat = struct('x', x, 'X', X);
    end

    % -------- build y_mat (estimates) --------
    estk = est{k};
    if isempty(estk)
        y_mat = struct('x', zeros(2,0), 'X', zeros(2,2,0));
    else
        M = numel(estk);
        y = zeros(2,M);
        Y = zeros(2,2,M);
        for m = 1:M
            y(:,m) = estk(m).mean(1:2);
            Y(:,:,m) = estk(m).extentMean;
        end
        y_mat = struct('x', y, 'X', Y);
    end

    % -------- GOSPA --------
    [gospa(k), ~, gospa_decomp(k)] = GOSPA_extended(x_mat, y_mat, p, c, alpha);
end
end