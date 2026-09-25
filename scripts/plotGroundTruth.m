function plotGroundTruth(data, par)
%PLOTGROUNDTRUTH Plot trajectories and one-standard-deviation extent ellipses.
figure('Name','Ground truth'); hold on; grid on;

theta  = linspace(0, 2*pi, 201);
circle = [cos(theta); sin(theta)];

numObjects  = data.numObjects;
appearance  = data.appearance;
truthTracks = data.truthTracks;
truthExtents= data.truthExtents;

colors = lines(numObjects);
for obj = 1:numObjects
    k0 = max(1, appearance(1, obj));
    k1 = min(data.numSteps, appearance(2, obj));
    if k0 > k1, continue; end

    xy = truthTracks(1:2, k0:k1, obj);

    plot(xy(1,1), xy(2,1), 'x', 'Color', colors(obj,:), 'MarkerSize', 10, 'LineWidth', 1.5);

    E0 = truthExtents(:, :, k0, obj);
    contour = chol(E0, 'lower') * circle;
    plot(contour(1,:) + xy(1,1), contour(2,:) + xy(2,1),  '-.', 'Color', colors(obj,:), 'LineWidth', 1.5);

    plot(xy(1,:), xy(2,:), 'Color', colors(obj,:), 'LineWidth', 1.5);
end
axis equal;
xlim(par.region(:,1));
ylim(par.region(:,2));

xticks([-150 -100 -50 0 50 100 150])
yticks([-150 -100 -50 0 50 100 150])

xlabel('x (m)', 'Interpreter','latex'); ylabel('y (m)', 'Interpreter','latex');

set(gca,'TickLabelInterpreter','latex','FontSize', 16)
end
