function tests = testGGIWExtent
% Independent regression checks for the centroid contribution to the extent.
tests = functiontests(localfunctions);
end

function setupOnce(testCase)
testCase.TestData.originalPath = path;
addpath(fileparts(fileparts(mfilename('fullpath'))));
setupPath;
end

function teardownOnce(testCase)
path(testCase.TestData.originalPath);
end

function testGeneralUpdateAgainstPrincipalRootReference(testCase)
[g, par, Z] = fixture();
for count = [1, 2, 5]
    [actual, ll] = ggiwUpdate(g, Z(:,1:count), par);
    [expected, expectedLL] = referenceUpdate(g, Z(:,1:count), par, false);
    verifyEqual(testCase, actual.V, expected.V, 'AbsTol', 1e-10);
    verifyEqual(testCase, actual.mean, expected.mean, 'AbsTol', 1e-12);
    verifyEqual(testCase, actual.cov, expected.cov, 'AbsTol', 1e-12);
    verifyEqual(testCase, ll, expectedLL, 'AbsTol', 1e-10);
    verifyEqual(testCase, actual.Xhat, expected.V/(expected.v-3), 'AbsTol', 1e-12);
    verifyEqual(testCase, actual.logdetV, ...
        log(det(expected.V + 1e-6*eye(2))), 'AbsTol', 1e-12);
end
end

function testSingletonUpdateAgainstPrincipalRootReference(testCase)
[g, par, Z] = fixture();
[actual, ll] = ggiwUpdate1(g, Z(:,1), par);
[expected, expectedLL] = referenceUpdate(g, Z(:,1), par, true);
verifyEqual(testCase, actual.V, expected.V, 'AbsTol', 1e-10);
verifyEqual(testCase, ll, expectedLL, 'AbsTol', 1e-10);
end

function testKnownPositionLimit(testCase)
[g, par] = fixture();
g.cov = zeros(4);
g = ggiwFinalizeCache(g);
innovation = [2; -3];
for count = [2, 5]
    Z = repmat(par.H*g.mean + innovation, 1, count);
    actual = ggiwUpdate(g, Z, par);
    % With H*P*H'=0 and S=Xhat/count, N=count*eps*eps', not count^2.
    verifyEqual(testCase, actual.V-g.V, count*(innovation*innovation.'), ...
        'AbsTol', 1e-6); % Allows the existing 1e-9 innovation regularizer.
end
end

function testBernoulliAndPPPPaths(testCase)
[g, par, Z] = fixture();
prior = struct('ggiw', g, 'existence', 0.7, 'logWeight', log(0.4));
ppp = rmfield(prior, 'existence');
second = ppp;
second.ggiw.mean = g.mean + [1; -2; 0; 0];
second.ggiw.V = 1.3*g.V;
second.ggiw = ggiwFinalizeCache(second.ggiw);
second.logWeight = log(0.2);
for count = [1, 5]
    cellZ = Z(:,1:count);
    [expected, ll] = referenceUpdate(g, cellZ, par, false);
    bern = bernoulliUpdate(prior, cellZ, par);
    verifyEqual(testCase, bern.ggiw.V, expected.V, 'AbsTol', 1e-10);
    verifyEqual(testCase, bern.logWeight, prior.logWeight+log(prior.existence)+ll, ...
        'AbsTol', 1e-10);
    [expected, ll] = referenceUpdate(g, cellZ, par, count == 1);
    [~, ll2] = referenceUpdate(second.ggiw, cellZ, par, count == 1);
    single = pppUpdate(ppp, cellZ, par);
    mixture = pppUpdate([ppp, second], cellZ, par);
    verifyEqual(testCase, single.ggiw.V, expected.V, 'AbsTol', 1e-10);
    targetMass = exp(ppp.logWeight+ll);
    mixtureMass = targetMass + exp(second.logWeight+ll2);
    clutter = (count == 1)*par.clutterIntensity;
    verifyEqual(testCase, single.logWeight, log(targetMass+clutter), 'AbsTol', 1e-10);
    verifyEqual(testCase, mixture.logWeight, log(mixtureMass+clutter), 'AbsTol', 1e-10);
    verifyEqual(testCase, single.existence, targetMass/(targetMass+clutter), 'AbsTol', 1e-12);
    verifyEqual(testCase, mixture.existence, mixtureMass/(mixtureMass+clutter), 'AbsTol', 1e-12);
end
end

function [g, par, Z] = fixture()
par.H = [eye(2), zeros(2)];
par.clutterIntensity = 0.02;
g.alpha = 8;
g.beta = 2;
g.v = 12;
g.V = 9*[4, 1.2; 1.2, 2];
g.mean = [1; -1; 0.5; -0.2];
g.cov = blkdiag([3, -0.7; -0.7, 1], [0.5, 0.1; 0.1, 0.8]);
g = ggiwFinalizeCache(g);
% Anisotropic, noncommuting extent/innovation covariances expose root-order errors.
Z = [4, 5, 3, 6, 2; -3, -1, -4, -2, -5];
end

function [up, ll] = referenceUpdate(g, Z, par, singleton)
count = size(Z, 2);
centroid = mean(Z, 2);
residuals = Z-centroid;
innovation = centroid-par.H*g.mean;
S = par.H*g.cov*par.H.' + g.Xhat/count + 1e-9*eye(2);
% Deliberately use general matrix roots, independent of the optimized helper.
A = sqrtm(g.Xhat)/sqrtm(S);
e = A*innovation;
up = g;
up.V = g.V + e*e.' + residuals*residuals.';
up.v = g.v+count;
up.alpha = g.alpha+count;
up.beta = g.beta+1;
K = (g.cov*par.H.')/S;
up.mean = g.mean+K*innovation;
up.cov = g.cov-K*par.H*g.cov;
logdetV = log(det(up.V + singleton*1e-9*eye(2)));
multiGamma2 = @(a) 0.5*log(pi)+gammaln(a)+gammaln(a-0.5);
ll = -count*log(pi)-log(count) ...
    + g.v/2*g.logdetV-up.v/2*logdetV ...
    + multiGamma2(up.v/2)-multiGamma2(g.v/2) ...
    + 0.5*(g.logdetXhat-log(det(S))) ...
    + gammaln(up.alpha)-gammaln(g.alpha) ...
    + g.alpha*log(g.beta)-up.alpha*log(up.beta);
end
