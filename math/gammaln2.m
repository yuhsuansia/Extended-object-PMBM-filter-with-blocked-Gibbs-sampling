function [gamma2log] = gammaln2(a)
%GAMMALN2 Evaluate the log multivariate gamma function in dimension two.
gamma2log = log(pi)/2 + gammaln(a) + gammaln(a - 1/2);
end
