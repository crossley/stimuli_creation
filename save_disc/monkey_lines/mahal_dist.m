function D = mahal_dist(point, mu, covar)

% Computes Mahalanobis distance for a point in a bivariate normal
% distribution where both the point and the mean are row vectors.

D = sqrt((point-mu)*inv(covar)*(point-mu)');