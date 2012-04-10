function X = pinvs(A,tol)
%PINVS   Pseudoinverse (sparse) utilisant une décomposition en valeurs singulières.
%   X = PINV(A) produces a matrix X of the same dimensions
%   as A' so that A*X*A = A, X*A*X = X and A*X and X*A
%   are Hermitian. The computation is based on SVD(A) and any
%   singular values less than a tolerance are treated as zero.
%   The default tolerance is MAX(SIZE(A)) * NORM(A) * EPS.
%
%   PINVS(A,TOL) uses the tolerance TOL instead of the default.
%
%   on the basis of the of the function pinv
	
% Thermique 1.0 - 30/04/01 - Olivier Vitrac


if ~issparse(A), A = sparse(A); end

[m,n] = size(A);
[U,S,V] = svds(A,n);
if m > 1, s = nonzeros(S);
   elseif m == 1, s = S(1);
   else s = 0;
end
if nargin < 2
   tol = max(m,n) * max(s) * eps;
end
r = sum(s > tol);
if (r == 0)
   X = zeros(size(A'));
else
   s = diag(ones(r,1)./s(1:r));
   X = V(:,1:r)*s*U(:,1:r)';
end