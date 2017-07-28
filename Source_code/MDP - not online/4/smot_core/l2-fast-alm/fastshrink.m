function A=fastshrink(D,th)

[U S V]=svd(D);

[M N] = size(S);
idx = (1:M+1:M*N);

% TODO: make this shrinkage with vector (faster)
% QUESTION: Isn't it better to reduced multiplication like
% U(:,1:r)*S(1:r,1:r)*V(:,1:r)'

% S(diag(vec))=max(0,S(diag(vec))-th);
S(idx) = max(0,S(idx)-th);
A=U*S*V';

