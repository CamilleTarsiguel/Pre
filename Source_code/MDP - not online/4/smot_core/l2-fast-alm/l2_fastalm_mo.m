% This code solves following robust PCA problem, using ALM method, with indexing
% matrices for fast operations. 

% min           norm(H,'nuclear') + 0.5*lamda*norm(P*h-P*d)^2
% subject to:   H = hankel(h)
% size(H)=[m n] , size(h) = [m+n-1 1], size(P) = [k m+n-1] 
% where k<=m+n-1. 

% P is intended to account for the possibility to choose a subset of h and
% require it to be close to data vector d, in l-2 norm sense. Hence, P 
% is derived from an eye(m+n-1) matrix, by deleting the entire row i if
% h(i) is not required to be close to corresponding data, or if there
% is not such data. 
% It may also be derived from any diagonal matrix, to account for weigthing
% of components of h in terms of penalizing non-closeness to data...etc.
% 
% Note: A generic P can be also used, but to keep the code fast and
% spesific, it is only allowed to have diagonal P, to give some flexibility
% which is commonly encountered. (like missing data, or fitting to future
% data which is absent...etc). 
% Also, for convergence, only the constraints are checked. Full check can
% be added, but this works nice in practice. 
% Author: Burak Yilmaz
% Last Update: 21 March 2013
% Change Log
% 02 April 2013 Caglayan Dicle added data range invariance
% 26 March 2013 Caglayan Dicle added multi-dimension support
% 26 March 2013 Caglayan Dicle changed namings for readibility
% 22 March 2013 Caglayan Dicle added options


function u_hat = l2_fastalm_mo(u,lambda,varargin)

[D N] = size(u);

% Default values
defOmega  = ones(1,N);
nr = ceil(N/(D+1))*D;
nc = N - ceil(N/(D+1))+1;
defJSize = [nr nc];
defMaxIter = 1000;
defTol = 1e-7;

% parse inputs
p = inputParser;
addParamValue(p,'omega',defOmega,@isnumeric);
addParamValue(p,'maxiter',defMaxIter,@isnumeric);
addParamValue(p,'tol',defTol,@isnumeric);
addParamValue(p,'size',defJSize,@isnumeric);
parse(p,varargin{:});

Omega = p.Results.omega;
Jsize = p.Results.size;
MaxIter = p.Results.maxiter;
Tol = p.Results.tol;

% % make u a vector
vec_u = u(:);

% form S
% Our matrix will be reshape(S*vec_u,[nr nc])
S = form_S(nr,nc,D);
% form P from omega
P = form_P(Omega,[D N]);


% coefficients (savings)
PtP = P'*P;
StS = S'*S;
diag_PtP = diag(PtP);
diag_StS = diag(StS);


% initializations
% mu = 0.05;
mu = 0.05 /( max(u(:))/10); % hack to adjust to data range
rho = 1.05;
h = 1.1*vec_u;
y = zeros(nr*nc,1);
R = zeros(Jsize);

for iter=1:MaxIter
    % j update
%     R = reshape(S*h+y/mu ,Jsize);     % slower
    R(:) = S*h+y/mu;                    % faster
    J = fastshrink(R, 1/mu);
    j = J(:);
    
    % h update
    % matrix form update
%     num = (lambda*PtP*vec_u + mu*S'*(j - y/mu));
%     denom = (lambda*PtP + mu*StS);    
%     h = denom\num;
    % faster update using structure of the problem
    h = (lambda*PtP*vec_u + mu*S'*(j-y/mu)) ./ (lambda*diag_PtP + mu*diag_StS);
    
    % lagrange updates
    y = y + mu*(S*h-j);
    
    % increase the mu so it will force the constraint
    mu = mu*rho;
    
    if norm(S*h-j) < Tol
        break;
    end    
    
end
u_hat = reshape(h,size(u));
% iter






function S = form_S(nr,nc,D)

% TODO: This functio may be faster. But it may not be necassary!
N = nr/D + nc-1;

S = sparse(nr*nc,N*D);
% s = sparse(nr,N*D);
% ind = 1:(nr+1):nr*nr;
% s(ind) = 1;

s = [eye(nr) zeros(nr,N*D-nr)];

for c=1:nc
    S((c-1)*nr+1:c*nr,:) = circshift(s,[0 (c-1)*D]);
end
35;

function P = form_P(Omega,DN)
D = DN(1);
N = DN(2);

P = sparse(N*D,N*D);
P(1:N*D+1:N^2*D^2) = reshape(repmat(Omega,[D 1]),[1 N*D]);