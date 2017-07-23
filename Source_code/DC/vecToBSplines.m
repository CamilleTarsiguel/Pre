function splines=vecToBSplines(stateVec, spl)

% Update splines coefficients from state vector
%


% 1. Assume all splines are cubic (order 3)
checkSplines(spl);

% 2. Assume number of splines in spl is correct


% How many splines?
N=length(spl);

% copy to new structures
splines = spl;

os=1;
for id=1:N
    n=spl(id).number;
    l=n*2;
    splines(id).coefs = ...
        reshape(stateVec(os:os+l-1),n,2)';
    os=os+l;
end
end
