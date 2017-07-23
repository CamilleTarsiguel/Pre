function stateVec=splinesToVec(spl)

% Convert splines structs to a state vector
%
% Take a vector of splines (pp structs)
% and extract the coefficients into one big vector.
% The order is Xn,...,X1,Cx,Yn,...,Y1,Cy,
% i.e. starting from the largest order to the constant.
% First all X coefficients, then all y coefficients.
%
% All splines must be cubic!
checkSplines(spl);

% How many splines do we have?
N=length(spl);

% preallocate space num pieces * 8
stateVec=zeros(sum([spl(:).pieces])*8,1);

os=1;
for id=1:N
    
    % extract coefficients and transpose to column vec
    xvec=spl(id).coefs'; xvec=xvec(:);
    l=length(xvec);             % current length
    stateVec(os:os+l-1)=xvec;   % assign
    os=os+l;                    % adjust offset
end


end