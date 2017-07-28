function itl = associate_itl(itl,param,horizon)

N = size(itl,2);

if nargin == 3
    hormin = horizon(1);
    hormax = horizon(2);
else
    hormin = -inf;
    hormax = +inf;
end


% Get the active tracklets in the horizon
itlh = get_itl_horizon(itl,horizon);

% eliminate very small itls
len = [itlh.length];
itlh(len<=2) = [];

if ~isempty(itlh)
    % Do merging in the horizon
    dN = 1;
    
    while dN > 0
        if param.debug
            figure(11)
            drawitl(itlh)
        end
        
        
        % Compute similarities
        [S,itlh] = compute_itl_similarity_matrix(itlh,param.similarity_method,param.eta_max);
        
        
        % Compute associations
        % 1) Use minus similarity, because lapjv solves minimization problem.
        assign = generalizedLinearAssignment(-S,-param.min_s);

        
        % Do associations
        itlh = process_itl_associations(itlh,assign,param);
        
        Nnew = size(itlh,2);
        dN = N - Nnew;
        N = Nnew;
        
        if param.debug
            figure(11)
            drawitl(itlh)
        end
        
    end
    
    % merge with real tracklets
    itl = set_itl_horizon(itl,itlh,horizon);
    

end