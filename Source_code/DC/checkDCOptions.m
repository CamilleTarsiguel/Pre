function checkDCOptions(opt)
% Perform a test to ensure that the opt struct is correct

%% Remarks
if opt.dataFunction==1 && (strcmpi(opt.conOpt.alg,'CGD') || strcmpi(opt.conOpt.alg,'fminunc'))
    warning('Absolute L2 distance will be approximated by the Pseudo-Huber loss for continuous minimization');
elseif  opt.dataFunction==3
    error('Lorentzian data loss not implemented for continuous case');
end

% if opt.readjustSE && ~opt.conOpt.enParEseg
%     warning('start/end points are beeing readjusted, but curvature not regularized');
% end
% if opt.fidelityFactor>0
%     warning('Occlusion gaps component not implemented for continuous case!');
% end

% if strcmpi(opt.conOpt.alg,'CGD')
%     opt.conOpt.maxIter=200;

end
