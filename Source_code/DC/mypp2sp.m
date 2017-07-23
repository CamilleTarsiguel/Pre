function sp=mypp2sp(pp)

% newpp=splinefit(pp.start:pp.end,ppval(pp,pp.start:pp.end),pp.pieces,pp.order);

% newpp=pp;
    fbr=pp.breaks(1); lbr=pp.breaks(end);
    spknots=[fbr fbr fbr pp.breaks lbr lbr lbr];
%     sp=spap2(spknots,sporder,xd,yd);
    
sp=spap2(spknots,pp.order,pp.start:pp.end,ppval(pp,pp.start:pp.end));


% sp=pp2sp(newpp);

sp.start=pp.start;
sp.end=pp.end;
sp.index=pp.index;
sp.ptindex=pp.ptindex;

end