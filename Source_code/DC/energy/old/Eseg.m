function [fx, dfx]=Eseg(coefs,splines, opt)

% To keep the splines consistent across
% consecutive pieces, a quadratic penalty is
% applied to coordinates and 1st derivative distances


% How many splines?
N=length(splines);

% state vector length
C=length(coefs);

% initialize value and derivative with 0
fx=0;
dfx=zeros(C,1);
offsetjump=8;

splos=0; % splines segments offset
% extraWt=0.0;
% extraWt=opt.conOpt.enParEseg(1);
extraWt=1;
Delta=opt.conOpt.enParEseg(1);

for id=1:N
    spl=splines(id);
    
    % how many segments?
    pieces=spl.pieces;
    % breaks=spl.breaks;

% % %     Ts=diff(spl.breaks);
% % % 
% % % 
% % %     % x,y coordinates
% % %     for seg=1:pieces-1
% % %         T=Ts(seg);
% % %         
% % %         os = (seg-1)*offsetjump + splos;
% % %         ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
% % %         ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);
% % %         d2x=coefs(12+os); d2y=coefs(16+os);
% % %         pxdx=(dx - d2x + T*cx + T^3*ax + T^2*bx);
% % %         pydy=(dy - d2y + T*cy + T^3*ay + T^2*by);
% % %         
% % %         obj=pxdx^2 + pydy^2;
% % %         fx=fx+obj;
% % %         
% % %         dfx(1+os) = dfx(1+os) + 2*T^3*pxdx;
% % %         dfx(2+os) = dfx(2+os) + 2*T^2*pxdx;
% % %         dfx(3+os) = dfx(3+os) + 2*T*pxdx;
% % %         dfx(4+os) = dfx(4+os) + 2*dx - 2*d2x + 2*T*cx + 2*T^3*ax + 2*T^2*bx;
% % %         dfx(5+os) = dfx(5+os) + 2*T^3*pydy;
% % %         dfx(6+os) = dfx(6+os) + 2*T^2*pydy;
% % %         dfx(7+os) = dfx(7+os) + 2*T*pydy;
% % %         dfx(8+os) = dfx(8+os) + 2*dy - 2*d2y + 2*T*cy + 2*T^3*ay + 2*T^2*by;
% % %         
% % %         % also take care of derivative for next piece
% % %         dfx(12+os) = dfx(12+os) + 2*d2x - 2*dx - 2*T*cx - 2*T^3*ax - 2*T^2*bx;
% % %         dfx(16+os) = dfx(16+os) + 2*d2y - 2*dy - 2*T*cy - 2*T^3*ay - 2*T^2*by;
% % %         
% % %     end
% % %     
% % %     % first derivative
% % %     for seg=1:pieces-1
% % %         T=Ts(seg);
% % %         
% % %         os = (seg-1)*offsetjump+splos;
% % %         ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);
% % %         ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);
% % %         c2x=coefs(11+os); c2y=coefs(15+os);
% % %         
% % %         pxdx=(cx - c2x + 2*T*bx + 3*T^2*ax);
% % %         pydy=(cy - c2y + 2*T*by + 3*T^2*ay);
% % %         
% % %         obj=pxdx^2 + pydy^2;
% % %         fx=fx+obj;
% % %         dfx(1 + os) = dfx(1 + os) + 6*T^2*pxdx;
% % %         dfx(2 + os) = dfx(2 + os) + 4*T*pxdx;
% % %         dfx(3 + os) = dfx(3 + os) + 2*cx - 2*c2x + 4*T*bx + 6*T^2*ax;
% % %         dfx(5 + os) = dfx(5 + os) + 6*T^2*pydy;
% % %         dfx(6 + os) = dfx(6 + os) + 4*T*pydy;
% % %         dfx(7 + os) = dfx(7 + os) + 2*cy - 2*c2y + 4*T*by + 6*T^2*ay;
% % %         dfx(11 + os) = dfx(11 + os) + 2*c2x - 2*cx - 4*T*bx - 6*T^2*ax;
% % %         dfx(15 + os) = dfx(15 + os) + 2*c2y - 2*cy - 4*T*by - 6*T^2*ay;
% % %     end
    
    %%%%%%%
    seg=1;
    tr=spl.start:spl.end;
    index = spl.index;
    breaks=spl.breaks;
    locc=tr-breaks(index);
    T=locc(1);
    T2=T-Delta;
        
    os = (seg-1)*offsetjump+splos;
    ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);
    ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);
    
    px_=cx + 2*T*bx + 3*T^2*ax;
    py_=cy + 2*T*by + 3*T^2*ay;

    px2_=cx + 2*T2*bx + 3*T2^2*ax;
    py2_=cy + 2*T2*by + 3*T2^2*ay;

    pxdx=(px_-px2_)^2;
    pydy=(py_-py2_)^2;
    obj=pxdx+pydy;
    
    fx=fx+extraWt*obj;
    dfx(1 + os) = dfx(1 + os) + extraWt*2*(3*T^2 - 3*T2^2)*(2*T*bx - 2*T2*bx + 3*T^2*ax - 3*T2^2*ax);
    dfx(2 + os) = dfx(2 + os) + extraWt*2*(2*T - 2*T2)*(2*T*bx - 2*T2*bx + 3*T^2*ax - 3*T2^2*ax);
    dfx(5 + os) = dfx(5 + os) + extraWt*2*(3*T^2 - 3*T2^2)*(2*T*by - 2*T2*by + 3*T^2*ay - 3*T2^2*ay);
    dfx(6 + os) = dfx(6 + os) + extraWt*2*(2*T - 2*T2)*(2*T*by - 2*T2*by + 3*T^2*ay - 3*T2^2*ay);

    %%%
    seg=pieces;
    T=locc(end);
    T2=T+Delta;

    os = (seg-1)*offsetjump+splos;
    ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);
    ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);
    
    px_=cx + 2*T*bx + 3*T^2*ax;
    py_=cy + 2*T*by + 3*T^2*ay;

    px2_=cx + 2*T2*bx + 3*T2^2*ax;
    py2_=cy + 2*T2*by + 3*T2^2*ay;

    pxdx=(px_-px2_)^2;
    pydy=(py_-py2_)^2;
    obj=pxdx+pydy;
    
    fx=fx+extraWt*obj;
    dfx(1 + os) = dfx(1 + os) + extraWt*2*(3*T^2 - 3*T2^2)*(2*T*bx - 2*T2*bx + 3*T^2*ax - 3*T2^2*ax);
    dfx(2 + os) = dfx(2 + os) + extraWt*2*(2*T - 2*T2)*(2*T*bx - 2*T2*bx + 3*T^2*ax - 3*T2^2*ax);
    dfx(5 + os) = dfx(5 + os) + extraWt*2*(3*T^2 - 3*T2^2)*(2*T*by - 2*T2*by + 3*T^2*ay - 3*T2^2*ay);
    dfx(6 + os) = dfx(6 + os) + extraWt*2*(2*T - 2*T2)*(2*T*by - 2*T2*by + 3*T^2*ay - 3*T2^2*ay);
    splos = splos + pieces*8;
end

end

