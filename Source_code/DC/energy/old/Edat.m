function [fx, dfx]=Edat(coefs,splines,parEdat, alldpoints, newLabeling)

% TODO
% Comment, docs

% How many splines?
N=length(splines);

% state vector length
C=length(coefs);

% initialize value and derivative with 0
fx=0;
dfx=zeros(C,1);
offsetjump=8;


% global alldpoints newLabeling;


splos=0; % splines segments offset

k=parEdat(3); ksq=k*k;
nf=parEdat(2);
dataFunction=parEdat(1);

exthead=2;
exttail=2;

% T=5; %%%% !!!! %%%%
fxcnt=0;
for id=1:N
    %     id
    spl=splines(id);
    
    
    splstart=spl.start-exttail;
    splend=spl.end+exthead;
    %     splineTimespan=max(1,splstart):min(T,splend);
    
    
    %
    idPts=find(newLabeling==id);
    
    alldpts.xp=alldpoints.xp(idPts);
    alldpts.yp=alldpoints.yp(idPts);
    alldpts.sp=alldpoints.sp(idPts);
    alldpts.tp=alldpoints.tp(idPts);
    npts=length(alldpts.xp);
    
    %     inSplineTimespan=find(alldpts.tp>=splstart & alldpts.tp<=splend);
    %     tr=alldpts.tp(inSplineTimespan);
    %     notinSplineTimespan=true(1,nPoints);
    %     notinSplineTimespan(inSplineTimespan)=0;
    %     allothers= allpind & notinSplineTimespan;
    
    %     idPts
    %     alldpts
    % how many segments?
    pieces=spl.pieces;
    breaks=spl.breaks;
    % breaks=[1 10];
    
    %     tr=spl.start:spl.end;
    %
    %     at=tr;
    %     [~,index] = histc(at,spl.breaks);
    %     index
    %     lastedge=(spl.breaks(end)==at);
    %     index(lastedge)=pieces;
    %     index(index>pieces)=0;
    % %
    %     locc=-1*ones(1,length(at));
    %     locc(~~index)=at(~~index)-breaks(index(~~index));
    %
    %     tr=splstart:splend;
    tr= alldpts.tp;
    %[~,index] = histc(tr,[-inf,spl.breaks(2:spl.pieces),inf]);
    index=spl.ptindex;
    
    locc=tr-breaks(index);
    
    %     index
    %     locc
    %     pause
    %     tr
    
    lcnt=0;
    
    for det=1:npts
        lcnt=lcnt+1;
        t=locc(det);
        %         if t>=0
        tt=tr(lcnt);
        
        os = (index(lcnt)-1)*offsetjump + splos;
        ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
        ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);
        
        %             inThisFrame=find(alldpts.tp==tt);
        %             [id tt]
        %                      [t        lcnt]
        %                      inThisFrame
        %                      pause
        
        %             for det=inThisFrame
        Dx=alldpts.xp(det);Dy=alldpts.yp(det); sc=alldpts.sp(det);
        px=dx + cx*t + ax*t^3 + bx*t^2;
        tsq=t^2; tcu=tsq*t;
        % Euclidian distance squared
        %                 deltaxsq=(dx - Dx + cx*t + ax*t^3 + bx*t^2)^2;
        %                 deltaysq=(dy - Dy + cy*t + ay*t^3 + by*t^2)^2;
        
        if dataFunction==2
            obj=(sc*((dx - Dx + cx*t + ax*t^3 + bx*t^2)^2 + (dy - Dy + cy*t + ay*t^3 + by*t^2)^2))/nf^2;
            
            fx=fx+obj;
            fxcnt=fxcnt+1;
            
            % obj
            % pause
            dfx(1 + os) = dfx(1 + os) + (2*sc*t^3*(dx - Dx + cx*t + ax*t^3 + bx*t^2))/nf^2;
            dfx(2 + os) = dfx(2 + os) + (2*sc*t^2*(dx - Dx + cx*t + ax*t^3 + bx*t^2))/nf^2;
            dfx(3 + os) = dfx(3 + os) + (2*sc*t*(dx - Dx + cx*t + ax*t^3 + bx*t^2))/nf^2;
            dfx(4 + os) = dfx(4 + os) + (sc*(2*dx - 2*Dx + 2*cx*t + 2*ax*t^3 + 2*bx*t^2))/nf^2;
            dfx(5 + os) = dfx(5 + os) + (2*sc*t^3*(dy - Dy + cy*t + ay*t^3 + by*t^2))/nf^2;
            dfx(6 + os) = dfx(6 + os) + (2*sc*t^2*(dy - Dy + cy*t + ay*t^3 + by*t^2))/nf^2;
            dfx(7 + os) = dfx(7 + os) + (2*sc*t*(dy - Dy + cy*t + ay*t^3 + by*t^2))/nf^2;
            dfx(8 + os) = dfx(8 + os) + (sc*(2*dy - 2*Dy + 2*cy*t + 2*ay*t^3 + 2*by*t^2))/nf^2;
            
        elseif (dataFunction==1 || dataFunction==4)
            % Pseudo Huber
%             obj=k*sc*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2) - 1);
%             fx=fx+obj;
%             dfx(1 + os) = dfx(1 + os) + (k*sc*tcu*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             dfx(2 + os) = dfx(2 + os) + (k*sc*tsq*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             dfx(3 + os) = dfx(3 + os) + (k*sc*t*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             dfx(4 + os) = dfx(4 + os) + (k*sc*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             dfx(5 + os) = dfx(5 + os) + (k*sc*tcu*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             dfx(6 + os) = dfx(6 + os) + (k*sc*tsq*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             dfx(7 + os) = dfx(7 + os) + (k*sc*t*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             dfx(8 + os) = dfx(8 + os) + (k*sc*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*ksq*nf^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/(ksq*nf^2) + 1)^(1/2));
%             
            % simple Charbonnier
            obj=sc*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2);
fx=fx+obj;
dfx(1 + os) = dfx(1 + os) + (sc*tcu*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));
dfx(2 + os) = dfx(2 + os) + (sc*tsq*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));
dfx(3 + os) = dfx(3 + os) + (sc*t*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));
dfx(4 + os) = dfx(4 + os) + (sc*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));
dfx(5 + os) = dfx(5 + os) + (sc*tcu*(dy - Dy + cy*t + ay*tcu + by*tsq))/(nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));
dfx(6 + os) = dfx(6 + os) + (sc*tsq*(dy - Dy + cy*t + ay*tcu + by*tsq))/(nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));
dfx(7 + os) = dfx(7 + os) + (sc*t*(dy - Dy + cy*t + ay*tcu + by*tsq))/(nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));
dfx(8 + os) = dfx(8 + os) + (sc*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*nf^2*(k + ((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/nf^2)^(1/2));

        end
        
    end
    
    splos = splos + pieces*8;
end

end