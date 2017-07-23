function [fx, dfx]=EdatBS(coefs,splines,parEdat, alldpoints, newLabeling)

% TODO
% Comment, docs

% How many splines?
N=length(splines);

% state vector length
C=length(coefs);

% initialize value and derivative with 0
fx=0;
dfx=zeros(C,1);
offsetjump=1;


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
    
    % how many segments?
%     pieces=spl.pieces;
%     breaks=spl.breaks;
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
    index=spl.ptindex;
    knots = spl.knots;
    sn=spl.number;
%     locc=tr-breaks(index);
    
    %     index
    %     locc
    %     pause
    %     tr
    
        
    lcnt=0;
    for det=1:npts
        
        lcnt=lcnt+1;
        ilcnt=index(lcnt);
        bslcnt=ilcnt+3;
        curt=tr(lcnt);
        
        tx=knots(bslcnt-2:bslcnt+3);%-tr(lcnt);
   
        
        os=splos;
        svpos=[os+ilcnt:os+ilcnt+3 os+ilcnt+sn:os+ilcnt+3+sn];
        clocx=coefs(svpos(1:4))';
        clocy=coefs(svpos(5:8))';
        
        t1=tx(1);t2=tx(2);t3=tx(3);t4=tx(4); t5=tx(5);t6=tx(6);
        c1x=clocx(1);c2x=clocx(2);c3x=clocx(3);c4x=clocx(4);
        c1y=clocy(1);c2y=clocy(2);c3y=clocy(3);c4y=clocy(4);
        
        cmt1=curt-t1;cmt2=curt-t2;cmt3=curt-t3;cmt4=curt-t4;cmt5=curt-t5;cmt6=curt-t6;
        t14=t1-t4;        t24=t2-t4;        t25=t2-t5;        t35=t3-t5;        t36=t3-t6;
        c2x_cmt1_c1x_cmt4=(c2x*cmt1 - c1x*cmt4);
        c2y_cmt1_c1y_cmt4=(c2y*cmt1 - c1y*cmt4);
        

        Dx=alldpts.xp(det);Dy=alldpts.yp(det); sc=alldpts.sp(det);
        
        % Euclidian distance squared        
        if dataFunction==2
            obj=(sc*((dx - Dx + cx*t + ax*t^3 + bx*t^2)^2 + (dy - Dy + cy*t + ay*t^3 + by*t^2)^2))/nf^2;
            
            fx=fx+obj;
            fxcnt=fxcnt+1;
            error('data L2 not implemented yet');
            
            % obj
            % pause

            
        elseif (dataFunction==1 || dataFunction==4)

            % simple Charbonnier
            obj=sc*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2);
            fx=fx+obj;
            dfx(svpos(1)) = dfx(svpos(1)) + -(sc*(Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*cmt4^3)/(nf^2*t14*t24*(t3 - t4)*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));
            dfx(svpos(2)) = dfx(svpos(2)) + (sc*(Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*((((cmt1*cmt4)/t14 + (cmt2*cmt5)/t25)*cmt4)/t24 + (cmt3*cmt5^2)/(t25*t35)))/(nf^2*(t3 - t4)*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));
            dfx(svpos(3)) = dfx(svpos(3)) + -(sc*(Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*((((cmt2*cmt5)/t25 + (cmt3*cmt6)/t36)*cmt3)/t35 + (cmt2^2*cmt4)/(t24*t25)))/(nf^2*(t3 - t4)*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));
            dfx(svpos(4)) = dfx(svpos(4)) + (sc*(Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*cmt3^3)/(nf^2*(t3 - t4)*t35*t36*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));
            dfx(svpos(5)) = dfx(svpos(5)) + -(sc*(Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*cmt4^3)/(nf^2*t14*t24*(t3 - t4)*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));
            dfx(svpos(6)) = dfx(svpos(6)) + (sc*(Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*((((cmt1*cmt4)/t14 + (cmt2*cmt5)/t25)*cmt4)/t24 + (cmt3*cmt5^2)/(t25*t35)))/(nf^2*(t3 - t4)*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));
            dfx(svpos(7)) = dfx(svpos(7)) + -(sc*(Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*((((cmt2*cmt5)/t25 + (cmt3*cmt6)/t36)*cmt3)/t35 + (cmt2^2*cmt4)/(t24*t25)))/(nf^2*(t3 - t4)*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));
            dfx(svpos(8)) = dfx(svpos(8)) + (sc*(Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))*cmt3^3)/(nf^2*(t3 - t4)*t35*t36*(k + ((Dx + ((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2 + (Dy + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4))^2)/nf^2)^(1/2));

        end
        
    end
    
    splos = splos + sn*2;
end

end