function [fx, dfx]=EbndBS(coefs,splines,opt)

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

splos=0; % splines segments offset
Delta=opt.conOpt.enParEseg(1);

sna=[splines(:).number]*2;
sna=cumsum(sna)-sna;

for id=1:N
    spl=splines(id);
    
    % how many segments?
    %     pieces=spl.pieces;
    %     breaks=spl.breaks;
    % breaks=[1 10];
    
    tr=spl.start:spl.end;
    index = spl.index;
    knots = spl.knots;
    sn=spl.number;
    
    os=splos;
    os=sna(id);
    %     locc=tr-breaks(index);
    
    lcnt=0;
    for t=1:length(tr)
        lcnt=lcnt+1;
        
        if t==1 || t==length(tr)
            ilcnt=index(lcnt);
            bslcnt=ilcnt+3;
            curt=tr(lcnt);
            curt_2=curt-Delta;
            if t==length(tr), curt_2=curt+Delta; end
            
            tx=knots(bslcnt-2:bslcnt+3);%-tr(lcnt);
            
            
            
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
            
            
            obj=((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)/t24 - (((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2x_cmt1_c1x_cmt4/t14 - (c3x*cmt2 - c2x*cmt5)/t25 - ((c1x - c2x)*cmt4)/t14 + ((c2x - c3x)*cmt2)/t25))/t24 - (cmt3*((c3x*cmt2 - c2x*cmt5)/t25 - (c4x*cmt3 - c3x*cmt6)/t36 - ((c2x - c3x)*cmt5)/t25 + ((c3x - c4x)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4x*(curt_2 - t3) - c3x*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))/t14 - (c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - ((c1x - c2x)*(curt_2 - t4))/t14 + ((c2x - c3x)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - (c4x*(curt_2 - t3) - c3x*(curt_2 - t6))/t36 - ((c2x - c3x)*(curt_2 - t5))/t25 + ((c3x - c4x)*(curt_2 - t3))/t36))/t35)/(t3 - t4))^2 + ((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)/t24 - (((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2y_cmt1_c1y_cmt4/t14 - (c3y*cmt2 - c2y*cmt5)/t25 - ((c1y - c2y)*cmt4)/t14 + ((c2y - c3y)*cmt2)/t25))/t24 - (cmt3*((c3y*cmt2 - c2y*cmt5)/t25 - (c4y*cmt3 - c3y*cmt6)/t36 - ((c2y - c3y)*cmt5)/t25 + ((c3y - c4y)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4y*(curt_2 - t3) - c3y*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))/t14 - (c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - ((c1y - c2y)*(curt_2 - t4))/t14 + ((c2y - c3y)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - (c4y*(curt_2 - t3) - c3y*(curt_2 - t6))/t36 - ((c2y - c3y)*(curt_2 - t5))/t25 + ((c3y - c4y)*(curt_2 - t3))/t36))/t35)/(t3 - t4))^2;
            fx=fx+obj;
            dfx(svpos(1)) = dfx(svpos(1)) + -2*((3*cmt4^2)/(t14*t24*(t3 - t4)) - (3*(curt_2 - t4)^2)/(t14*t24*(t3 - t4)))*((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)/t24 - (((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2x_cmt1_c1x_cmt4/t14 - (c3x*cmt2 - c2x*cmt5)/t25 - ((c1x - c2x)*cmt4)/t14 + ((c2x - c3x)*cmt2)/t25))/t24 - (cmt3*((c3x*cmt2 - c2x*cmt5)/t25 - (c4x*cmt3 - c3x*cmt6)/t36 - ((c2x - c3x)*cmt5)/t25 + ((c3x - c4x)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4x*(curt_2 - t3) - c3x*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))/t14 - (c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - ((c1x - c2x)*(curt_2 - t4))/t14 + ((c2x - c3x)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - (c4x*(curt_2 - t3) - c3x*(curt_2 - t6))/t36 - ((c2x - c3x)*(curt_2 - t5))/t25 + ((c3x - c4x)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
            dfx(svpos(2)) = dfx(svpos(2)) + 2*((((cmt1*cmt4)/t14 + (cmt2*cmt5)/t25)/t24 + (cmt4*(cmt1/t14 + cmt2/t25 + cmt4/t14 + cmt5/t25))/t24 + cmt5^2/(t25*t35) + (2*cmt3*cmt5)/(t25*t35))/(t3 - t4) - ((((curt_2 - t1)*(curt_2 - t4))/t14 + ((curt_2 - t2)*(curt_2 - t5))/t25)/t24 + ((curt_2 - t4)*((curt_2 - t1)/t14 + (curt_2 - t2)/t25 + (curt_2 - t4)/t14 + (curt_2 - t5)/t25))/t24 + (curt_2 - t5)^2/(t25*t35) + (2*(curt_2 - t3)*(curt_2 - t5))/(t25*t35))/(t3 - t4))*((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)/t24 - (((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2x_cmt1_c1x_cmt4/t14 - (c3x*cmt2 - c2x*cmt5)/t25 - ((c1x - c2x)*cmt4)/t14 + ((c2x - c3x)*cmt2)/t25))/t24 - (cmt3*((c3x*cmt2 - c2x*cmt5)/t25 - (c4x*cmt3 - c3x*cmt6)/t36 - ((c2x - c3x)*cmt5)/t25 + ((c3x - c4x)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4x*(curt_2 - t3) - c3x*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))/t14 - (c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - ((c1x - c2x)*(curt_2 - t4))/t14 + ((c2x - c3x)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - (c4x*(curt_2 - t3) - c3x*(curt_2 - t6))/t36 - ((c2x - c3x)*(curt_2 - t5))/t25 + ((c3x - c4x)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
            dfx(svpos(3)) = dfx(svpos(3)) + -2*((((cmt2*cmt5)/t25 + (cmt3*cmt6)/t36)/t35 + (cmt3*(cmt2/t25 + cmt3/t36 + cmt5/t25 + cmt6/t36))/t35 + cmt2^2/(t24*t25) + (2*cmt2*cmt4)/(t24*t25))/(t3 - t4) - ((((curt_2 - t2)*(curt_2 - t5))/t25 + ((curt_2 - t3)*(curt_2 - t6))/t36)/t35 + ((curt_2 - t3)*((curt_2 - t2)/t25 + (curt_2 - t3)/t36 + (curt_2 - t5)/t25 + (curt_2 - t6)/t36))/t35 + (curt_2 - t2)^2/(t24*t25) + (2*(curt_2 - t2)*(curt_2 - t4))/(t24*t25))/(t3 - t4))*((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)/t24 - (((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2x_cmt1_c1x_cmt4/t14 - (c3x*cmt2 - c2x*cmt5)/t25 - ((c1x - c2x)*cmt4)/t14 + ((c2x - c3x)*cmt2)/t25))/t24 - (cmt3*((c3x*cmt2 - c2x*cmt5)/t25 - (c4x*cmt3 - c3x*cmt6)/t36 - ((c2x - c3x)*cmt5)/t25 + ((c3x - c4x)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4x*(curt_2 - t3) - c3x*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))/t14 - (c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - ((c1x - c2x)*(curt_2 - t4))/t14 + ((c2x - c3x)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - (c4x*(curt_2 - t3) - c3x*(curt_2 - t6))/t36 - ((c2x - c3x)*(curt_2 - t5))/t25 + ((c3x - c4x)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
            dfx(svpos(4)) = dfx(svpos(4)) + 2*((3*cmt3^2)/((t3 - t4)*t35*t36) - (3*(curt_2 - t3)^2)/((t3 - t4)*t35*t36))*((((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)/t24 - (((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2x_cmt1_c1x_cmt4/t14 - (c3x*cmt2 - c2x*cmt5)/t25 - ((c1x - c2x)*cmt4)/t14 + ((c2x - c3x)*cmt2)/t25))/t24 - (cmt3*((c3x*cmt2 - c2x*cmt5)/t25 - (c4x*cmt3 - c3x*cmt6)/t36 - ((c2x - c3x)*cmt5)/t25 + ((c3x - c4x)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4x*(curt_2 - t3) - c3x*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2x*(curt_2 - t1) - c1x*(curt_2 - t4))/t14 - (c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - ((c1x - c2x)*(curt_2 - t4))/t14 + ((c2x - c3x)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3x*(curt_2 - t2) - c2x*(curt_2 - t5))/t25 - (c4x*(curt_2 - t3) - c3x*(curt_2 - t6))/t36 - ((c2x - c3x)*(curt_2 - t5))/t25 + ((c3x - c4x)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
            dfx(svpos(5)) = dfx(svpos(5)) + -2*((3*cmt4^2)/(t14*t24*(t3 - t4)) - (3*(curt_2 - t4)^2)/(t14*t24*(t3 - t4)))*((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)/t24 - (((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2y_cmt1_c1y_cmt4/t14 - (c3y*cmt2 - c2y*cmt5)/t25 - ((c1y - c2y)*cmt4)/t14 + ((c2y - c3y)*cmt2)/t25))/t24 - (cmt3*((c3y*cmt2 - c2y*cmt5)/t25 - (c4y*cmt3 - c3y*cmt6)/t36 - ((c2y - c3y)*cmt5)/t25 + ((c3y - c4y)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4y*(curt_2 - t3) - c3y*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))/t14 - (c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - ((c1y - c2y)*(curt_2 - t4))/t14 + ((c2y - c3y)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - (c4y*(curt_2 - t3) - c3y*(curt_2 - t6))/t36 - ((c2y - c3y)*(curt_2 - t5))/t25 + ((c3y - c4y)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
            dfx(svpos(6)) = dfx(svpos(6)) + 2*((((cmt1*cmt4)/t14 + (cmt2*cmt5)/t25)/t24 + (cmt4*(cmt1/t14 + cmt2/t25 + cmt4/t14 + cmt5/t25))/t24 + cmt5^2/(t25*t35) + (2*cmt3*cmt5)/(t25*t35))/(t3 - t4) - ((((curt_2 - t1)*(curt_2 - t4))/t14 + ((curt_2 - t2)*(curt_2 - t5))/t25)/t24 + ((curt_2 - t4)*((curt_2 - t1)/t14 + (curt_2 - t2)/t25 + (curt_2 - t4)/t14 + (curt_2 - t5)/t25))/t24 + (curt_2 - t5)^2/(t25*t35) + (2*(curt_2 - t3)*(curt_2 - t5))/(t25*t35))/(t3 - t4))*((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)/t24 - (((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2y_cmt1_c1y_cmt4/t14 - (c3y*cmt2 - c2y*cmt5)/t25 - ((c1y - c2y)*cmt4)/t14 + ((c2y - c3y)*cmt2)/t25))/t24 - (cmt3*((c3y*cmt2 - c2y*cmt5)/t25 - (c4y*cmt3 - c3y*cmt6)/t36 - ((c2y - c3y)*cmt5)/t25 + ((c3y - c4y)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4y*(curt_2 - t3) - c3y*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))/t14 - (c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - ((c1y - c2y)*(curt_2 - t4))/t14 + ((c2y - c3y)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - (c4y*(curt_2 - t3) - c3y*(curt_2 - t6))/t36 - ((c2y - c3y)*(curt_2 - t5))/t25 + ((c3y - c4y)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
            dfx(svpos(7)) = dfx(svpos(7)) + -2*((((cmt2*cmt5)/t25 + (cmt3*cmt6)/t36)/t35 + (cmt3*(cmt2/t25 + cmt3/t36 + cmt5/t25 + cmt6/t36))/t35 + cmt2^2/(t24*t25) + (2*cmt2*cmt4)/(t24*t25))/(t3 - t4) - ((((curt_2 - t2)*(curt_2 - t5))/t25 + ((curt_2 - t3)*(curt_2 - t6))/t36)/t35 + ((curt_2 - t3)*((curt_2 - t2)/t25 + (curt_2 - t3)/t36 + (curt_2 - t5)/t25 + (curt_2 - t6)/t36))/t35 + (curt_2 - t2)^2/(t24*t25) + (2*(curt_2 - t2)*(curt_2 - t4))/(t24*t25))/(t3 - t4))*((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)/t24 - (((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2y_cmt1_c1y_cmt4/t14 - (c3y*cmt2 - c2y*cmt5)/t25 - ((c1y - c2y)*cmt4)/t14 + ((c2y - c3y)*cmt2)/t25))/t24 - (cmt3*((c3y*cmt2 - c2y*cmt5)/t25 - (c4y*cmt3 - c3y*cmt6)/t36 - ((c2y - c3y)*cmt5)/t25 + ((c3y - c4y)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4y*(curt_2 - t3) - c3y*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))/t14 - (c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - ((c1y - c2y)*(curt_2 - t4))/t14 + ((c2y - c3y)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - (c4y*(curt_2 - t3) - c3y*(curt_2 - t6))/t36 - ((c2y - c3y)*(curt_2 - t5))/t25 + ((c3y - c4y)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
            dfx(svpos(8)) = dfx(svpos(8)) + 2*((3*cmt3^2)/((t3 - t4)*t35*t36) - (3*(curt_2 - t3)^2)/((t3 - t4)*t35*t36))*((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)/t24 - (((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)/t35 + (cmt4*(c2y_cmt1_c1y_cmt4/t14 - (c3y*cmt2 - c2y*cmt5)/t25 - ((c1y - c2y)*cmt4)/t14 + ((c2y - c3y)*cmt2)/t25))/t24 - (cmt3*((c3y*cmt2 - c2y*cmt5)/t25 - (c4y*cmt3 - c3y*cmt6)/t36 - ((c2y - c3y)*cmt5)/t25 + ((c3y - c4y)*cmt3)/t36))/t35)/(t3 - t4) - ((((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))*(curt_2 - t4))/t14 - ((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t2))/t25)/t24 - (((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))*(curt_2 - t5))/t25 - ((c4y*(curt_2 - t3) - c3y*(curt_2 - t6))*(curt_2 - t3))/t36)/t35 + ((curt_2 - t4)*((c2y*(curt_2 - t1) - c1y*(curt_2 - t4))/t14 - (c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - ((c1y - c2y)*(curt_2 - t4))/t14 + ((c2y - c3y)*(curt_2 - t2))/t25))/t24 - ((curt_2 - t3)*((c3y*(curt_2 - t2) - c2y*(curt_2 - t5))/t25 - (c4y*(curt_2 - t3) - c3y*(curt_2 - t6))/t36 - ((c2y - c3y)*(curt_2 - t5))/t25 + ((c3y - c4y)*(curt_2 - t3))/t36))/t35)/(t3 - t4));
        end
        
    end
    splos = splos + sn*2;
end

end