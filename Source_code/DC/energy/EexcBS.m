function [fx, dfx]=EexcBS(coefs,splines,parEexc)

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

%%%%%%%% !!!!!!!!!!!!!
nf=1000;
nf=parEexc(3);
targetSize=parEexc(2);

siga=parEexc(1);
% siga=0.0001;

sigb=(targetSize/2)*siga; %%%%
% global opt sceneInfo
% sigscale=opt.proxcostFactor;
% sigscale=1;

sn=[splines(:).number]*2;
sn=cumsum(sn)-sn;


splos1=0; % splines segments offset

for id1=1:N
    spl1=splines(id1);
        index1 = spl1.index;
        knots1 = spl1.knots;
        sn1 = spl1.number;
    splos2=0;

    for id2=1:N                
        if id1>=id2, continue; end
        spl2=splines(id2);
        
        if spl1.bbox(1)>spl2.bbox(2), continue; end
        if spl1.bbox(2)<spl2.bbox(1), continue; end
        if spl1.bbox(3)>spl2.bbox(4), continue; end
        if spl1.bbox(4)<spl2.bbox(3), continue; end
        
        
        
        s1=spl1.start; s2=spl2.start;
        e1=spl1.end; e2=spl2.end;
        
        s1=round(s1);s2=round(s2);e1=round(e1);e2=round(e2);
%         [s1 e1]
%         [s2 e2]
        
        % no need to proceed if no temporal overlap
        if s1>e2 || s2>e1, continue; end
        
        timeoverlap=max(s1,s2):min(e1,e2);   
        tr=timeoverlap;    
    
    
        index2 = spl2.index;
        knots2 = spl2.knots;
        sn2 = spl2.number;
    
    
        lcnt=0;
        for t=1:length(tr)
            lcnt=lcnt+1;
            
            
            talong1=tr(lcnt)-s1+1;
            talong2=tr(lcnt)-s2+1;
            
            
            ilcnt1=index1(talong1);            
            bslcnt1=ilcnt1+3;
            
            ilcnt2=index2(talong2);
            bslcnt2=ilcnt2+3;
%             [tr(lcnt) talong1 talong2 bslcnt1 bslcnt2]
            
            curt=tr(lcnt);
            tx=knots1(bslcnt1-2:bslcnt1+3);
            tx_2=knots2(bslcnt2-2:bslcnt2+3);
        
            os1=splos1;
            os2=splos2;
            os1=sn(id1);os2=sn(id2);
            
            svpos1=[os1+ilcnt1:os1+ilcnt1+3 os1+ilcnt1+sn1:os1+ilcnt1+3+sn1];
            svpos2=[os2+ilcnt2:os2+ilcnt2+3 os2+ilcnt2+sn2:os2+ilcnt2+3+sn2];
            
            clocx=coefs(svpos1(1:4))';            clocy=coefs(svpos1(5:8))';
            clocx_2=coefs(svpos2(1:4))';            clocy_2=coefs(svpos2(5:8))';
        
            t1=tx(1);t2=tx(2);t3=tx(3);t4=tx(4); t5=tx(5);t6=tx(6);
            c1x=clocx(1);c2x=clocx(2);c3x=clocx(3);c4x=clocx(4);
            c1y=clocy(1);c2y=clocy(2);c3y=clocy(3);c4y=clocy(4);

            cmt1=curt-t1;cmt2=curt-t2;cmt3=curt-t3;cmt4=curt-t4;cmt5=curt-t5;cmt6=curt-t6;
            t14=t1-t4;        t24=t2-t4;        t25=t2-t5;        t35=t3-t5;        t36=t3-t6;
            c2x_cmt1_c1x_cmt4=(c2x*cmt1 - c1x*cmt4);
            c2y_cmt1_c1y_cmt4=(c2y*cmt1 - c1y*cmt4);
            
            t1_2=tx_2(1);
            t2_2=tx_2(2);
            t3_2=tx_2(3);
            t4_2=tx_2(4);
            t5_2=tx_2(5);
            t6_2=tx_2(6);
            
            c1x_2=clocx_2(1);
            c2x_2=clocx_2(2);
            c3x_2=clocx_2(3);
            c4x_2=clocx_2(4);
            
            c1y_2=clocy_2(1);
            c2y_2=clocy_2(2);
            c3y_2=clocy_2(3);
            c4y_2=clocy_2(4);
           
            
            px=-(((((c2x*(curt - t1) - c1x*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3x*(curt - t2) - c2x*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4x*(curt - t3) - c3x*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
            py=-(((((c2y*(curt - t1) - c1y*(curt - t4))*(curt - t4))/(t1 - t4) - ((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t2))/(t2 - t5))*(curt - t4))/(t2 - t4) - ((((c3y*(curt - t2) - c2y*(curt - t5))*(curt - t5))/(t2 - t5) - ((c4y*(curt - t3) - c3y*(curt - t6))*(curt - t3))/(t3 - t6))*(curt - t3))/(t3 - t5))/(t3 - t4);
            p2x=-(((((c2x_2*(curt - t1_2) - c1x_2*(curt - t4_2))*(curt - t4_2))/(t1_2 - t4_2) - ((c3x_2*(curt - t2_2) - c2x_2*(curt - t5_2))*(curt - t2_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2) - ((((c3x_2*(curt - t2_2) - c2x_2*(curt - t5_2))*(curt - t5_2))/(t2_2 - t5_2) - ((c4x_2*(curt - t3_2) - c3x_2*(curt - t6_2))*(curt - t3_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2))/(t3_2 - t4_2);
            p2y=-(((((c2y_2*(curt - t1_2) - c1y_2*(curt - t4_2))*(curt - t4_2))/(t1_2 - t4_2) - ((c3y_2*(curt - t2_2) - c2y_2*(curt - t5_2))*(curt - t2_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2) - ((((c3y_2*(curt - t2_2) - c2y_2*(curt - t5_2))*(curt - t5_2))/(t2_2 - t5_2) - ((c4y_2*(curt - t3_2) - c3y_2*(curt - t6_2))*(curt - t3_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2))/(t3_2 - t4_2);

            
            px_p2x=px-p2x; px_p2xsq = px_p2x*px_p2x;
            py_p2y=py-p2y; py_p2ysq = py_p2y*py_p2y;
            dist = sqrt(px_p2xsq + py_p2ysq);
            if dist>targetSize
%                 continue;
            end
                      
            tmp1=((((c3x*cmt2 - c2x*cmt5)*cmt5)/t25 - ((c4x*cmt3 - c3x*cmt6)*cmt3)/t36)*cmt3);
            tmp2=(((c2x_cmt1_c1x_cmt4*cmt4)/t14 - ((c3x*cmt2 - c2x*cmt5)*cmt2)/t25)*cmt4)/t24;
            tmp3=((c2x_2*(curt - t1_2) - c1x_2*(curt - t4_2))*(curt - t4_2))/(t1_2 - t4_2);
            tmp4=((tmp3 - ((c3x_2*(curt - t2_2) - c2x_2*(curt - t5_2))*(curt - t2_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2);
            tmp10=((((c3x_2*(curt - t2_2) - c2x_2*(curt - t5_2))*(curt - t5_2))/(t2_2 - t5_2) - ((c4x_2*(curt - t3_2) - c3x_2*(curt - t6_2))*(curt - t3_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2);
            tmp5=((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2))^2;
            tmp6=((((c2y_cmt1_c1y_cmt4*cmt4)/t14 - ((c3y*cmt2 - c2y*cmt5)*cmt2)/t25)*cmt4)/t24 - ((((c3y*cmt2 - c2y*cmt5)*cmt5)/t25 - ((c4y*cmt3 - c3y*cmt6)*cmt3)/t36)*cmt3)/t35)/(t3 - t4);
            tmp7=((((c2y_2*(curt - t1_2) - c1y_2*(curt - t4_2))*(curt - t4_2))/(t1_2 - t4_2) - ((c3y_2*(curt - t2_2) - c2y_2*(curt - t5_2))*(curt - t2_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2);
            tmp8=((((c3y_2*(curt - t2_2) - c2y_2*(curt - t5_2))*(curt - t5_2))/(t2_2 - t5_2) - ((c4y_2*(curt - t3_2) - c3y_2*(curt - t6_2))*(curt - t3_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2);
            expins=sigb - siga*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2);
            expexp=exp(expins);
            expponesq=(expexp + 1)^2;
            
            tmp9=siga*expexp;
            
            
            tmp11=(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3 - t4));
            
            obj=1 - 1/(expexp + 1);
            fx=fx+obj;

            dfx(svpos1(1)) = dfx(svpos1(1)) + (tmp9*cmt4^3*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*t14*t24*(t3 - t4));
            dfx(svpos1(2)) = dfx(svpos1(2)) + -(tmp9*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2))*((((cmt1*cmt4)/t14 + (cmt2*cmt5)/t25)*cmt4)/t24 + (cmt3*cmt5^2)/(t25*t35)))/tmp11;
            dfx(svpos1(3)) = dfx(svpos1(3)) + (tmp9*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2))*((((cmt2*cmt5)/t25 + (cmt3*cmt6)/t36)*cmt3)/t35 + (cmt2^2*cmt4)/(t24*t25)))/tmp11;
            dfx(svpos1(4)) = dfx(svpos1(4)) + -(tmp9*cmt3^3*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3 - t4)*t35*t36);
            dfx(svpos1(5)) = dfx(svpos1(5)) + (tmp9*cmt4^3*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*t14*t24*(t3 - t4));
            dfx(svpos1(6)) = dfx(svpos1(6)) + -(tmp9*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))*((((cmt1*cmt4)/t14 + (cmt2*cmt5)/t25)*cmt4)/t24 + (cmt3*cmt5^2)/(t25*t35)))/tmp11;
            dfx(svpos1(7)) = dfx(svpos1(7)) + (tmp9*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))*((((cmt2*cmt5)/t25 + (cmt3*cmt6)/t36)*cmt3)/t35 + (cmt2^2*cmt4)/(t24*t25)))/tmp11;
            dfx(svpos1(8)) = dfx(svpos1(8)) + -(tmp9*cmt3^3*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3 - t4)*t35*t36);
            dfx(svpos2(1)) = dfx(svpos2(1)) + -(tmp9*(curt - t4_2)^3*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t1_2 - t4_2)*(t2_2 - t4_2)*(t3_2 - t4_2));
            dfx(svpos2(2)) = dfx(svpos2(2)) + (tmp9*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2))*(((((curt - t1_2)*(curt - t4_2))/(t1_2 - t4_2) + ((curt - t2_2)*(curt - t5_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2) + ((curt - t3_2)*(curt - t5_2)^2)/((t2_2 - t5_2)*(t3_2 - t5_2))))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3_2 - t4_2));
            dfx(svpos2(3)) = dfx(svpos2(3)) + -(tmp9*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2))*(((((curt - t2_2)*(curt - t5_2))/(t2_2 - t5_2) + ((curt - t3_2)*(curt - t6_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2) + ((curt - t2_2)^2*(curt - t4_2))/((t2_2 - t4_2)*(t2_2 - t5_2))))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3_2 - t4_2));
            dfx(svpos2(4)) = dfx(svpos2(4)) + (tmp9*(curt - t3_2)^3*((tmp2 - tmp1/t35)/(t3 - t4) - (tmp4 - tmp10)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3_2 - t4_2)*(t3_2 - t5_2)*(t3_2 - t6_2));
            dfx(svpos2(5)) = dfx(svpos2(5)) + -(tmp9*(curt - t4_2)^3*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t1_2 - t4_2)*(t2_2 - t4_2)*(t3_2 - t4_2));
            dfx(svpos2(6)) = dfx(svpos2(6)) + (tmp9*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))*(((((curt - t1_2)*(curt - t4_2))/(t1_2 - t4_2) + ((curt - t2_2)*(curt - t5_2))/(t2_2 - t5_2))*(curt - t4_2))/(t2_2 - t4_2) + ((curt - t3_2)*(curt - t5_2)^2)/((t2_2 - t5_2)*(t3_2 - t5_2))))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3_2 - t4_2));
            dfx(svpos2(7)) = dfx(svpos2(7)) + -(tmp9*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))*(((((curt - t2_2)*(curt - t5_2))/(t2_2 - t5_2) + ((curt - t3_2)*(curt - t6_2))/(t3_2 - t6_2))*(curt - t3_2))/(t3_2 - t5_2) + ((curt - t2_2)^2*(curt - t4_2))/((t2_2 - t4_2)*(t2_2 - t5_2))))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3_2 - t4_2));
            dfx(svpos2(8)) = dfx(svpos2(8)) + (tmp9*(curt - t3_2)^3*(tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2)))/(expponesq*(tmp5 + (tmp6 - (tmp7 - tmp8)/(t3_2 - t4_2))^2)^(1/2)*(t3_2 - t4_2)*(t3_2 - t5_2)*(t3_2 - t6_2));


        end
%         splos2 = splos2 + sn2*2;
%         splos1 = splos1 + sn1*2;
    end
    
end

end