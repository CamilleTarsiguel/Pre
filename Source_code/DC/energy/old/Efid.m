function [fx, dfx]=Efid(coefs,splines,parEfid,alldpoints)

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


global oldAlldpoints newLabeling;
% alldpoints=oldAlldpoints;

splos=0; % splines segments offset

tau=parEfid(1);


exthead=2;
exttail=2;

        siga=parEfid(1);
        sigb=parEfid(2)*siga;
        k=parEfid(3);
        ksq=k*k;

% T=5; %%%% !!!! %%%%
fxcnt=0;
for id=1:N
    spl=splines(id);
    
    % how many segments?
    pieces=spl.pieces;
    breaks=spl.breaks;
    % breaks=[1 10];
    
    tr=spl.start:spl.end;
    index = spl.index;
    locc=tr-breaks(index);
    
    lcnt=0;
    
%     fprintf('\nspline %i\n',id);
    
    for t=locc
        lcnt=lcnt+1;
        os = (index(lcnt)-1)*offsetjump + splos;
        ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
        ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);

        px=dx + cx*t + ax*t^3 + bx*t^2;
        py=dy + cy*t + ay*t^3 + by*t^2;
        tsq=t^2; tcu=tsq*t;
        
        % find all detections in current frame
        tPts=find(alldpoints.tp==tr(lcnt));
        
        alldpts.xp=alldpoints.xp(tPts);
        alldpts.yp=alldpoints.yp(tPts);
%         alldpts.sp=alldpoints.sp(tPts);
%         alldpts.tp=alldpoints.tp(tPts);
        npts=length(alldpts.xp);
        
%             fprintf('%i %i,',tr(lcnt),npts);
        
        tmpfx=ones(1,npts); tmpdfx=zeros(8,npts);
        for det=1:npts
		
%             lcnt=lcnt+1;
                        

            Dx=alldpts.xp(det);Dy=alldpts.yp(det);
            pxDx=(px-Dx)^2;            
                pyDy=(py-Dy)^2;
            dist=sqrt(pxDx+pyDy);
            disthu=sqrt((pxDx + pyDy)/ksq + 1);
%             if id==2 && tr(lcnt)==5 && det==11
%                 [t tr(lcnt) px py Dx Dy dist]
%             end
								
%             [tr(lcnt) px py Dx Dy dist]
            tmpexp=exp(sigb - siga * dist);
            tmpexp=exp(sigb - k * siga*(disthu - 1));
            obj=1/(tmpexp + 1);
%             obj=1/(tmpexp + 1);
%             if (id==1)
%                 fprintf('%i %f %f %f, %i %f %f %f\n',id,t,px,py,det,Dx,Dy,obj);
%             end
					%             obj
            tmpfx(det) = obj;
            tmp1=(tmpexp + 1)^2;
%             tmpfx
%             pause

tmpdfx(1,det) = (k*siga*tcu*tmpexp*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*disthu*tmp1);
tmpdfx(2,det) = (k*siga*tsq*tmpexp*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*disthu*tmp1);
tmpdfx(3,det) = (k*siga*t*tmpexp*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*disthu*tmp1);
tmpdfx(4,det) = (k*siga*tmpexp*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*ksq*disthu*tmp1);
tmpdfx(5,det) = (k*siga*tcu*tmpexp*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*disthu*tmp1);
tmpdfx(6,det) = (k*siga*tsq*tmpexp*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*disthu*tmp1);
tmpdfx(7,det) = (k*siga*t*tmpexp*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*disthu*tmp1);
tmpdfx(8,det) = (k*siga*tmpexp*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*ksq*disthu*tmp1);

% Charbonnier
% obj=1/(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1);
% fx=fx+obj;
% dfx(1 + os) = dfx(1 + os) + (k^2*siga*tcu*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));
% dfx(2 + os) = dfx(2 + os) + (k^2*siga*tsq*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));
% dfx(3 + os) = dfx(3 + os) + (k^2*siga*t*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(dx - Dx + cx*t + ax*tcu + bx*tsq))/(ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));
% dfx(4 + os) = dfx(4 + os) + (k^2*siga*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(2*dx - 2*Dx + 2*cx*t + 2*ax*tcu + 2*bx*tsq))/(2*ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));
% dfx(5 + os) = dfx(5 + os) + (k^2*siga*tcu*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));
% dfx(6 + os) = dfx(6 + os) + (k^2*siga*tsq*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));
% dfx(7 + os) = dfx(7 + os) + (k^2*siga*t*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(dy - Dy + cy*t + ay*tcu + by*tsq))/(ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));
% dfx(8 + os) = dfx(8 + os) + (k^2*siga*exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2))*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)*(2*dy - 2*Dy + 2*cy*t + 2*ay*tcu + 2*by*tsq))/(2*ksq*(exp(sigb - siga*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2)) + 1)^2*(((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2)*(k + k^2*((((dx - Dx + cx*t + ax*tcu + bx*tsq)^2 + (dy - Dy + cy*t + ay*tcu + by*tsq)^2)/ksq + 1)^(1/2) - 1)^2)^(1/2));


%             tmpdfx
%             pause
            
%             if any(isnan(tmpdfx(:)))
%                 dist
%                 obj
%                 tmp1
%                 (tmp1*dist)
%                 pause
%             end
            
            
        end
%         if id==2
%             [tr(lcnt) prod(tmpfx)]
%         end
%         if id==2 && tr(lcnt)==5
%             tmpfx
%             [tr(lcnt) prod(tmpfx)]
%         end
        fx=fx+prod(tmpfx);
        
        
    
        for co=1:8
            
% % % %             tmpfx
% % %             ttmpfx=repmat(tmpfx,npts,1);
% % % %             ttmpfx
% % %             inds=1:npts+1:npts*npts;
% % % %             inds
% % %             ttmpfx(inds)=tmpdfx(co,:);
% % % %             ttmpfx
% % %             ttmpfx=prod(ttmpfx);
% % %             fac=sum(ttmpfx);
% % % %             fac
            
            fac=0;
            for det=1:npts
                indobj=[1:det-1 det+1:npts];
%                 indobj=ones(1,npts);
%                 indobj(det)=0;
%                 indder=det;
%                 det
%                 tmpdfx(co,det)
%                 indobj.*(tmpfx)
%                 prod(indobj.*(tmpfx))
%                 pause
                fac=fac + tmpdfx(co,det)*prod(tmpfx(indobj));
            end
% % %             if isnan(fac)
% % %                 fac=0;
% % %                 for det=1:npts
% % %                     det
% % %                     indobj=[1:det-1 det+1:npts];
% % %                     
% % %                     Dx=alldpts.xp(det);Dy=alldpts.yp(det);
% % %                     pxDx=(px-Dx)^2;   pyDy=(py-Dy)^2;
% % %                     dist=sqrt(pxDx+pyDy);
% % % 
% % %                     tmpexp=exp(sigb - siga*dist);
% % %                     obj=1/(tmpexp + 1);
% % % 
% % %                     tmpfx(det) = obj;
% % %                     tmp1=(tmpexp + 1)^2; 
% % % %                     dist
% % % %                     obj
% % % %                     tmp1
% % % %             
% % %                     tmpdfx
% % %                     prod(tmpfx(indobj))
% % %                     pause
% % %                     fac=fac + tmpdfx(co,det)*prod(tmpfx(indobj));
% % %                 end
% % %             end

            dfx(co+os)=dfx(co+os)+fac;
        end
%         fx
%         pause
    end
%     [id fx]
    splos = splos + pieces*8;
end

end