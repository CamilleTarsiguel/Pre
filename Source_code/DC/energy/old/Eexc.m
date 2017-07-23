function [fx, dfx]=Eexc(coefs,splines)

% Energy component for mutual exclusion
%


% How many splines?
N=length(splines);

% state vector length
C=length(coefs);

% initialize value and derivative with 0
fx=0;
dfx=zeros(C,1);

offsetjump=8;

splos=0; % splines segments offset

Npieces=sum([splines(:).pieces]);

allbreaks=[splines(:).breaks];
broffset=[1 cumsum([splines(1:end-1).pieces]+1)+1];
% assert(Npieces==length(broffset),'breaks offset is wrong');

% which piece belongs to which spline?
% where do we find the breaks?
pcindex=[];
brindex=[]; brindxcnt=1;
for id=1:N
    pcindex=[pcindex id*ones(1,splines(id).pieces)];
    brindex=[brindex brindxcnt:brindxcnt+splines(id).pieces-1];
    brindxcnt=brindxcnt+splines(id).pieces+1;
end
% pcindex
% brindex

%
nf=1000;

siga=.05;
% siga=0.0001;
sigb=(350/2)*siga; %%%%
global opt
sigscale=opt.proxcostFactor;
sigscale=1;

global proxEexc proxPieces
proxEexc=zeros(N);
proxPieces=zeros(Npieces);
omitting=false(N,N,max([splines(:).end]));
for p1=1:Npieces
    id1=pcindex(p1);
    spl1=splines(id1);
    s1=allbreaks(brindex(p1)); e1=allbreaks(brindex(p1)+1);
    
    if spl1.start>s1, s1=spl1.start; end
    if spl1.end<e1, e1=spl1.end; end
    
    % only if we are in the first piece of spline
    if (p1>1 && pcindex(p1-1)~=id1) || p1==1
    if spl1.start<spl1.breaks(1)-.5
        s1=spl1.start;
    end
    end
    
    % only if we are in the last piece of spline
    if (p1<Npieces && pcindex(p1+1)~=id1) || p1==Npieces
    if spl1.end>spl1.breaks(end)+.5
        e1=spl1.end;
    end
    end
    s1=round(s1); e1=round(e1);
    
    
    
    for p2=p1+1:Npieces
        % what splines are we in?
        
        id2=pcindex(p2);
        
        
        spl2=splines(id2);
        
        
        s2=allbreaks(brindex(p2)); e2=allbreaks(brindex(p2)+1);
        %                 [s1 e1 s2 e2]
        % splines (breaks) sometimes start or end beyond defined points
        
        if spl2.start>s2, s2=spl2.start; end
        
        if spl2.end<e2, e2=spl2.end; end
        
            % only if we are in the first piece of spline
    if (p2>1 && pcindex(p2-1)~=id2) || p2==1
        if spl2.start<spl2.breaks(1)-.5
            s2=spl2.start;
        end
    end
            % only if we are in the last piece of spline
    if (p2<Npieces && pcindex(p2+1)~=id2) || p2==Npieces
        if spl2.end>spl2.breaks(end)+.5
            e2=spl2.end;
        end
    end
        
        %                 s1=spl1.start;
        
        %                 s1=max(s1,spl1.start); e1=min(e1,spl1.end);
        %                 s2=max(s2,spl2.start); e2=min(e2,spl2.end);
        % [s1 e1 s2 e2]
        
        s2=round(s2); e2=round(e2);
%                         [s1 e1 s2 e2]
        % 							if (p1==1 && p2==2)
        % fprintf('%i %i %i %i %f %f %f %f\n',p1,p2,id1,id2,s1,e1,s2,e2);
        % 			end
        %         s1=spl1.start; e1=spl1.end;
        %         s2=spl2.start; e2=spl2.end;
        
        % temporal overlap exists if
        % s1 starts before s2 ends ...and... s1 ends after s2 starts
        
        leaveout=find(omitting(id1,id2,:))';
        
        tmpol = (s1<=e2) && (e1>=s2);
        timeoverlap=max(s1,s2):min(e1,e2);        
        timeoverlap=setdiff(timeoverlap,leaveout);
        
%         leaveout=find(omitting(id2,id1,:))';
%         timeoverlap=setdiff(timeoverlap,leaveout);
        %         omit = length(timeoverlap)==1;
        
        
        %         if omit
        %             if spl1.start==spl2.end || spl2.start == spl1.end
        %                 omit=0;
        %             end
        %             (brindex(p1)==spl1.pieces || brindex(p2)==spl2.pieces )
        %             omit=0;
        %         end
        %         omit
%                 timeoverlap
%                 [s1 e1 s2 e2]
        %         brindex(p1)
        %         spl1.pieces
        %         brindex(p2)
        %         spl2.pieces
        %         omit
        %         fprintf('%i %i %i %i, %i. %i %i %i %i, %i %i %i %i\n',id1,id2,p1,p2,omit,s1,s2,e1,e2,spl1.start,spl2.start,spl1.end,spl2.end);
        if tmpol && id1~=id2
            %             timeoverlap=intersect(s1:e1,s2:e2);
            
%                         timeoverlap
            %             timeoverlap=round(timeoverlap);
            t1locc=timeoverlap-allbreaks(brindex(p1));
            t2locc=timeoverlap-allbreaks(brindex(p2));
%                         [id1 s1 e1]
            
%                         [id2 s2 e2]
            %             [p1 p2]
            %             t1locc
            %             t2locc
            %             pause
            
            os1=(p1-1)*8; os2=(p2-1)*8;
            
            ax=coefs(1+os1);bx=coefs(2+os1);cx=coefs(3+os1);dx=coefs(4+os1);
            ay=coefs(5+os1);by=coefs(6+os1);cy=coefs(7+os1);dy=coefs(8+os1);
            
            a2x=coefs(1+os2);b2x=coefs(2+os2);c2x=coefs(3+os2);d2x=coefs(4+os2);
            a2y=coefs(5+os2);b2y=coefs(6+os2);c2y=coefs(7+os2);d2y=coefs(8+os2);
            
            dcnt=1;
            tmpprox=0;
            for tl=1:length(timeoverlap)
                t1=t1locc(tl);t2=t2locc(tl);
                px = t1^3*ax + t1^2*bx + t1*cx + dx; py = t1^3*ay + t1^2*by + t1*cy + dy;
                p2x = t2^3*a2x + t2^2*b2x + t2*c2x + d2x; p2y = t2^3*a2y + t2^2*b2y + t2*c2y + d2y;
%                 p2x = t2*(t2*(t2*a2x+b2x)+c2x)+d2x;
                
                dist(dcnt)=sqrt((px-p2x)^2+(py-p2y)^2);
                dcnt=dcnt+1;
                %                 dist
                %                 if dist<700
                %                 [id1 id2 t1]
                %                 [px py p2x p2y]
                tmp1=(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2;
                tmp2=exp(sigb - siga*(tmp1 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                tmp3=(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2);
                
                dist=sqrt((px-p2x)^2+(py-p2y)^2);
                % obj=sigscale - sigscale/(tmp2 + 1);
                obj=1-1*1./(1+exp(-siga*dist+sigb));
                %                 [id1 id2 timeoverlap(tl)]
                %                 (py-p2y)^2
                %                 (px-p2x)^2
                %                 dist
                %                 obj
                tmpprox=tmpprox+obj;
                
                %                 tmpprox
                
                %  				if p1==1 && p2==9
                % 					fprintf('%.15f %.15f %.15f %.15f %.15f\n',sigb,-siga*dist,-siga*dist+sigb,exp(-siga*dist+sigb),obj);
%                 if timeoverlap(tl)==3
%                 px
%                 py
%                 p2x
%                 p2y
%                 end
                
%                                 fprintf('%i %i %i %f %f %.15f %.15f\n',id1, id2, timeoverlap(tl),t1, t2, dist,obj);
                %  				end
                
                % obj=1 - 1/(exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1);
                fx=fx+obj;
                % fx
                % tmpprox
                dfx(1 + os1) = dfx(1 + os1) + (siga*t1^3*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(2 + os1) = dfx(2 + os1) + (siga*t1^2*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(3 + os1) = dfx(3 + os1) + (siga*t1*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(4 + os1) = dfx(4 + os1) + (siga*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(2*d2x - 2*dx + 2*c2x*t2 - 2*cx*t1 + 2*a2x*t2^3 - 2*ax*t1^3 + 2*b2x*t2^2 - 2*bx*t1^2))/(2*(exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(5 + os1) = dfx(5 + os1) + (siga*t1^3*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(6 + os1) = dfx(6 + os1) + (siga*t1^2*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(7 + os1) = dfx(7 + os1) + (siga*t1*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(8 + os1) = dfx(8 + os1) + (siga*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(2*d2y - 2*dy + 2*c2y*t2 - 2*cy*t1 + 2*a2y*t2^3 - 2*ay*t1^3 + 2*b2y*t2^2 - 2*by*t1^2))/(2*(exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(1 + os2) = dfx(1 + os2) + -(siga*t2^3*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(2 + os2) = dfx(2 + os2) + -(siga*t2^2*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(3 + os2) = dfx(3 + os2) + -(siga*t2*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(4 + os2) = dfx(4 + os2) + -(siga*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(2*d2x - 2*dx + 2*c2x*t2 - 2*cx*t1 + 2*a2x*t2^3 - 2*ax*t1^3 + 2*b2x*t2^2 - 2*bx*t1^2))/(2*(exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(5 + os2) = dfx(5 + os2) + -(siga*t2^3*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(6 + os2) = dfx(6 + os2) + -(siga*t2^2*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(7 + os2) = dfx(7 + os2) + -(siga*t2*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2))/((exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                dfx(8 + os2) = dfx(8 + os2) + -(siga*exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2))*(2*d2y - 2*dy + 2*c2y*t2 - 2*cy*t1 + 2*a2y*t2^3 - 2*ay*t1^3 + 2*b2y*t2^2 - 2*by*t1^2))/(2*(exp(sigb - siga*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2)) + 1)^2*((d2x - dx + c2x*t2 - cx*t1 + a2x*t2^3 - ax*t1^3 + b2x*t2^2 - bx*t1^2)^2 + (d2y - dy + c2y*t2 - cy*t1 + a2y*t2^3 - ay*t1^3 + b2y*t2^2 - by*t1^2)^2)^(1/2));
                
                
            end
            proxPieces(p1,p2)=proxPieces(p1,p2)+tmpprox;
            proxEexc(id1,id2)=proxEexc(id1,id2)+tmpprox;
            %             dist
            %             sum(dist)
            %             pause
        
%         omitting(id2,id1,timeoverlap)=1;
            
        end
        
        %         [id1 id2]
        omitting(id1,id2,timeoverlap)=1;
    end
    
end
%
%
%
%     % how many segments?
%     pieces=spl.pieces;
%     breaks=spl.breaks;
%     % breaks=[1 10];
%
%     tr=spl.start:spl.end;
%     index = spl.index;
%
%     locc=tr-breaks(index);
%
%     lcnt=0;
%
%     for t=locc
%         lcnt=lcnt+1;
%
%         os = (index(lcnt)-1)*offsetjump + splos;
%         ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
%         ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);
%         a2x=coefs(1+os2);b2x=coefs(2+os2);c2x=coefs(3+os2);d2x=coefs(4+os2);
%         a2y=coefs(5+os2);b2y=coefs(6+os2);c2y=coefs(7+os2);d2y=coefs(8+os2);
%
%
%         obj=sigscale - sigscale/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1);
%         fx=fx+obj;
%         dfx(1 + os) = dfx(1 + os) + (2*siga*sigscale*t^3*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(2 + os) = dfx(2 + os) + (2*siga*sigscale*t^2*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(3 + os) = dfx(3 + os) + (2*siga*sigscale*t*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(4 + os) = dfx(4 + os) + (siga*sigscale*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(2*d2x - 2*dx + 2*c2x*t - 2*cx*t + 2*a2x*t^3 - 2*ax*t^3 + 2*b2x*t^2 - 2*bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(5 + os) = dfx(5 + os) + (2*siga*sigscale*t^3*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(6 + os) = dfx(6 + os) + (2*siga*sigscale*t^2*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(7 + os) = dfx(7 + os) + (2*siga*sigscale*t*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(8 + os) = dfx(8 + os) + (siga*sigscale*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(2*d2y - 2*dy + 2*c2y*t - 2*cy*t + 2*a2y*t^3 - 2*ay*t^3 + 2*b2y*t^2 - 2*by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(9 + os) = dfx(9 + os) + -(2*siga*sigscale*t^3*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(10 + os) = dfx(10 + os) + -(2*siga*sigscale*t^2*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(11 + os) = dfx(11 + os) + -(2*siga*sigscale*t*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(12 + os) = dfx(12 + os) + -(siga*sigscale*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(2*d2x - 2*dx + 2*c2x*t - 2*cx*t + 2*a2x*t^3 - 2*ax*t^3 + 2*b2x*t^2 - 2*bx*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(13 + os) = dfx(13 + os) + -(2*siga*sigscale*t^3*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(14 + os) = dfx(14 + os) + -(2*siga*sigscale*t^2*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(15 + os) = dfx(15 + os) + -(2*siga*sigscale*t*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%         dfx(16 + os) = dfx(16 + os) + -(siga*sigscale*exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2))*(2*d2y - 2*dy + 2*c2y*t - 2*cy*t + 2*a2y*t^3 - 2*ay*t^3 + 2*b2y*t^2 - 2*by*t^2))/(exp(sigb - siga*((d2x - dx + c2x*t - cx*t + a2x*t^3 - ax*t^3 + b2x*t^2 - bx*t^2)^2 + (d2y - dy + c2y*t - cy*t + a2y*t^3 - ay*t^3 + b2y*t^2 - by*t^2)^2)) + 1)^2;
%
%     end
%     splos = splos + pieces*8;
% end

end