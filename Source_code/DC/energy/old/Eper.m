function [fx, dfx]=Eper(coefs,splines, parEper)

% TODO
% Comment, docs

global opt sceneInfo T;
% How many splines?
N=length(splines);

% state vector length
C=length(coefs);

% initialize value and derivative with 0
fx=0;
dfx=zeros(C,1);

offsetjump=8;

splos=0; % splines segments offset

nf=parEper(1); % normalize to from mm to meter
k=parEper(2); % parameter for enter buffer
ksqnfsq=k^2*nf^2;

imOnGP=sceneInfo.imOnGP;
x1=imOnGP(1);y1=imOnGP(2);
x2=imOnGP(3);y2=imOnGP(4);
x3=imOnGP(5);y3=imOnGP(6);
x4=imOnGP(7);y4=imOnGP(8);

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
    
    if spl.start>1.5
        
        % persistence start
        t=locc(1);
        
        os = (index(1)-1)*offsetjump+splos;
        %         [id os]
        ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
        ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);
        px=dx + cx*t + ax*t^3 + bx*t^2; py=dy + cy*t + ay*t^3 + by*t^2;
        
        dl = ((y1-y2)*px + (x2-x1)*py + (x1*y2-x2*y1) ) / (sqrt((x2-x1)^2 + (y2-y1)^2) );
        dt = ((y2-y3)*px + (x3-x2)*py + (x2*y3-x3*y2) ) / (sqrt((x3-x2)^2 + (y3-y2)^2) );
        dr = ((y3-y4)*px + (x4-x3)*py + (x3*y4-x4*y3) ) / (sqrt((x4-x3)^2 + (y4-y3)^2) );
        db = ((y4-y1)*px + (x1-x4)*py + (x4*y1-x1*y4) ) / (sqrt((x1-x4)^2 + (y1-y4)^2) );
        
        % left
        
        x0=imOnGP(1);y0=imOnGP(2); x1=imOnGP(3);y1=imOnGP(4);
        dl = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); dl=dl/nf;
        ol = (k^2/6)*(1-(1-(dl/k)^2)^3);
        derxl=  (6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        deryl= -(6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        if abs(dl)>=k, derxl=0; deryl=0; ol=k^2/6; end
        ol=ol*6/k;
        
        
        % top
        x0=imOnGP(3);y0=imOnGP(4); x1=imOnGP(5);y1=imOnGP(6);
        dt = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); dt=dt/nf;
        ot = (k^2/6)*(1-(1-(dt/k)^2)^3);
        derxt=  (6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        deryt= -(6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        if abs(dt)>=k, derxt=0; deryt=0; ot=k^2/6;  end
        ot=ot*6/k;
        
        % right
        x0=imOnGP(5);y0=imOnGP(6); x1=imOnGP(7);y1=imOnGP(8);
        dr = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); dr=dr/nf;
        or = (k^2/6)*(1-(1-(dr/k)^2)^3);
        derxr=  (6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        deryr= -(6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        if abs(dr)>=k, derxr=0; deryr=0; or=k^2/6;  end
        or=or*6/k;
        
        % bottom
        x0=imOnGP(7);y0=imOnGP(8); x1=imOnGP(1);y1=imOnGP(2);
        db = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); db=db/nf;
        ob = (k^2/6)*(1-(1-(db/k)^2)^3);
        derxb=  (6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(y0 - y1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        deryb= -(6*((x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1))^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*(x0 - x1)*(x0*y1 - x1*y0 - dy*(x0 - x1) + dx*(y0 - y1)))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        if abs(db)>=k, derxb=0; deryb=0; ob=k^2/6;  end
        ob=ob*6/k;
        
%         [dl dr dt db]
        % [ol ot or ob]
        obj=ol*ot*or*ob;
        derx=derxl*ot*or*ob + ol*derxt*or*ob + ol*ot*derxr*ob + ol*ot*or*derxb;
        dery=deryl*ot*or*ob + ol*deryt*or*ob + ol*ot*deryr*ob + ol*ot*or*deryb;
        
        fx=fx+obj;
        dfx(4 + os) = derx;
        dfx(8 + os) = dery;
%         disp([derx derxl derxt derxr derxb])
        
    end
    

    derl=zeros(8,1);dert=zeros(8,1);derr=zeros(8,1);derb=zeros(8,1);
    % persistence end
    if spl.end<T-.5
%         [id os]
        t=locc(end);
        tsq=t*t; tcu=tsq*t;
        
        os = (index(length(locc))-1)*offsetjump+splos;
        %         [id os]
        ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
        ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);
        px = t^3*ax + t^2*bx + t*cx + dx; py = t^3*ay + t^2*by + t*cy + dy;
        
        % left
        x0=imOnGP(1);y0=imOnGP(2); x1=imOnGP(3);y1=imOnGP(4);
        dl = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); dl=dl/nf;
        ol = (k^2/6)*(1-(1-(dl/k)^2)^3);
        derl(1) = -(6*tcu*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derl(2) = -(6*tsq*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derl(3) = -(6*t*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derl(4) = -(6*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derl(5) = (6*tcu*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derl(6) = (6*tsq*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derl(7) = (6*t*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derl(8) = (6*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));       
        if abs(dl)>=k, derl(:)=0; ol=k^2/6; end
        ol=ol*6/k;

        % top
        x0=imOnGP(3);y0=imOnGP(4); x1=imOnGP(5);y1=imOnGP(6);
        dt = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); dt=dt/nf;
        ot = (k^2/6)*(1-(1-(dt/k)^2)^3);
        dert(1) = -(6*tcu*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        dert(2) = -(6*tsq*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        dert(3) = -(6*t*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        dert(4) = -(6*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        dert(5) = (6*tcu*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        dert(6) = (6*tsq*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        dert(7) = (6*t*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        dert(8) = (6*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));       
        if abs(dt)>=k, dert(:)=0; ot=k^2/6; end
        ot=ot*6/k;
        
        % right
        x0=imOnGP(5);y0=imOnGP(6); x1=imOnGP(7);y1=imOnGP(8);
        dr = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); dr=dr/nf;
        or = (k^2/6)*(1-(1-(dr/k)^2)^3);
        derr(1) = -(6*tcu*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derr(2) = -(6*tsq*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derr(3) = -(6*t*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derr(4) = -(6*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derr(5) = (6*tcu*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derr(6) = (6*tsq*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derr(7) = (6*t*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derr(8) = (6*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));       
        if abs(dr)>=k, derr(:)=0; or=k^2/6; end
        or=or*6/k;
        
        % bottom
        x0=imOnGP(7);y0=imOnGP(8); x1=imOnGP(1);y1=imOnGP(2);
        db = ((y0-y1)*px + (x1-x0)*py + (x0*y1-x1*y0) ) / (sqrt((x1-x0)^2 + (y1-y0)^2) ); db=db/nf;
        ob = (k^2/6)*(1-(1-(db/k)^2)^3);
        derb(1) = -(6*tcu*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derb(2) = -(6*tsq*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derb(3) = -(6*t*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derb(4) = -(6*(y0 - y1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derb(5) = (6*tcu*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derb(6) = (6*tsq*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derb(7) = (6*t*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));
        derb(8) = (6*(x0 - x1)*(((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0)^2/(ksqnfsq*((x0 - x1)^2 + (y0 - y1)^2)) - 1)^2*((x0 - x1)*py - (y0 - y1)*px - x0*y1 + x1*y0))/(k*nf^2*((x0 - x1)^2 + (y0 - y1)^2));       
        if abs(db)>=k, derb(:)=0; ob=k^2/6; end
        ob=ob*6/k;
        
        
        obj=ol*ot*or*ob;
%         obj=ol*ot;
        
        allder=derl*ot*or*ob + ol*dert*or*ob + ol*ot*derr*ob + ol*ot*or*derb;
%         dery=deryl*ot*or*ob + ol*deryt*or*ob + ol*ot*deryr*ob + ol*ot*or*deryb;
%         derx=derl;
%         allder=derl*ot + ol*dert;


        fx=fx+obj;
        
        dfx(os+1:os+8) = allder;
        
        
        
    end
    
    splos = splos + pieces*8;
end

end