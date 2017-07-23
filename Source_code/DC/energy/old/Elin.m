function [fx, dfx]=Elin(coefs,splines,parElin)

% Energy component for linear velocity
% 


% How many splines?
N=length(splines);

% state vector length
C=length(coefs);

% initialize value and derivative with 0
fx=0;
dfx=zeros(C,1);

offsetjump=8;
nf=parElin(1); 
speed=parElin(2);
frR=parElin(3);

splos=0; % splines segments offset
% global sceneInfo;
% frR=sceneInfo.frameRate;

allcnt=0;
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

    for t=locc
        lcnt=lcnt+1;

        os = (index(lcnt)-1)*offsetjump + splos;
        ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
        ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);
        tsq=t*t;
        px_=(cx + 2*bx*t + 3*ax*t^2);
        py_=(cy + 2*by*t + 3*ay*t^2);
        

obj=(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)^2;
fx=fx+obj;
dfx(1 + os) = dfx(1 + os) + -(6*frR*tsq*(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)*px_)/(nf*(px_^2 + py_^2)^(1/2));
dfx(2 + os) = dfx(2 + os) + -(4*frR*t*(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)*px_)/(nf*(px_^2 + py_^2)^(1/2));
dfx(3 + os) = dfx(3 + os) + -(frR*(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)*(2*cx + 4*bx*t + 6*ax*tsq))/(nf*(px_^2 + py_^2)^(1/2));
dfx(4 + os) = dfx(4 + os) + 0;
dfx(5 + os) = dfx(5 + os) + -(6*frR*tsq*(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)*py_)/(nf*(px_^2 + py_^2)^(1/2));
dfx(6 + os) = dfx(6 + os) + -(4*frR*t*(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)*py_)/(nf*(px_^2 + py_^2)^(1/2));
dfx(7 + os) = dfx(7 + os) + -(frR*(speed - (frR*(px_^2 + py_^2)^(1/2))/nf)*(2*cy + 4*by*t + 6*ay*tsq))/(nf*(px_^2 + py_^2)^(1/2));
dfx(8 + os) = dfx(8 + os) + 0;

    end
    splos = splos + pieces*8; 
end

end