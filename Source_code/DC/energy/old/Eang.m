function [fx, dfx]=Eang(coefs,splines)

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

splos=0; % splines segments offset
global sceneInfo;
frR=sceneInfo.frameRate;

k=.1; % parameter for Pseudo-Huber
epsil=.1; % parameter for Charbonnier

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
        
        os = (index(lcnt)-1)*offsetjump+splos;
        ax=coefs(1+os);bx=coefs(2+os);cx=coefs(3+os);dx=coefs(4+os);
        ay=coefs(5+os);by=coefs(6+os);cy=coefs(7+os);dy=coefs(8+os);
        
        tsq=t*t;
        px_=(cx + 2*bx*t + 3*ax*t^2);
        py_=(cy + 2*by*t + 3*ay*t^2);
        px__=(2*bx + 6*ax*t);
        py__=(2*by + 6*ay*t);
        
        % simple curvature
        % obj=frR*((py__*px_ - px__*py_)^2/(px_^2 + py_^2)^2)^(1/2);
        % fx=fx+obj;
        % dfx(1 + os) = dfx(1 + os) + -(frR*((2*(6*t*py_ - 3*t^2*py__)*(py__*px_ - px__*py_))/(px_^2 + py_^2)^2 + (12*t^2*(py__*px_ - px__*py_)^2*px_)/(px_^2 + py_^2)^3))/(2*((py__*px_ - px__*py_)^2/(px_^2 + py_^2)^2)^(1/2));
        % dfx(2 + os) = dfx(2 + os) + -(frR*((2*(py__*px_ - px__*py_)*(2*cy - 2*t*py__ + 4*by*t + 6*ay*t^2))/(px_^2 + py_^2)^2 + (8*t*(py__*px_ - px__*py_)^2*px_)/(px_^2 + py_^2)^3))/(2*((py__*px_ - px__*py_)^2/(px_^2 + py_^2)^2)^(1/2));
        % dfx(3 + os) = dfx(3 + os) + -(frR*((2*(py__*px_ - px__*py_)^2*(2*cx + 4*bx*t + 6*ax*t^2))/(px_^2 + py_^2)^3 - (2*py__*(py__*px_ - px__*py_))/(px_^2 + py_^2)^2))/(2*((py__*px_ - px__*py_)^2/(px_^2 + py_^2)^2)^(1/2));
        % dfx(4 + os) = dfx(4 + os) + 0;
        % dfx(5 + os) = dfx(5 + os) + (frR*((2*(6*t*px_ - 3*t^2*px__)*(py__*px_ - px__*py_))/(px_^2 + py_^2)^2 - (12*t^2*(py__*px_ - px__*py_)^2*py_)/(px_^2 + py_^2)^3))/(2*((py__*px_ - px__*py_)^2/(px_^2 + py_^2)^2)^(1/2));
        % dfx(6 + os) = dfx(6 + os) + (frR*((2*(py__*px_ - px__*py_)*(2*cx - 2*t*px__ + 4*bx*t + 6*ax*t^2))/(px_^2 + py_^2)^2 - (8*t*(py__*px_ - px__*py_)^2*py_)/(px_^2 + py_^2)^3))/(2*((py__*px_ - px__*py_)^2/(px_^2 + py_^2)^2)^(1/2));
        % dfx(7 + os) = dfx(7 + os) + -(frR*((2*(py__*px_ - px__*py_)^2*(2*cy + 4*by*t + 6*ay*t^2))/(px_^2 + py_^2)^3 + (2*px__*(py__*px_ - px__*py_))/(px_^2 + py_^2)^2))/(2*((py__*px_ - px__*py_)^2/(px_^2 + py_^2)^2)^(1/2));
        % dfx(8 + os) = dfx(8 + os) + 0;
        
        % lorentzian
        % obj=-log(1/((frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2)/(px_^2 + py_^2)^2 + 1));
        % fx=fx+obj;
        % dfx(1 + os) = dfx(1 + os) + -((12*frR^2*k^4*tsq*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2*px_)/(px_^2 + py_^2)^3 + (2*frR^2*k^2*(6*t*py_ - 3*tsq*py__)*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)*(py__*px_ - px__*py_))/((px_^2 + py_^2)^2*((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2)))/((frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2)/(px_^2 + py_^2)^2 + 1);
        % dfx(2 + os) = dfx(2 + os) + -((8*frR^2*k^4*t*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2*px_)/(px_^2 + py_^2)^3 + (2*frR^2*k^2*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)*(py__*px_ - px__*py_)*(2*cy - 2*t*py__ + 4*by*t + 6*ay*tsq))/((px_^2 + py_^2)^2*((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2)))/((frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2)/(px_^2 + py_^2)^2 + 1);
        % dfx(3 + os) = dfx(3 + os) + -((2*frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2*(2*cx + 4*bx*t + 6*ax*tsq))/(px_^2 + py_^2)^3 - (2*frR^2*k^2*py__*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)*(py__*px_ - px__*py_))/((px_^2 + py_^2)^2*((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2)))/((frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2)/(px_^2 + py_^2)^2 + 1);
        % dfx(4 + os) = dfx(4 + os) + 0;
        % dfx(5 + os) = dfx(5 + os) + -((12*frR^2*k^4*tsq*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2*py_)/(px_^2 + py_^2)^3 - (2*frR^2*k^2*(6*t*px_ - 3*tsq*px__)*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)*(py__*px_ - px__*py_))/((px_^2 + py_^2)^2*((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2)))/((frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2)/(px_^2 + py_^2)^2 + 1);
        % dfx(6 + os) = dfx(6 + os) + -((8*frR^2*k^4*t*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2*py_)/(px_^2 + py_^2)^3 - (2*frR^2*k^2*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)*(py__*px_ - px__*py_)*(2*cx - 2*t*px__ + 4*bx*t + 6*ax*tsq))/((px_^2 + py_^2)^2*((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2)))/((frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2)/(px_^2 + py_^2)^2 + 1);
        % dfx(7 + os) = dfx(7 + os) + -((2*frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2*(2*cy + 4*by*t + 6*ay*tsq))/(px_^2 + py_^2)^3 + (2*frR^2*k^2*px__*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)*(py__*px_ - px__*py_))/((px_^2 + py_^2)^2*((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2)))/((frR^2*k^4*(((py__*px_ - px__*py_)^2/k^2 + 1)^(1/2) - 1)^2)/(px_^2 + py_^2)^2 + 1);
        % dfx(8 + os) = dfx(8 + os) + 0;
        % fprintf('%i %.2f, %f\n',id, t, obj);
        
        % Lorentzian Charbonnier
        obj=-log(1/((frR^2*(epsil + (py__*px_ - px__*py_)^2))/(epsil + px_^2 + py_^2)^2 + 1));
        fx=fx+obj;
        dfx(1 + os) = dfx(1 + os) + -((2*frR^2*(6*t*py_ - 3*tsq*py__)*(py__*px_ - px__*py_))/(epsil + px_^2 + py_^2)^2 + (12*frR^2*tsq*(epsil + (py__*px_ - px__*py_)^2)*px_)/(epsil + px_^2 + py_^2)^3)/((frR^2*(epsil + (py__*px_ - px__*py_)^2))/(epsil + px_^2 + py_^2)^2 + 1);
        dfx(2 + os) = dfx(2 + os) + -((2*frR^2*(py__*px_ - px__*py_)*(2*cy - 2*t*py__ + 4*by*t + 6*ay*tsq))/(epsil + px_^2 + py_^2)^2 + (8*frR^2*t*(epsil + (py__*px_ - px__*py_)^2)*px_)/(epsil + px_^2 + py_^2)^3)/((frR^2*(epsil + (py__*px_ - px__*py_)^2))/(epsil + px_^2 + py_^2)^2 + 1);
        dfx(3 + os) = dfx(3 + os) + -((2*frR^2*(epsil + (py__*px_ - px__*py_)^2)*(2*cx + 4*bx*t + 6*ax*tsq))/(epsil + px_^2 + py_^2)^3 - (2*frR^2*py__*(py__*px_ - px__*py_))/(epsil + px_^2 + py_^2)^2)/((frR^2*(epsil + (py__*px_ - px__*py_)^2))/(epsil + px_^2 + py_^2)^2 + 1);
        dfx(4 + os) = dfx(4 + os) + 0;
        dfx(5 + os) = dfx(5 + os) + ((2*frR^2*(6*t*px_ - 3*tsq*px__)*(py__*px_ - px__*py_))/(epsil + px_^2 + py_^2)^2 - (12*frR^2*tsq*(epsil + (py__*px_ - px__*py_)^2)*py_)/(epsil + px_^2 + py_^2)^3)/((frR^2*(epsil + (py__*px_ - px__*py_)^2))/(epsil + px_^2 + py_^2)^2 + 1);
        dfx(6 + os) = dfx(6 + os) + ((2*frR^2*(py__*px_ - px__*py_)*(2*cx - 2*t*px__ + 4*bx*t + 6*ax*tsq))/(epsil + px_^2 + py_^2)^2 - (8*frR^2*t*(epsil + (py__*px_ - px__*py_)^2)*py_)/(epsil + px_^2 + py_^2)^3)/((frR^2*(epsil + (py__*px_ - px__*py_)^2))/(epsil + px_^2 + py_^2)^2 + 1);
        dfx(7 + os) = dfx(7 + os) + -((2*frR^2*(epsil + (py__*px_ - px__*py_)^2)*(2*cy + 4*by*t + 6*ay*tsq))/(epsil + px_^2 + py_^2)^3 + (2*frR^2*px__*(py__*px_ - px__*py_))/(epsil + px_^2 + py_^2)^2)/((frR^2*(epsil + (py__*px_ - px__*py_)^2))/(epsil + px_^2 + py_^2)^2 + 1);
        dfx(8 + os) = dfx(8 + os) + 0;
        
        
    end
    splos = splos + pieces*8;
end

end