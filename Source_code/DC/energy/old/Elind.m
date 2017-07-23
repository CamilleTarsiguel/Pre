function [fx, dfx]=Elind(x,spl)

% TODO
% Comment, docs



tr=spl.start:spl.end;
T=length(tr);


X=x(1:2:end);
Y=x(2:2:end);

% [X Y]

% initialize value and derivative with 0
fx=0;
dfx=zeros(length(x),1);

id=1;
xind=1;
yind=xind+1;

t=2;
a=X(t-1,id);        b=Y(t-1,id);
c=X(t,id);          d=Y(t,id);

%         dfx(xind) = dfx(xind) + (2*a - 2*c)/(2*((a - c)^2 + (b - d)^2)^(1/2));
%         dfx(yind) = dfx(yind) + (2*b - 2*d)/(2*((a - c)^2 + (b - d)^2)^(1/2));


    xind=xind+2;
    yind=xind+1;
    
for t=1:T-1
		% (a,b) = past frame position
		% (c,d) = current frame position
		% (e,f) = next frame position        
        a=X(t,id);        b=Y(t,id);
        c=X(t+1,id);          d=Y(t+1,id);
        e=X(t+1,id);        f=Y(t+1,id);
        obj=((a - c)^2 + (b - d)^2)^(1/2);
        fx = fx + obj;
    
        dfx(xind-2) = dfx(xind-2) + (2*a - 2*c)/(2*((a - c)^2 + (b - d)^2)^(1/2));
        dfx(yind-2) = dfx(yind-2) + (2*b - 2*d)/(2*((a - c)^2 + (b - d)^2)^(1/2));
        dfx(xind) = dfx(xind) -(2*a - 2*c)/(2*((a - c)^2 + (b - d)^2)^(1/2));
        dfx(yind) = dfx(yind) -(2*b - 2*d)/(2*((a - c)^2 + (b - d)^2)^(1/2));
        
        xind=xind+2;
        yind=xind+1;
end

% dfx(xind-2) = dfx(xind-2) + (2*a - 2*c)/(2*((a - c)^2 + (b - d)^2)^(1/2));
% dfx(yind-2) = dfx(yind-2) + (2*b - 2*d)/(2*((a - c)^2 + (b - d)^2)^(1/2));



end