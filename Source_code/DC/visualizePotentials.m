%% used to create Fig. 5 PAMI
x=linspace(-5,5,100);
k=1;
epsil=.1;
pseudohuber=k^2*(sqrt((x./k).^2+1)-1);
huber=k*(abs(x)-k/2); huber(abs(x)<k)=x(abs(x)<k).^2*(1/2);
charbonnier=k^2 * sqrt(epsil + (x./k).^2);
absval=abs(x);

clf; hold on
plot(x,absval,x,huber,x,pseudohuber,x,charbonnier);

legend('abs','hub','ps-hub','charb')