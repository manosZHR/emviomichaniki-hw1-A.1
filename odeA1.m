function [dcdt]=odeA1(t,c,par,nnr)

Cinf=par(1); k=par(2); %rn_origin=par(3); rn_last=par(4);

rn=linspace(0,1,nnr);
drn = (1-0)/nnr;

dcdt=zeros(nnr,1); 

dcdt(1)=0;
dcdt(nnr)=0;

for i=2:nnr-1
    
    dcdt(i) = k/drn^2 * ( (drn+rn(i))/rn(i)*c(i+1) - 2*c(i) + (rn(i)-drn)/rn(i)*c(i-1) );

end
