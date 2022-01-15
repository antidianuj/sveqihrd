clear all
close all 
clc

num_data=100;
for lambdy=1:num_data
OF=zeros(1,100);
para=zeros(length(OF),9);
in=10;
a=0.1;
b=0.1;
InitConds = [(b-a)*rand + a;(b-a)*rand + a;(b-a)*rand + a;(b-a)*rand + a;(b-a)*rand + a;(b-a)*rand + a;(b-a)*rand + a;(b-a)*rand + a];  

for mm=1:length(OF)
    a=0;
    b=2;
    P= (b-a)*rand + a;
    nu=(b-a)*rand + a;
    beta= (b-a)*rand + a;
    omegaq= (b-a)*rand + a;
    etas=(b-a)*rand+a;
    alphae= (b-a)*rand + a;
    deltai= (b-a)*rand + a;
    gamma= (b-a)*rand + a;
    epsilonh= (b-a)*rand + a;
    para(mm,:)=[P,nu,beta,omegaq,etas,alphae,deltai,gamma,epsilonh];

    dt=0.05; 
    T=10; 
    t=0:dt:T;

    SVEQIHRD = @(t,x) ([ P-(x(4)-omegaq)*x(4)-etas*x(1)-beta*(x(1)-etas)*x(1)
                        etas*x(1)-nu*x(2)
                        beta*(x(1)-etas)*x(1)-x(3)^2
                        nu*x(2)+(x(3)-alphae)*x(3)-x(4)^2
                        alphae*x(3)+omegaq*x(4)-x(5)^2
                        deltai*x(5)-x(6)^2
                        (x(5)-deltai-gamma*(x(5)-deltai))*x(5)+(x(6)-epsilonh)*x(6)
                        (gamma*(x(5)-deltai))*x(5)+epsilonh*x(6)]);
    




    [t,y] = ode45(SVEQIHRD, t, InitConds);
    OF(mm)=sum((y(:,1)-abs(y(:,1))).^2)+sum((y(:,2)-abs(y(:,2))).^2)+sum((y(:,3)-abs(y(:,3))).^2)+sum((y(:,4)-abs(y(:,4))).^2)+sum((y(:,5)-abs(y(:,5))).^2)+sum((y(:,6)-abs(y(:,6))).^2)+sum((y(:,7)-abs(y(:,7))).^2)+sum((y(:,8)-abs(y(:,8))).^2);
end
[CC,K]=min(OF);
vect_par=para(K,:);

P=vect_par(1);
nu=vect_par(2);
beta=vect_par(3);
omegaq=vect_par(4);
etas=vect_par(5);
alphae=vect_par(6);
deltai=vect_par(7);
gamma=vect_par(8);
epsilonh=vect_par(9);
dt=0.05; 
T=10; 
t=0:dt:T;

SVEQIHRD = @(t,x) ([ P-(x(4)-omegaq)*x(4)-etas*x(1)-beta*(x(1)-etas)*x(1)
                    etas*x(1)-nu*x(2)
                    beta*(x(1)-etas)*x(1)-x(3)^2
                    nu*x(2)+(x(3)-alphae)*x(3)-x(4)^2
                    alphae*x(3)+omegaq*x(4)-x(5)^2
                    deltai*x(5)-x(6)^2
                    (x(5)-deltai-gamma*(x(5)-deltai))*x(5)+(x(6)-epsilonh)*x(6)
                    (gamma*(x(5)-deltai))*x(5)+epsilonh*x(6)]);
                        
                        

[t,y] = ode45(SVEQIHRD, t, InitConds);
matty_boy=[t,y];
writematrix(matty_boy,join(['Input/dat',num2str(lambdy),'.csv']));
writematrix(vect_par,join(['Labels/lab',num2str(lambdy),'.csv']));
end

