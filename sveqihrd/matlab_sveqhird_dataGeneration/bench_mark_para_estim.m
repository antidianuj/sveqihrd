close all 
clear all 
clc

rtrue = [1,1,1,1,0,1,1,1,1];
y0 = [1 1 0 1 0 0,0,0];
tspan = linspace(0,5);
soltrue = ode45(@(t,y)sveqhird(t,y,rtrue),tspan,y0);
yvalstrue = deval(soltrue,tspan);
for i = 1:8
    subplot(4,2,i)
    plot(tspan,yvalstrue(i,:))
    title(['y(',num2str(i),')'])
end

r = optimvar('r',9,"LowerBound",-10,"UpperBound",10);

myfcn = fcn2optimexpr(@RtoODE,r,tspan,y0);
obj = sum(sum((myfcn - yvalstrue).^2));
prob = optimproblem("Objective",obj);
show(prob)

r0.r = [0 0 0 0 0 0 0 0 0];
[rsol,sumsq] = solve(prob,r0);
disp("Estimated:");
disp(rsol.r)
disp("True");
disp(rtrue)