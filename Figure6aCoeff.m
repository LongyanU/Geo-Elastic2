clear;
clc;
close all
global M v dt h ratio

options = optimset('TolFun',10^-20,'TolX',10^-20,'MaxFunEvals',8000,'MaxIter',200);

for ii=1:3
    if ii==1
        M=3;
        ratio=0.4;
    elseif ii==2
        M=5;
        ratio=0.585;
    else
        M=7;
        ratio=0.69;
    end
    x0=0.01*ones(1,M);
    lb=-5*ones(M,1);
    ub=5*ones(M,1);
    [x,fval,out,iteration]= fmincon(@myfun3,x0,[],[],[],[],lb,ub,[],options) ;   % Invoke optimizer
    digits(6)
    vpa(x)
    
    
    
    k=linspace(1/10000,pi/h,100);
    temp1=0;
    for m=1:M
        temp1=temp1+2*x(m)*sin((m-1/2)*k*h);
    end
    
    a1=temp1.^2-k.^2*h^2;
    
    if ii==1
        figure;plot(a1,'k','LineWidth',2)
    elseif ii==2
        hold on;plot(a1,'m--','LineWidth',2);
    elseif ii==3
        hold on;plot(a1,'r:','LineWidth',2);
    end
end
xlabel('percentage of kh')
ylabel('E')
grid on
axis([0 100 -6*10^-4  2*10^-4])
legend('M=3','M=5','M=7')