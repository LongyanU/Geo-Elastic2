% Please first run Garvin2(0.5,0.55) to get GarvinResult.mat
% and then run this program so as to compare with analytic results
clear;
clc
close all;

nt=3346/2;
isnap=20;    % snapshot sampling

dx=10;
h=10;
%nx=430;
nx=700;
%nz=220;
nz=400;

x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
% dt=0.003; % calculate time step from stability criterion
dt=0.001; % calculate time step from stability criterion
tau=dt;


%f0=13;
t=(1:nt)*dt;
%t0=2/f0;                       % initialize time axis

f0=15;  %f0=1/tp %the peak frequency of the seismic source
t0=1.5/f0;     %ts

%src=10^7*exp(-(pi*f0*(t-t0)).^2);              % source time function
%src=-diff((src))/dx^2;				% time derivative to obtain gaussian
src=10^2*((pi*f0*(t-t0)).^2-0.5).*exp(-(pi*f0*(t-t0)).^2);

for i=1:nz
    for j=1:nx
        %vp(i,j)=4000;
        vp(i,j)=1732.10;
    end
end


for i=1:nz
    for j=1:nx
        %vs(i,j)=2310;
        vs(i,j)=1000.0;
    end
end

for i=1:nz
    for j=1:nx
        %rho(i,j)=2.5*10^6;
        rho(i,j)=1.0*10^3;
    end
end
vp(1:46,:)=1000.0*2^0.5;
rho(1:46,:)=0.5*10^3;
p=zeros([nz nx]); Vx=p; Vz=p;
Txxx=p;
Txzz=p;
Tzzz=p;
Txzx=p;

Txx=p;
Txz=p;
Tzz=p;


Vxx=p;
Vzz=p;
Vxz=p;
Vzx=p;

% coeff=[ 1.53927, -0.25946, 0.066042, -0.0168935, 0.00363274, -0.000556345, 0.0000444248];  % space domain k Smaller

coeff=[ 1.58988, -0.303127, 0.0976618, -0.0355769, 0.0122981, -0.00340316, 0.000538606]% time-space domain

ixs=nx/2;
izs=146; %%%%%%


Seismic_Vx=zeros(nt,nx);
Seismic_Vz=zeros(nt,nx);
Seismic_Txx=zeros(nt,nx);
Seismic_Tzz=zeros(nt,nx);
Seismic_Txz=zeros(nt,nx);

tic
for it=1:nt-2,
    
    %Txx/x
    Txxx=coeff(1)*( (Txx)-circshift(Txx,[0,1]))+...
        coeff(2)*(circshift(Txx,[0,-1])-circshift(Txx,[0,2]))+...
        coeff(3)*( circshift(Txx,[0,-2])-circshift(Txx,[0,3]))+...
        coeff(4)*( circshift(Txx,[0,-3])-circshift(Txx,[0,4]))+...
        coeff(5)*( circshift(Txx,[0,-4])-circshift(Txx,[0,5]))+...
        coeff(6)*( circshift(Txx,[0,-5])-circshift(Txx,[0,6]))+...
        coeff(7)*( circshift(Txx,[0,-6])-circshift(Txx,[0,7]));
    
    %Txz/z
    %Txzz(i,j)=(Txz(i+1,j)-Txz(i,j));
    Txzz=(circshift(Txz,[ -1])-circshift(Txz,[ 0]));
    
    %Tzz/z
    Tzzz=coeff(1)*( (Tzz)-circshift(Tzz,[1]))+...
        coeff(2)*(circshift(Tzz,[-1])-circshift(Tzz,[2]))+...
        coeff(3)*( circshift(Tzz,[-2])-circshift(Tzz,[3]))+...
        coeff(4)*( circshift(Tzz,[-3])-circshift(Tzz,[4]))+...
        coeff(5)*( circshift(Tzz,[-4])-circshift(Tzz,[5]))+...
        coeff(6)*( circshift(Tzz,[-5])-circshift(Tzz,[6]))+...
        coeff(7)*( circshift(Tzz,[-6])-circshift(Tzz,[7]));
    
    %Txz/x
    %Txzx(i,j)=(Txz(i,j+1)-Txz(i,j));
    Txzx=(circshift(Txz,[0 -1])-Txz);
    
    
    
    Vx=Vx+1./(rho).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rho).*dt.*(Tzzz+Txzx)/h;
    
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    
    seismogramVx(it)=Vx(101,400);
    seismogramVz(it)=Vz(101,400);
    
    
    Seismic_Vx(it,:)=Vx(46,:);
    Seismic_Vz(it,:)=Vz(46,:);
    
    
    Vxz=coeff(1)*(circshift(Vx,[0])-circshift(Vx,[1]))+...
        coeff(2)*(circshift(Vx,[-1])-circshift(Vx,[2]))+...
        coeff(3)*(circshift(Vx,[-2])-circshift(Vx,[3]))+...
        coeff(4)*(circshift(Vx,[-3])-circshift(Vx,[4]))+...
        coeff(5)*(circshift(Vx,[-4])-circshift(Vx,[5]))+...
        coeff(6)*(circshift(Vx,[-5])-circshift(Vx,[6]))+...
        coeff(7)*(circshift(Vx,[-6])-circshift(Vx,[7]));
    
    Vzz=(circshift(Vz,[-1])-circshift(Vz,[0]));
    Vxx=(circshift(Vx,[0 -1])-circshift(Vx,[0 0]));
    Vzx=coeff(1)*(circshift(Vz,[0 0])-circshift(Vz,[0 1]))+...
        coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
        coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
        coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
        coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
        coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
        coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
    
    Txx=Txx+dt*rho.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rho.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rho.*(vs.^2).*(Vxz+Vzx)/h;
    
    Seismic_Txx(it,:)=Txx(46,:);
    Seismic_Tzz(it,:)=Tzz(46,:);
    Seismic_Txz(it,:)=Txz(46,:);
    
    
    Txx(izs,ixs)=Txx(izs,ixs)+src(it);
    Tzz(izs,ixs)=Tzz(izs,ixs)+src(it);
    
    
    Tzz(1:46,:)=0;
    Txz(1:46,:)=0;
    if rem(it,isnap)== 0,
        imagesc(-Txx(46:end,:),[-12 12]), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
end
toc

save('NonBalancedScheme1msChangePosition.mat')

% figure;imagesc(Vx(46:end,46:end))
% colormap gray
% hold on; plot(400,100,'r*','linewidth',2)
% hold on; plot(350,140,'bo','linewidth',2)
% grid on
% legend('The red is the receiver location(55,400)', 'the blue is the seismic location(145-45,350)');
% title('the grid size is 10 m')


figure;imagesc(vs(46:end,46:end))

% hold on; plot(350,145,'r*','linewidth',2)
% hold on; plot(400,100,'bo','linewidth',2)

hold on; plot(350,100,'r*','linewidth',2)
hold on; plot(400,55,'bo','linewidth',2)

grid on
legend('the red is the seismic location','The blue is the receiver location' );
xlabel('x/dx')
ylabel('z/dz')

% load('GarvinResult.mat') %obtain by  Garvin2(0.5,0.55)
load('GarvinResult.mat')
figure;plot(2.4*10^8/400*seismogramVx(101:end),'r')
temp=-conv(diff(src),u)/400;
hold on; plot(temp(101:end),'b')
axis([0 1500 -400/400 600/400 ])
legend('Non-Balanced FD method','Analytical Method')
xlabel('Travel time(ms)')
ylabel('Amp');

figure;plot(2.4*10^8/400*seismogramVz(101:end),'r')
temp=-conv(diff(src),w)/400;
hold on; plot(temp(101:end),'b')
axis([0 1500 -500/400 1 ])
legend('Non-Balanced FD method','Analytical Method')
xlabel('Travel time(ms)')
ylabel('Amp');
