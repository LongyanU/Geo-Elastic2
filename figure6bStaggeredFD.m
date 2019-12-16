% first run Figure11aStaggeredFDTra.m,figure11bStaggeredFD.m,figure11cKspace to get seismic records,
% then run figure11compareSeimogramsVz.m, Figure11CompareSeis_recordVz.m; figure12compareSeimogramsTxx.m,figure12CompareSeis_recordTxx
% to get the figures in figure 11 and figure 12.
% This is only for the convenient of the reviewers.
% 时间已过 647.480234 秒。

% 时间已过 433.010355 秒。
clear
clc %%%%%%%
close all
nt=4002;    % number of time steps
eps=.6;     % stability
isnap=20;    % snapshot sampling
load('vv')

c1=flipud(c);

v=c1;
nx=800;
nx=nx+45*2;
nz=475;
nz=nz+45*2;

vv=zeros(nz,nx);
for ii=1:nz-90
    for jj=1:nx-90
        vv(ii+45,jj+45)=v(ii,jj);
    end
end

for ii=1:nz-90  %%left
    for jj=1:45
        vv(ii+45,jj)=v(ii,1);
    end
end

for ii=1:nz-90  %%right
    for jj=nx-45:nx
        vv(ii+45,jj)=v(ii,800);
    end
end


for ii=1:45  %%top
    for jj=1:nx
        vv(ii,jj)=vv(46,jj);
    end
end

for ii=nz-44:nz  %%bottom
    for jj=1:nx
        vv(ii,jj)=vv(nz-45,jj);
    end
end




clear v
v=vv;
vp=v;
vs=vp/sqrt(3);
rou=ones(nz,nx);

dx=15;  %calculate space increment
h=dx;
x=(0:(nx-1))*dx;
z=(0:(nz-1))*dx;  % initialize space coordinates
dt=0.001; % calculate time step from stability criterion
tau=dt;


f0=23;
t=(1:nt)*dt;
t0=4/f0;                       % initialize time axis
src=10^2*exp(-f0^2*(t-t0).*(t-t0));              % source time function
src=-diff((src))/dx^2;				% time derivative to obtain gaussian


xs=floor(nx/2);
zs=46;

seis_recordVx=zeros(nt,nx);
seis_recordVz=zeros(nt,nx);
seis_recordTxx=zeros(nt,nx);
seis_recordTzz=zeros(nt,nx);
seis_recordTxz=zeros(nt,nx);

%  ******************************************************************
%  txz-------------vx
%  |                |
%  |                |
%  |                |
%  |                |
%  |                |
%  vz---------------txx,tzz
% Virieux 1986 Geophysics

dx=h;
dz=h;
tic

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


coeff=[1.57665 -0.290813 0.0879337 -0.0293731 0.00924336 -0.00239015 0.000409115];
for it=1:nt-2,
    
    %Txx/x
    Txxx=coeff(1)*(Txx-circshift(Txx,[0 1]))+...
        coeff(2)*(circshift(Txx,[0 -1])-circshift(Txx,[0 2]))+...
        coeff(3)*(circshift(Txx,[0 -2])-circshift(Txx,[0 3]))+...
        coeff(4)*(circshift(Txx,[0 -3])-circshift(Txx,[0 4]))+...
        coeff(5)*(circshift(Txx,[0 -4])-circshift(Txx,[0 5]))+...
        coeff(6)*(circshift(Txx,[0 -5])-circshift(Txx,[0 6]))+...
        coeff(7)*(circshift(Txx,[0 -6])-circshift(Txx,[0 7]));
    
    %Txz/z
    Txzz=(circshift(Txz,[ 1])-circshift(Txz,[ 0]));
       
    
    %Tzz/z
    Tzzz=coeff(1)*(circshift(Tzz,[ 0])-circshift(Tzz,[ -1]))+...
        coeff(2)*(circshift(Tzz,[ 1])-circshift(Tzz,[ -2]))+...
        coeff(3)*(circshift(Tzz,[ 2])-circshift(Tzz,[ -3]))+...
        coeff(4)*(circshift(Tzz,[ 3])-circshift(Tzz,[ -4]))+...
        coeff(5)*(circshift(Tzz,[ 4])-circshift(Tzz,[ -5]))+...
        coeff(6)*(circshift(Tzz,[ 5])-circshift(Tzz,[ -6]))+...
        coeff(7)*(circshift(Tzz,[ 6])-circshift(Tzz,[ -7]));
    
    
    %Txz/x
    Txzx=(circshift(Txz,[0 -1])-Txz);
      
    
    
    Vx=Vx+1./(rou).*dt.*(Txxx+Txzz)/h;
    Vz=Vz+1./(rou).*dt.*(Tzzz+Txzx)/h;
 
    [Vx,Vz]=spongeABC(Vx,Vz,nx,nz,45,45,0.009);
    seis_recordVx(it,:)=Vx(zs,:);
    seis_recordVz(it,:)=Vz(zs,:);
    
   Vxx=(circshift(Vx,[0 -1])-Vx);
     
    
    Vxz=coeff(1)*(circshift(Vx,[ 0])-circshift(Vx,[ -1]))+...
        coeff(2)*(circshift(Vx,[ 1])-circshift(Vx,[ -2]))+...
        coeff(3)*(circshift(Vx,[ 2])-circshift(Vx,[ -3]))+...
        coeff(4)*(circshift(Vx,[ 3])-circshift(Vx,[ -4]))+...
        coeff(5)*(circshift(Vx,[ 4])-circshift(Vx,[ -5]))+...
        coeff(6)*(circshift(Vx,[ 5])-circshift(Vx,[ -6]))+...
        coeff(7)*(circshift(Vx,[ 6])-circshift(Vx,[ -7]));
    
      Vzz=(circshift(Vz,[ 1])-circshift(Vz,[ 0]));
       
    
    Vzx=coeff(1)*(Vz-circshift(Vz,[0 1]))+...
        coeff(2)*(circshift(Vz,[0 -1])-circshift(Vz,[0 2]))+...
        coeff(3)*(circshift(Vz,[0 -2])-circshift(Vz,[0 3]))+...
        coeff(4)*(circshift(Vz,[0 -3])-circshift(Vz,[0 4]))+...
        coeff(5)*(circshift(Vz,[0 -4])-circshift(Vz,[0 5]))+...
        coeff(6)*(circshift(Vz,[0 -5])-circshift(Vz,[0 6]))+...
        coeff(7)*(circshift(Vz,[0 -6])-circshift(Vz,[0 7]));
    
    
    Txx=Txx+dt*rou.*(vp.^2.*Vxx+(vp.^2-2*vs.^2).*Vzz)/h;
    Tzz=Tzz+dt*rou.*(vp.^2.*Vzz+(vp.^2-2*vs.^2).*Vxx)/h;
    Txz=Txz+dt*rou.*(vs.^2).*(Vxz+Vzx)/h;
    
       Txx(zs,xs)=Txx(zs,xs)+src(it);
    Tzz(zs,xs)=Tzz(zs,xs)+src(it);
    
    seis_recordTxx(it,:)=Txx(zs,:);
    seis_recordTzz(it,:)=Tzz(zs,:);
    seis_recordTxz(it,:)=Txz(zs,:);
    
    if rem(it,isnap)== 0,
        imagesc(x,z,Vx), axis equal
        colormap gray
        xlabel('x'),ylabel('z')
        title(sprintf(' Time step: %i - Max ampl: %g ',it,max(max(Vx))))
        drawnow
    end
end

toc
save('StaggeredGridFD4.mat')