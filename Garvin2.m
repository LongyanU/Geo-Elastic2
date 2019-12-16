function Garvin2 (x,z,pois,h,Cs,rho)
% Solves the generalized Garvin problem of a line blast
% source applied at depth h within a homogeneous half-space.
% The response is sought within that same half-space
% at a receiver at range x and depth z.
%
% Written by Eduardo Kausel, MIT, Room 1-271, Cambridge, MA
% Version 1, July 19, 2011
%
%  Input arguments:
%     x    = range of receiver >0
%     z    = depth of receiver >=0
%     pois = Poisson's ratio        Defaults to 0.25 if not given
%     h    = Depth of source   > 0      "     "  1   "   "    "
%     Cs   = Shear wave velocity        "     "  1   "   "    "
%     rho  = mass density               "     "  1   "   "    "
%
% Sign convention:
%   x from left to right, z=0 at the surface, z points down
%   Displacements are positive down and to the right.
%
% References:
%   W.W. Garvin, Exact transient solution of the buried line source problem,
%                Proceedings of the Royal Society of London, Series A
%                Vol. 234, No. 1199, March 1956, pp. 528-541
%   Z.S. Alterman and D. Loewenthal, Algebraic Expressions for the
%                impulsive motion of an elastic half-space, Israel Journal
%                of Technology, Vol. 7, No. 6, 1969, pp. 495-504

% default data
if nargin<6, rho=1; end     % mass density
if nargin<5, Cs=1; end      % shear wave velocity
if nargin<4, h=1; end       % depth of source
if nargin<3, pois=0.25; end % Poisson's ratio       
N = 1180;   % number of time steps to arrival of PS waves

mu = rho*Cs^2;      % shear modulus
r1 = sqrt(x^2+(z-h)^2); % source-receiver distance
r2 = sqrt(x^2+(z+h)^2); % image source-receiver distance
s1 = x/r1;          % sin(theta1)
c1 = (z-h)/r1;      % cos(theta1)
s2 = x/r2;          % sin(theta2)
c2 = (h+z)/r2;      % cos(theta2)

a2 = (0.5-pois)/(1-pois);
a = sqrt(a2);                   % Cs/Cp
Cp = Cs/a;                      % P-wave velocity
tS = r1/Cs;                     % S-wave arrival (none here)
tP =  r1/Cp;                    % time of arrival of direct P waves
tPP = r2/Cp;                    % time of arrival of PP waves
tPS = t_PS (x,z,h,Cs,Cp);       % time of arrival or PS waves
dt =  tPS/N;                    % time step
t1 = (tP+dt):dt:tPP;            % time before reflections
t2 = (tPP+dt):dt:tPS;           % time from PP to PS reflection
t3 = (tPS+dt):dt:3*tPS;         % time after arrival of PS waves
% Find q(tau) from  tau(q) by solving quartic
X=x/r2; Z=z/r2; H=h/r2;
A = ((H+Z)^2+X^2)*((H-Z)^2+X^2);
B1 = X*(X^2+H^2+Z^2);
C1 = X^2*(a2*H^2+Z^2)+(H^2-Z^2)*(a2*H^2-Z^2);
C2 = 3*X^2+H^2+Z^2;
D1 = X*(a2*H^2+Z^2);
E1 = (a*H+Z)^2;
E2 = (a*H-Z)^2;
tau = t3*Cs/r2;                 % dimensionless time for PS waves
q3 = [];
for j=1:length(t3)
    tau2 = tau(j)^2;
    B = tau(j)*B1;
    C = tau2*C2-C1;
    D = tau(j)*(tau2*X-D1);
    E = (tau2-E1)*(tau2-E2);
    q = Ferrari(A,B,C,D,E);
    %q = roots([A,-4*i*B,-2*C,4*i*D,E]); % use in lieu of Ferrari?
    q = q(find(real(q)>=0 & imag(q)>=0)); % discard negative roots 
    [q1,I] = min(imag(q));  % find position of true root
    q3 = [q3,q(I)];         % choose that root
end
% Paco's approximation
R = (h+z/a)/(h+z);  % r_eq/r2
r3 = R*r2;          % equivalent radius
tapp = tau/R;
T = conj(sqrt(tapp.^2-a2)); % conj --> T must have neg. imag part
q3app = R*(c2*T+i*tapp*s2);
% Compare exact vs. Paco's 
plot(t3,real(q3));
hold on
plot(t3,imag(q3),'r');
plot(t3,real(q3app),'--');
plot(t3,imag(q3app),'r--');
grid on
title ('q3 --> exact vs. Paco''s approximation')
xlabel('Time')



% Find and plot the time histories
% ****************************

% a) From tP to tPP
T1 = sqrt(t1.^2-tP^2);
f1 = (0.5/r1)*t1./T1;
u1 = f1*s1;
w1 = f1*c1;

% b) From tPP to tPS
T1 = sqrt(t2.^2-tP^2);
T2 = sqrt(t2.^2-tPP^2);
f1 = (0.5/r1)*t2./T1;
f2 = (0.5/r2)*t2./T2;
q2 = (c2*T2+i*s2*t2)*Cs/r2;
dq2 = c2*t2./T2+i*s2;  % derivative
Q2 = q2.^2;
Q2S = sqrt(Q2+1);
Q2P = sqrt(Q2+a2);
S2 = (1+2*Q2).^2;
D2 = S2-4*Q2.*Q2S.*Q2P; % Rayleigh function
u2 = f1*s1-f2*s2-(4/r2)*imag(q2.^3.*Q2S.*dq2./D2);
w2 = f1*c1+f2*c2-(1/r2)*real(S2.*dq2./D2);

% c) From tPS and on
T1 = sqrt(t3.^2-tP^2);
T2 = sqrt(t3.^2-tPP^2);
f1 = (0.5/r1)*t3./T1;
f2 = (0.5/r2)*t3./T2;
% Contribution of PP waves
q2 = (c2*T2+i*s2*t3)*Cs/r2;
dq2 = c2*t3./T2+i*s2;  % derivative
Q2 = q2.^2;
Q2S = sqrt(Q2+1);
Q2P = sqrt(Q2+a2);
S2 = (1+2*Q2).^2;
D2 = S2-4*Q2.*Q2S.*Q2P; % Rayleigh function
f3 = (4/r2)*imag(q2.^3.*Q2S.*dq2./D2);
f5 = (1/r2)*real(S2.*dq2./D2);
% Contribution of PS waves
Q3 = q3.^2;
Q3S = sqrt(Q3+1);
Q3P = sqrt(Q3+a2);
S = 1+2*Q3;
S3 = S.^2;
D3 = S3-4*Q3.*Q3S.*Q3P; % Rayleigh function
dq3 = 1./((h/r2./Q3P+z/r2./Q3S).*q3-i*x/r2);
f4 = (2/r2)*imag(q3.*S.*Q3S.*dq3./D3);
f6 = (2/r2)*real(Q3.*S.*dq3./D3);
u3 = f1*s1-f2*s2-f3+f4;
w3 = f1*c1+f2*c2-f5+f6;

% Combine the results and plot
temp=[0:t1(2)-t1(1):tP];
time = [temp,t1,t2,t3]*Cs;
u = [zeros(1,length(temp)),u1,u2,u3]*(pi);
w = [zeros(1,length(temp)),w1,w2,w3]*(pi);
figure;
plot(time,u);
tit = sprintf('Horizontal displacements at x =%f z =%f' , x, z);
title(tit);
xlabel('t*Cs/r1');
ylabel('Ux*r1*\mu');
grid on
figure;
plot(time,w);
tit = sprintf('Vertical displacements at x =%f z =%f' , x, z);
title(tit);
xlabel('t*Cs/r1');
ylabel('Uz*r1*\mu');
grid on

save('GarvinResult.mat','w','u','time')
%--------------------------------------------------------------

function [tPS,xP,xS] = t_PS (x,z,h,Cs,Cp)
% Determines the total travel time of the PS reflection
% Arguments
% *********
%  x = range of receiver
%  z = depth of receiver
%  h = depth of source
%  Cs = S-wave velocity
%  Cp = P-wave velocity

if z==0, tPS=sqrt(x^2+h^2)/Cp; return; end
a = Cs/Cp;
% Bracket the S point
xP = x*h/(h+z);    % point of reflection of PP ray
ang1 = atan(xP/h); % minimum angle of incident ray
ang2 = atan(x/h);  % maximum angle
% Find the S point by search within bracket
dang = (ang2-ang1)/10;
TOL = 1.e-8;
TRUE = 1;
while TRUE
    angP = ang1+dang;
    angS = asin(a*sin(angP));
    L = h*tan(angP)+z*tan(angS);
    if L>x
        if L-x<TOL*x, break, else, dang = dang/10; end
    else
        ang1 = angP;
    end
end
tPS = h/cos(angP)/Cp+z/cos(angS)/Cs;
if nargout<3, return, end
xS = h*tan(angP);  % point of reflection of PS ray
return

%--------------------------------------------------------------

function [q] = Ferrari(A,B,C,D,E)
% Solves quartic equation
%   A*q^4 - 4*i*B*q^3 - 2*C*q^2 + 4*i*D*q + E
% by Ferrari's method

B = B/A; C=C/A; D=D/A; E=E/A;
a = 2*(3*B^2-C);
b = 4*i*(D-B*C+2*B^3);
c = E-4*B*D+2*C*B^2-3*B^4;
if b==0
    s2 = sqrt(a^2-4*c);
    s1 = sqrt((s2-a)/2);
    s2 = sqrt((-s2-a)/2);
    q1 = i*B+s1;
    q2 = i*B-s1;
    q3 = i*B+s2;
    q4 = i*B-s2;
else
    P = -(a^2/12+c)/3;
    Q = a/6*(c-a^2/36)-b^2/16;
    R = sqrt(Q^2+P^3)-Q;
    U = R^(1/3);
    V = U-5/6*a;
    if U==0
        V = V-(2*Q)^(1/3);
    else
        V = V-P/U;
    end
    W = sqrt(a+2*V);
    s2 = -(3*a+2*V);
    s1 = sqrt(s2-2*b/W);
    s2 = sqrt(s2+2*b/W);
    q1 = i*B+0.5*(W-s1);
    q2 = i*B+0.5*(-W+s2);
    q3 = i*B+0.5*(W+s1);
    q4 = i*B+0.5*(-W-s2);
end
q = [q1, q2, q3, q4].';
return
