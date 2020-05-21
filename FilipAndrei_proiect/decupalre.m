m=0.32; %weght of the wheel
R=0.15/2; %diameter of the wheel
Jw=0.0013;%Moment of inertia of the wheel
M=5.41;%weight of the vehicle
W=0.4;%width of the vehicle
L=0.102;%height of the mass center of vehicle
Jpsi=0.104;%moment of inertia of the tilt axis vehicle
Jphi=0.0484;%moment of inertia of the vehicle related to the axis of rotation
Jm=0.00119;%Identified moment of inertia of the DC motor and gearbox taking into account gearbox ratio
RDC=1;%Resistance of the winding of the DC motor
Kt=0.025;%Torque constant of the DC motor
Ke=0.025;%Voltage constant of the DC motor
fm=0.00024;%Identified friction coefficient between the vehicleand DC motor
c1=(2*m+M)*R*R + 2*Jw+2*Jm;
c2=M*L*L+Jpsi+2*Jm;
f10=M*R*L-2*Jm;
d1=Kt/RDC;
d2=Kt*Ke/RDC+fm;
beta=d2;
alfa=d1;
d=c1;
c=c2;
g=9.81;
numitor=d*c-f10*f10;
a22=-2*beta*(c+f10)/numitor;
a23=-M*g*L*f10/numitor;
a24=2*beta*(c+f10)/numitor;
a42=2*beta*(d+f10)/numitor;
a43=d*M*g*L/numitor;
a44=-2*beta*(d+f10)/numitor;
b21=alfa*(c+f10)/numitor;
b22=alfa*(c+f10)/numitor;
b41=-alfa*(d+f10)/numitor;
b42=-alfa*(d+f10)/numitor;
a66=((-W*W*beta)/(2*R*R))/(m*W*W/2+W*W/(2*R*R)*(Jw+Jm)+Jphi);
b61=((-W*alfa)/(2*R))/(m*W*W/2+W*W/(2*R*R)*(Jw+Jm)+Jphi);
b62=((-W*alfa)/(2*R))/(m*W*W/2+W*W/(2*R*R)*(Jw+Jm)+Jphi);
A1=[0 1 0 0; 0 a22 a23 a24; 0 0 0 1; 0 a42 a43 a44];
B1=[0 0;b21 b22; 0 0; b41 b42];
A2=[0 1; 0 a66]; 
B2=[0 0; b61 b62];



%% simulare si identificare prima intrare
A1=[0 1 0 0; 0 a22 a23 a24; 0 0 0 1; 0 a42 a43 a44];B1=[0;b21;0;b41];C1=[1 0 0 0];D1=[0];
k1=acker(A1,B1,[-5 -5 -5 -5]);
A0=A1-B1*k1;
F1=inv(C1/(-A0)*B1);

A2=[0 1; 0 a66]; B2=[0;b61];C2=[1 0];D2=[0];
k3=acker(A2,B2,[-1 -1]);
A0=A2-B2*k3;
F3=inv(C2/(-A0)*B2);

A=[0 1 0 0 0 0; 0 a22 a23 a24 0 0; 0 0 0 1 0 0; 0 a42 a43 a44 0 0; 0 0 0 0 0 1; 0 0 0 0 0 a66];
B=[0;b21;0;b41;0; b61];
C=[1 0 0 0 0 0];
C1=[0 0 0 0 1 0];
D=[0];
k=[k1 k3];
F=F1;

[num11,dem11]=ss2tf(A-B*k,B*F*0.25,C,D*F*0.25);
G11=tf(num11,dem11);
[num11,dem11]=tfdata(G11,'v');

[num21,dem21]=ss2tf(A-B*k,B*F*0.25,C1,D*F*0.25);
G21=tf(num21,dem21);
[num21,dem21]=tfdata(G21,'v');


t=0:0.05:10;
y11=step(G11,t);
yst=0.26;
ust=1;
k=yst/ust;
ymax=0.286;
sigma=(ymax-yst)/yst;
zeta=-log(sigma)/(sqrt(pi*pi+log(sigma)*log(sigma)));
tr=2.65;
wn=4/zeta/tr;
H11=tf([-0.35 k*wn*wn],[1 2*zeta*wn wn*wn]);
hold on
ysim11=step(H11,t);
plot(t,[y11 ysim11]);
J11=norm(y11-ysim11);
EMPN11=J11/(norm(y11-mean(y11)))

figure
t=0:0.05:50;
y21=step(G21,t);
yst=-0.0964;
H21=tf([4 -23 -0.97],[1 15 50 100 10]);
hold on
ysim21=step(H21,t);
plot(t,[y21 ysim21]);
J21=norm(y21-ysim21);
EMPN21=J21/(norm(y21-mean(y21)))
%% simulare si identificare a doua intrare
A1=[0 1 0 0; 0 a22 a23 a24; 0 0 0 1; 0 a42 a43 a44];B1=[0;b21;0;b41];C1=[1 0 0 0];D1=[0];
k1=acker(A1,B1,[-10 -10 -10 -10]);
A0=A1-B1*k1;
F1=inv(C1/(-A0)*B1);

A2=[0 1; 0 a66]; B2=[0;b61];C2=[1 0];D2=[0];
k3=acker(A2,B2,[-1 -1]);
A0=A2-B2*k3;


A=[0 1 0 0 0 0; 0 a22 a23 a24 0 0; 0 0 0 1 0 0; 0 a42 a43 a44 0 0; 0 0 0 0 0 1; 0 0 0 0 0 a66];
B=[0;b21;0;b41;0; b61];
C=[1 0 0 0 0 0];
C1=[0 0 0 0 1 0];
D=[0];
k=[k1 k3];
F=F1;

[num12,dem12]=ss2tf(A-B*k,B*F*0.75,C,D*F*0.75);
G12=tf(num12,dem12);
[num12,dem12]=tfdata(G12,'v');

[num22,dem22]=ss2tf(A-B*k,B*F*0.75,C1,D*F*0.75);
G22=tf(num22,dem22);
[num22,dem22]=tfdata(G22,'v');

t=0:0.05:2;
y12=step(G12,t);
yst=0.75;
ust=1;
k=yst/ust;
ymax=0.79;
sigma=(ymax-yst)/yst;
zeta=-log(sigma)/(sqrt(pi*pi+log(sigma)*log(sigma)));
tr=0.82;
wn=4/zeta/tr;
H12=tf([-7.5 k*wn*wn],[1 2*zeta*wn wn*wn]);
hold on
ysim12=step(H12,t);
plot(t,[y12 ysim12]);
J12=norm(y12-ysim12);
EMPN12=J12/(norm(y12-mean(y12)))

figure
t=0:0.05:100;
y22=step(G22,t);
yst=-0.282;
H22=tf([160 -1000 -42.3],[1 35 360 1450 150]);
hold on
ysim22=step(H22,t);
plot(t,[y22 ysim22]);
J22=norm(y22-ysim22);
EMPN22=J22/(norm(y22-mean(y22)))

%% Rereuri dem interdependenta
GG11=tf([-0.3 1.5],[1 3 6.25]);
GG21=tf([-0.6],[1 3 6.25]); % 1 3.6 9
GG12=tf([-7.5 37.5],[1 10 50]);
GG22=tf([-14.25],[1 10 50]); %1 13.6 65.6
[numgg11,demgg11]=tfdata(GG11,'v');
[numgg12,demgg12]=tfdata(GG12,'v');
[numgg21,demgg21]=tfdata(GG21,'v');
[numgg22,demgg22]=tfdata(GG22,'v');


detG=GG11*GG22-GG12*GG21;
[numdetG,demdetG]=tfdata(detG,'v');
D11=1;D22=1;
D12=minreal(-GG12/GG11);
D21=minreal(-GG21/GG22);
[numD21,demD21]=tfdata(D21,'v');
[numD12,demD12]=tfdata(D12,'v');

M11=minreal(detG/GG22);
[numM11,demM11]=tfdata(M11,'v');
M22=minreal(detG/GG21);
[numM22,demM22]=tfdata(M22,'v');

step(M11*1000/(-12.5))
hold on 
step(M22*1000/15)
%% DECUPALRE pe H-uri
detG=H11*H22-H12*H21;
[numdetG,demdetG]=tfdata(detG,'v');
D11=1;D22=1;
D12=minreal(-H12/H11);
D21=minreal(-H21/H22);
[numD21,demD21]=tfdata(D21,'v');
[numD12,demD12]=tfdata(D12,'v');

M11=minreal(series(1/H22,detG));
[numM11,demM11]=tfdata(M11,'v');
M22=minreal(series(1/H11,detG));
[numM22,demM22]=tfdata(M22,'v');


t=0:0.01:3;
Ho=tf(1,[0.1 1]);
Hc11=minreal(Ho/(1-Ho)/M11);
[numcon11,demcon11]=tfdata(Hc11,'v');
Hd11=minreal(Hc11*M11);
H011=feedback(Hd11,1);
step(H011,t)

figure
Ho=tf(400,[1 20 400]);
Hc22=minreal(Ho/(1-Ho)/M22);
[numcon22,demcon22]=tfdata(Hc22,'v');
Hd22=minreal(Hc22*M22);
H022=feedback(Hd22,1);
step(H022,t)
%% num/dem Hs
[numh11,demh11]=tfdata(H11,'v');
[numh12,demh12]=tfdata(H12,'v');
[numh21,demh21]=tfdata(H21,'v');
[numh22,demh22]=tfdata(H22,'v');
%% DECUPALRE pe G-uri
detG=G11*G22-G12*G21;
[numdetG,demdetG]=tfdata(detG,'v');
D11=1;D22=1;
D12=-G12/G11;
D21=-G21/G22;
[numD21,demD21]=tfdata(D21,'v');
[numD12,demD12]=tfdata(D12,'v');

M11=series(1/G22,detG);
[numM11,demM11]=tfdata(M11,'v');
M22=series(1/G11,detG);
[numM22,demM22]=tfdata(M22,'v');

Ho=tf([625],[1 20 150 500 625]);
t=0:0.01:5;

Hc11=Ho/(1-Ho)/M11;
[numcon11,demcon11]=tfdata(Hc11,'v');
Hd11=minreal(Hc11*M11);
H011=feedback(Hd11,1);
step(H011,t)

figure

Hc22=Ho/(1-Ho)/M22;
[numcon22,demcon22]=tfdata(Hc22,'v');
Hd22=minreal(Hc22*M22);
H022=feedback(Hd22,1);
step(H022,t)


