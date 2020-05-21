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
a44=2*beta*(d+f10)/numitor;
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

%% theta si psi
A1=[0 1 0 0; 0 a22 a23 a24; 0 0 0 1; 0 a42 a43 a44];
B1=[0;b21;0;b41];
C1=[1 0 0 0];
D1=[0];
eig(A1)
Co=ctrb(A1,B1)
rank(Co)
k1=acker(A1,B1,[-50 -50 -50 -50])
step(A1,B1,C1,D1);
figure
step(A1-B1*k1,B1,C1,D1);
A0=A1-B1*k1
F1=inv(C1/(-A0)*B1)
figure
step(A1-B1*k1,B1*F1,C1,D1*F1);
F11=2*F1;
F12=F1;

Ob=obsv(A1,C1)
rank(Ob)
Lt1=acker(A1',C1',[-1000 -1000 -1000 -1000])

%% reglare phi
A2=[0 1; 0 a66]; 
B21=[0;b61];
B2=B21;
C2=[1 0];
D2=[0];
eig(A2)
Co=ctrb(A2,B2)
rank(Co)

k3=acker(A2,B2,[-50 -50])
step(A2,B2,C2,D2);
figure
step(A2-B2*k3,B2,C2,D2);
A0=A2-B2*k3
F3=inv(C2/(-A0)*B2)
figure
step(A2-B2*k3,B2*F3,C2,D2*F3);

Ob=obsv(A2,C2)
rank(Ob)
Lt3=acker(A2',C2',[-1000 -1000])
%%
A=[0 1 0 0 0 0; 0 a22 a23 a24 0 0; 0 0 0 1 0 0; 0 a42 a43 a44 0 0; 0 0 0 0 0 1; 0 0 0 0 0 a66];
B=[0;b21;0;b41;0; b61];
C=[1 0 0 0 0 0];
C1=[0 0 0 0 1 0];
D=[0];
k=[k1 k3];
F=F1;
t=0:100;
sys=ss(A-B*k,B*F,C,D*F);
step(sys,t);
figure
sys1=ss(A-B*k,B*F,C1,D*F);
step(sys1,t);


