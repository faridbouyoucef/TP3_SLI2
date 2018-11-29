clear all 
close all
clc

S=0.0154;
Sn=5*10^-5;
g=9.81;
H10=0.27474;
H20=0.0299;
H30=0.1368;
H00=0;
a13=0.4753*Sn*sqrt(2*g);
a32=0.4833*Sn*sqrt(2*g);
a20=0.9142*Sn*sqrt(2*g);

R13=(2*sqrt(H10-H30))/a13;
R32=(2*sqrt(H30-H20))/a32;
R20=(2*sqrt(H20-H00))/a20;

A=[-1/(S*R13) 1/(S*R13) 0;
    1/(S*R13) -(1/S)*((1/R13)+(1/R32)) 1/(S*R32);
    0 1/(S*R32) -(1/S)*((1/R32)+(1/R20))]
B=[ 1/S; 0; 0]
C=[1 0 0];
D=[0];

sys=ss(A,B,C,D);
Vp=eig(A);
Co=ctrb(sys)
rang_co=rank(Co);
% le rang=3=dimension du système, alors le système est commandable
P=[-0.03333 -0.03777 -0.0444]; 

K=acker(A,B,P)
Abf=A-B*K;
N=0.98/(-(Abf+B*K)*B*1)

sys=ss(Abf,B,C,D);
Vp=eig(A)




%********partie 4 *********

%%%Q1%%%%%%
OBSV=obsv(sys)
rang_co=rank(Co)


%%Q2%%%
P1=[-0.0333*5 -0.3777*5 -0.0444*5];
G=acker(A',C',P1)'

A22=[A(2,2) A(2,3);A(3,2) A(3,3)] 
A12=[A(1,2) A(1,3)]   



 F=[-0.3777*5 ; -0.0444*5]
 %inverse= inv(A12);

 Gind=acker(A22',A12',[-0.3777*5 -0.0444*5])'
 A21=[A(2,1); A(3,1) ]
%Gtild=(F*Gind)-(Gind*A(1,1))+A21











