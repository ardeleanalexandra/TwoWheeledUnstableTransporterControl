clear all
clc
A=[0 1 0 0 0 0;
    0 -0.4733 -29.1464 0.4733 0 0;
    0 0 0 1 0 0;
    0 0.1401 29.1651 -0.1401 0 0;
    0 0 0 0 0 1;
    0 0 0 0 0 -2.6349]
B=[0 0;
    109.3664 109.3664;
    0 0;
    -32.3752 -32.3752;
    0 0;
    -0.5055 0.5055]

C=eye(6)
D=zeros(6,2)
sys=ss(A,B,C,D);
Hp=zpk(tf(sys))

C0=ctrb(A,B)
ctrb=rank(C0)
O0=obsv(A,C)
obs=rank(O0)
G11=Hp(1,1);
G12=Hp(2,1);
G13=Hp(3,1);
G14=Hp(4,1);
G15=Hp(5,1);
G16=Hp(6,1);
G21=Hp(1,2);
G22=Hp(2,2);
G23=Hp(3,2);
G24=Hp(4,2);
G25=Hp(5,2);
G26=Hp(6,2);

G_11=minreal(zpk(G11))
G_12=minreal(zpk(G13))
G_21=minreal(zpk(G23))
G_22=minreal(zpk(G21))

%alegeri bune:
%g11 g15 g25 g21
%g11 g13 g23 g21

[num_11,den_11]=tfdata(G_11,'V');
[num_12,den_12]=tfdata(G_12,'V');
[num_21,den_21]=tfdata(G_21,'V');
[num_22,den_22]=tfdata(G_22,'V');

G=[G_11 G_12;
   G_21 G_22]

%% 
D11=1;
D12=-G_12/G_11;
D21=-G_21/G_22;
D22=1;
D=[D11 D12;
   D21 D22]
[num_d11,den_d11]=tfdata(G_11,'V');
[num_d12,den_d12]=tfdata(G_12,'V');
[num_d21,den_d21]=tfdata(G_21,'V');
[num_d22,den_d22]=tfdata(G_22,'V');
H=minreal(zpk(G*D))

H11=H(1,1)
H22=H(2,2)

G1=minreal(zpk(H11),0.1)
G2=minreal(zpk(H22),0.1)

[num_G1,den_G1]=tfdata(G1,'V');
[num_G2,den_G2]=tfdata(G2,'V');
%% Calcul controller 1 Guillemin-Truxal
sigma=0.1;
tr=2;

% =>
zeta=(-log(sigma))/(sqrt(pi^2+(log(sigma))^2))
wn=4/(zeta*tr)
cv=wn/(2*zeta);
wb=wn*sqrt(1-2*zeta^2+sqrt(2-4*zeta^2+4*zeta^4));

G01=tf(wn^2,[1 2*zeta*wn wn^2])
Gc1=minreal(zpk(1/G1*G01/(1-G01)))


[num_gc1,den_gc1]=tfdata(Gc1,'V')

Gs1=series(Gc1,G1)
figure
bode(Gs1)
Gc01=minreal(feedback(Gs1,1))
[num_gc01,den_gc01]=tfdata(Gc01,'V')
figure
step(Gc01)
%overshoot=10 settling time=1.75
%% Calcul controller 2 Guillemin Truxal
sigma=0.1;
tr=2;

% =>
zeta=(-log(sigma))/(sqrt(pi^2+(log(sigma))^2))
wn=4/(zeta*tr)
cv=wn/(2*zeta);
wb=wn*sqrt(1-2*zeta^2+sqrt(2-4*zeta^2+4*zeta^4));

G02=tf(wn^2,[1 2*zeta*wn wn^2])

Gc2=minreal(zpk(1/G2*G02/(1-G02)))

[num_gc2,den_gc2]=tfdata(Gc2,'V');

Gs2=series(Gc2,G2)
figure
bode(Gs2)
Gc02=minreal(feedback(Gs2,1))
figure
step(Gc02)



