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
poles=[-0.5 -0.4 -0.3 -0.2 -0.6 -0.7];
K=place(A,B,poles*10)
K_opt=lqr(A,B,eye(6),eye(1));

L=place(A',C',poles)'
%% kalman filter
Ts=0.001
Ad = eye(6) + Ts*A;
Bd = Ts*B;
Cd=C;

Kd_optimal = dlqr(Ad, Bd, eye(6), eye(1));

xd = [0.9;0.8;0.7;0.6;0.5;0.4]; 
x_hat = [0;0;0;0;0;0]; 
xPred = [0;0;0;0;0;0];
u=[0;0];

alfa_1 = 0.9;
Q = alfa_1*eye(6)*1e-5;
alfa_2 = 0.1;
R = alfa_2*eye(1)*1e-5;

P=99999*eye(6);
I=eye(6);

for k=2:100
    w = sqrt(Q)*randn(6,1);  
    v = sqrt(R)*randn(1,1); 
    
    u(:,k-1)=-Kd_optimal*x_hat(:,k-1);
    xd(:,k)=Ad*xd(:,k-1)+Bd*u(:,k-1)+w;
    y(:,k)=Cd*xd(:,k-1)+v;
    
    xPred(:,k)=Ad*xPred(:,k-1)+Bd*u(:,k-1);
    pPred=Ad*P*Ad'+Q;
    
    Kk=pPred*Cd'*inv(Cd*pPred*Cd'+R);
    x_hat(:,k)=xPred(:,k)+Kk*(y(:,k)-Cd*xPred(:,k));
    P=(I-Kk*Cd)*pPred*(I-Kk*Cd)'+Kk*R*Kk';
end
figure, plot(x_hat'-xd')
%figure, plot(xd')
%figure, plot(x_hat')
%% discretizare sist liniar
Ts=0.001;
Ad = eye(6) + Ts*A
Bd = Ts*B
Cd=C

Kd_optimal = dlqr(Ad, Bd, eye(6), eye(1));
Ld_optimal = dlqr(Ad',Cd',eye(6),1)';

eig(Ad-Bd*Kd_optimal)

xd = [0.9;0.8;0.7;0.6;0.5;0.4]; 
x_hat = [0;0;0;0;0;0];
u(:,1)=-Kd_optimal*xd(:,1);
y(:,1)=Cd*xd(:,1);
y_hat(:,1)=Cd*x_hat(:,1);

for k=2:500
    u(:,k) = -Kd_optimal*xd(:,k-1);
    xd(:,k) = Ad*xd(:,k-1) + Bd*u(:,k);
    y(:,k) = Cd*xd(:,k);
    y_hat(:,k) = Cd*x_hat(:,k-1);
    x_hat(:,k) = Ad*x_hat(:,k-1) + Bd*u(:,k) + Ld_optimal*(y(:,k) - y_hat(:,k)); 
end
figure, plot(xd'-x_hat')

%% unknown observer

E = [1;0;0;0;0;0];
Ae = [A E;0 0 0 0 0 0 0]
Be = [B; 0 0]
Ce=[C [0;0;0;0;0;0]]

Co_ext=ctrb(Ae,Be)
rank(Co_ext)

Ob_ext=obsv(Ae,Ce);
rank(Ob_ext);

poles = [-0.1 -0.2 -0.3 -0.4 -0.5 -0.6 -0.7];
%K=place(A,B,[-3, -2, -5, -4, -6, -2]);
L_extins=place(Ae',Ce',poles)'

% decoupling of an unknown input state
rank(E)
rank(C*E)
H=E*pinv(C*E)
T=eye(6)-H*C

Ob_dec=obsv(T*A,C)
rank(Ob_dec) %observabila

K1=place((T*A)',C',[-50, -15, -11, -13, -16, -12])'

F=T*A-K1*C
K2=F*H

K_total=K1+K2