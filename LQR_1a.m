% MATLAB program for implementing linear quadratic regulator (LQR)
function main
clear all;
close all;
% initialize system and simulation parameters
global A B   
N=5;
Qf=[1 0;0 1];
Q=[1 0;0 1]; J0=0;
R=1;s=1; n=2;
x0 = [10;5]; 
A=[0.5 0;-1 1.5];
B=[0.5;0.1];


% Calling the LQR function
[K0,P0] = LQR(A,B,Qf,Q,R,N,n);

% Initializing the vectors and matrices
x = zeros(2,N+1);
x(:,1)=x0;
u = zeros(1,N);
  
% Simulating the system with LQR
for j=1:N   
      K=K0(j,:);
      %K=[2.735 -2.747];
      u(j)=-K*x(:,j);
      x(:,j+1)=A*x(:,j)+B*u(1,j);
      J0=J0+x(:,j)'*Q*x(:,j)+u(j)'*R*u(j);
   end      
     J1=J0+ x(:,N+1)'*Qf*x(:,N+1)     
   
% Plotting the responces     
figure(1)
time = (0:N);
subplot(4,1,1)
plot(time,x(1,:),'k.-','LineWidth',1) 
hold on
plot(time,x(2,:),'r.-','LineWidth',1) 
legend('$x_{1}$','$x_{2}$','Interpreter','latex');
axis([0 5 -10 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{x}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(4,1,2)
plot(time(1:end-1),u,'k.-','LineWidth',1)
axis([0 5 -13 2])
xlabel('k','Interpreter','latex');ylabel('$\textbf{u}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(4,1,3)
plot(time(1:end-1),K0(:,1),'k.-','LineWidth',1) 
hold on
plot(time(1:end-1),K0(:,2),'r.-','LineWidth',1) 
legend('$K_{1}$','$K_{2}$','Interpreter','latex');
axis([0 5 -3 3])
xlabel('k','Interpreter','latex');ylabel('$\textbf{K}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(4,1,4)
plot(time,P0(:,1),'k.-','LineWidth',1) 
hold on
plot(time,P0(:,4),'r.-','LineWidth',1)
legend('$P_{11}$','$P_{22}$','Interpreter','latex');
axis([0 5 0 20])
xlabel('k','Interpreter','latex');ylabel('$\textbf{P}_{k}$','Interpreter','latex');
grid on
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg fig5
end


% LQR function
function [K0,P0] = LQR(A,B,Qf,Q,R,N,n)
P=Qf;
P0(N+1,1:2)=P(1,1:2);
P0(N+1,3:4)=P(2,1:2);
    for k=N-1:-1:0
                                    
     K = inv(R +B'* P*B)*(B'* P*A);
     P = Q + K'*R*K+(A-B*K)'* P*(A-B*K); 
     K0(k+1,:)=K;
     P0(k+1,1:2)=P(1,1:2);
     P0(k+1,3:4)=P(2,1:2);
    end   

end