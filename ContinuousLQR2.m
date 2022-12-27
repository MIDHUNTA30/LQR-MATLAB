% MATLAB program for implementing Continuous LQR
function main
clear all;
close all;
% initialize system and simulation parameters
A=[0 1;-2 1];
B=[0;1];
Tf=2;dt=0.1;N=Tf/dt;n=2;m=1;s=1;J0=0;
Qf=[1 0;0 1];Q=[1 0;0 1];R=1;
x0 = [10;5];

[Ka,Sa,Pa]=lqr(A,B,Q,R);


% Calling the LQR function
[K0,P0] = LQR(A,B,Qf,Q,R,N,dt,n,m);

% Initializing the vectors and matrices
x = zeros(n,N+1);x(:,1)=x0;
u = zeros(m,N);
  
% Simulating the system with LQR
for j=1:N   
      K=K0(j,:);
      %K=[2.735 -2.747];
      u(j)=-K*x(:,j);
      x(:,j+1)=x(:,j)+(A*x(:,j)+B*u(:,j))*dt;
      J0=J0+(x(:,j)'*Q*x(:,j)+u(:,j)'*R*u(:,j))*dt;
   end      
     J1=J0+ x(:,N+1)'*Qf*x(:,N+1)     
   
% Plotting the responces     
figure(1)
time = (0:N)*dt;
subplot(2,2,1)
plot(time,x(1,:),'k','LineWidth',1) 
hold on
plot(time,x(2,:),'r','LineWidth',1) 
legend('$x_{1}$','$x_{2}$','Interpreter','latex');
axis([0 Tf -15 15])
xlabel('t','Interpreter','latex');ylabel('$\textbf{x}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:Tf/4:Tf])
set(gca,'ytick',[-15:7.5:15])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,2,2)
plot(time(1:end-1),u,'k','LineWidth',1)
axis([0 Tf -30 30])
xlabel('t','Interpreter','latex');ylabel('$\textbf{u}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:Tf/4:Tf])
set(gca,'ytick',[-30:15:30])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,2,3)
plot(time(1:end-1),K0(:,1),'k','LineWidth',1) 
hold on
plot(time(1:end-1),K0(:,2),'r','LineWidth',1) 
legend('$K_{1}$','$K_{2}$','Interpreter','latex');
axis([0 Tf -1 5])
xlabel('t','Interpreter','latex');ylabel('$\textbf{K}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:Tf/4:Tf])
set(gca,'ytick',[-1:1.5:5])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,2,4)
plot(time,P0(:,1),'k','LineWidth',1) 
hold on
plot(time,P0(:,2),'r','LineWidth',1)
legend('$P_{11}$','$P_{22}$','Interpreter','latex');
axis([0 Tf 0 10])
xlabel('t','Interpreter','latex');ylabel('$\textbf{P}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:Tf/4:Tf])
set(gca,'ytick',[0:2.5:10])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg fig2
end


% LQR function
function [K0,P0] = LQR(A,B,Qf,Q,R,N,dt,n,m)
P=Qf;
P0(N+1,:)=diag(P);

  for k=N:-1:1    
     P = P + (A'*P+P*A-P*B*inv(R)*B'*P+Q)*dt; 
     K = inv(R)*B'*P;
     K0(k,:)=K;
     P0(k,:)=diag(P);
  end   

end