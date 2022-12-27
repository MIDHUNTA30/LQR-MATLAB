% MATLAB program for implementing Discrete LQR
function main
clear all;
close all;
% initialize system and simulation parameters
A=[0.5 0;-1 1.5];
B=[0.5;0.1];
N=50;n=2;m=1;s=1;J0=0;
Qf=[1 0;0 1];Q=[1 0;0 1];R=1;
x0 = [10;5];
 


% Calling the LQR function
[K0,P0] = LQR(A,B,Qf,Q,R,N,n,m);

% Initializing the vectors and matrices
x = zeros(n,N+1);x(:,1)=x0;
u = zeros(m,N);
  
% Simulating the system with LQR
for j=1:N   
      K=K0(j,:);
      %K=[2.735 -2.747];
      u(j)=-K*x(:,j);
      x(:,j+1)=A*x(:,j)+B*u(:,j);
      J0=J0+x(:,j)'*Q*x(:,j)+u(:,j)'*R*u(:,j);
   end      
     J1=J0+ x(:,N+1)'*Qf*x(:,N+1)     
   
% Plotting the responces     
figure(1)
time = (0:N);
subplot(2,2,1)
plot(time,x(1,:),'k.-','LineWidth',1) 
hold on
plot(time,x(2,:),'r.-','LineWidth',1) 
legend('$x_{1}$','$x_{2}$','Interpreter','latex');
axis([0 N -10 10])
xlabel('k','Interpreter','latex');ylabel('$\textbf{x}_{k}$','Interpreter','latex');
grid on
set(gca,'xtick',[0:N/5:N])
set(gca,'ytick',[-10:5:10])
ax = gca;
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,2,2)
plot(time(1:end-1),u,'k.-','LineWidth',1)
axis([0 N -15 5])
xlabel('k','Interpreter','latex');ylabel('$\textbf{u}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:N/5:N])
set(gca,'ytick',[-15:5:5])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,2,3)
plot(time(1:end-1),K0(:,1),'k.-','LineWidth',1) 
hold on
plot(time(1:end-1),K0(:,2),'r.-','LineWidth',1) 
legend('$K_{1}$','$K_{2}$','Interpreter','latex');
axis([0 N -4 4])
xlabel('k','Interpreter','latex');ylabel('$\textbf{K}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:N/5:N])
set(gca,'ytick',[-4:2:4])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
subplot(2,2,4)
plot(time,P0(:,1),'k.-','LineWidth',1) 
hold on
plot(time,P0(:,2),'r.-','LineWidth',1)
legend('$P_{11}$','$P_{22}$','Interpreter','latex');
axis([0 N 0 40])
xlabel('k','Interpreter','latex');ylabel('$\textbf{P}_{k}$','Interpreter','latex');
grid on
ax = gca;
set(gca,'xtick',[0:N/5:N])
set(gca,'ytick',[0:10:40])
ax.GridAlpha = 1
ax.GridLineStyle = ':'
print -dsvg fig1
end


% LQR function
function [K0,P0] = LQR(A,B,Qf,Q,R,N,n,m)
P=Qf;
P0(N+1,:)=diag(P);

  for k=N-1:-1:0                                    
     K = inv(R +B'* P*B)*(B'* P*A);
     P = Q + K'*R*K+(A-B*K)'* P*(A-B*K); 
     K0(k+1,:)=K;
     P0(k+1,:)=diag(P);
  end   

end