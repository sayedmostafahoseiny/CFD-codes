clear all ; close all ; clc ;

% definition of parameters %
L = 0.2 ; % length of fin
R = 0.01 ; % radius of fin
h = 25 ; % convection heat transfer coefficient of ambient 
k = 200 ; % conduction heat transfer coefficient
Tinf = 300 ; % ambient temperature in kelvin
Tbase = 500 ; % temperature of base of fin
teta_base = 1 ;
P = 2*pi*R ; % perimeter of fin
A = pi*R^2 ; % cross section area of fin
m = (h*P/(k*A))^0.5 ;
n = 2560; % number of divisions
deltax = 1/n ;
S=(2+(m*L*deltax)^2);
B=-(S+(2*h*L*deltax/k));


% assembling coefficient matrix A
A=zeros(n+1,n+1);% since we have n+1 nodes, we need n+1 equations
A(1,1)=1;% left boundary node
for i=2:n ; % interior nodes
    A(i,i)= -S ;
    A(i,i-1)=1 ;
    A(i,i+1)=1 ;
end
% right boundary node
A(n+1,n)= 2 ;
A(n+1,n+1)=B ;

% assembling right hand side matrix
B = zeros(n+1,1) ;
B(1,1)=teta_base ;% left boundary condition

% solving system of equations
solution=gmres(A,B,[],1e-7,1000);
theta=solution' ;
x=linspace(0,L,n+1);
figure (1)
plot1=plot(x,theta,'o')

% analytical solution
xx=linspace(0,L,length(x));
teta_Analytic=(cosh(m*(L-xx))+(h/(m*k))*sinh(m*(L-xx)))/(cosh(m*L)+(h/(m*k))*sinh(m*L));
hold on
plot2=plot(xx,teta_Analytic,'r','linewidth',2)
grid on
title( ' Temperature distribution ' )
xlabel( ' x (m) ' )
ylabel( ' theta ' )
legend ([plot1 plot2],{'numerical (n=50)','Analitical solution'}) % validation

% Error Analysis
L1_norm=sum(abs(teta_Analytic-theta))/(n+1) % calculation of first norm of Error

% post processing
T=theta*(Tbase-Tinf)+Tinf; % temperature distribution in (C) in fin
q=-k*((T(2)-T(1))/(x(2)-x(1))) % the flux which enter the fin from base
eta=(pi*q*R^2)/(h*2*pi*R*L*(Tbase-Tinf))*100 % fin efficiency
eps = (pi*q*R^2)/(pi*R^2*h*(Tbase-Tinf)) % fin performance coefficient
T_ave=sum(T)/length(T)