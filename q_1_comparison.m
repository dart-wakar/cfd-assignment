%claring the workspaces
clc;%clear all;close all;

%defining the variables
xmin=0;    %minimum value of x
xmax=1;    %maximum value of x
N=100;     %no of nodes 
dt=0.001;  %timestep
t=0;       %time
tmax=1;    %maximum value of t 
a=1;       %velocity


%descritize the domain 
dx=(xmax-xmin)/(N-1);
x=xmin: dx: xmax;

A=load('lax_friedrich.txt');
[m,n]=size(A);
B=load('exact.txt');
C=load('lax_wendroff.txt');
for i=1:m

    u=A(i,:);
    wend_u=C(i,:);
    exact_u=B(i,:);
    %plot solution 
    plot(x,u,'ro','markerfacecolor','r');
    axis([0,1,-1,1]);
    hold on
    plot(x,exact_u,'b-','markerfacecolor','b');
    %hold on
    %plot(x,wend_u,'go','markerfacecolor','g');
    hold off
    shg;
    pause(dt);
    
end