clear class;clear all;clc;
intv=intvclass2();

% a=intv.def(2,3)
% b=intv.def(4)
% c=intv.def([1 3;4 5],[0 2;-7 2])
% d=intv.def([3 6;-4 5],[2 -2;9 2])
% intv.isinterval(c)
% intv.def(c)
% 
% aunionb=intv.union(a,b)
% cuniond=intv.union(c,d)
% aplusb=intv.add(a,b)
% cplusd=intv.add(c,d)
% aminusb=intv.sub(a,b)
% cminusd=intv.sub(c,d)
% atimesb=intv.mult(a,b)
% ctimesd=intv.mult(c,d)
% adivb=intv.div(a,b)
% % acdivd=intv.div(c,d)
% % apower=intv.pow(a,2)
% % cpower=intv.pow(c,2)
% a=intv.int2mat(a)
% 
% A=[1 0 1;
%    4 5 6;
%    7 8 9];
% A=intv.def(A,A+5)
% B=[1 2;
%    3 4;
%    5 6]
% C=intv.mult(A,B)
% 
% BMI=intv.def(18.5,25)
% height=intv.def(1.68,1.71)
% weight=intv.mult(intv.mult(BMI,height),height)

a=intv.def(1,3)
b=intv.def(2,4)
x=intv.div(b,a)