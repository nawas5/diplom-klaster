close all
clear all
clc

syms A B C D E F x y p q s

E1 = 2*x*C - (A^2*p+B^2*q);
E2 = (E*x+C*y) - (A*B*p+B*D*q);
E3 = x*F - (A*C*p+B*E*q);
E4 = 2*y*E - (B^2*p+D^2*q);
E5 = y*F - (B*C*p+D*E*q);
E6 = s - (C^2*p+E^2*q);

vars = [A, B, C, D, E, F];
E = [E1, E2, E3, E4, E5, E6];

[S1, S2, S3, S4, S5, S6] = solve(E, vars)