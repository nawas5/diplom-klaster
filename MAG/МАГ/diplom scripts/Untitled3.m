close all
clear all
clc

syms A B C D E F G H I J x y z p q s

E1 = 2*(y*C-x*z*D) - (A^2*p+B^2*q);
E2 = (F*y-x*z*G)+(C*x+D*y*z)-(A*B*p+B*E*q);
E3 = (H*y-x*z*I) - (A*C*p+B*F*q);
E4 = (I*y-x*z*J) - (A*D*p+B*G*q);

E5 = 2*(x*F+y*z*G) - (B^2*p+E^2*q);
E6 = (H*x+y*z*I) - (C*B*p+E*F*q);
E7 = (I*x+y*z*J) - (B*D*p+E*G*q);

E8 = s - (C^2*p+F^2*q);
E9 = -(C*D*p+F*G*q);

E10 = -(D^2*p+G^2*q);

vars = [A, B, C, D, E, F, G, H, I, J];
E = [E1, E2, E3, E4, E5, E6, E7, E8, E9, E10];

[S1, S2, S3, S4, S5, S6, S7, S8, S9, S10] = solve(E, vars)