close all
clear all
clc

syms A B C D E F G H I J x y z p q r s m

E1 = 2*(x*C-y*z*D) - (A^2*s+B^2*r);
E2 = (F*x-y*z*G) + (C*x+D*y*z)-(A*B*s+B*E*r);
E3 = (H*x-y*z*I) - (A*C*s+B*F*r);
E4 = (I*x-y*z*J) - (A*D*s+B*G*r);

E5 = 2*(y*F+x*z*G) - (B^2*s+E^2*r);
E6 = (H*y+x*z*I) - (C*B*s+E*F*r);
E7 = (I*y+x*z*J) - (B*D*s+E*G*r);

E8 = p - (C^2*s+F^2*r);
E9 = -(C*D*s+F*G*r);

E10 = q - (D^2*s+G^2*r);

vars = [A, B, C, D, E, F, G, H, I, J];
E = [E1, E2, E3, E4, E5, E6, E7, E8, E9, E10];

[S1, S2, S3, S4, S5, S6, S7, S8, S9, S10] = solve(E, vars)