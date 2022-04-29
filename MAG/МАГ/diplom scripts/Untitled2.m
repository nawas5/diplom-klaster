close all
clear all
clc

syms A B C D E F G H I J K L M N O x y z p q r s t
%eqns = [2*(C*y-D*x*z)-(A^2*p+B^2*q+E^2*r) == 0,    2*(G*y-H*x*z)-(A*B*p+B*F*q+E*I*r) == 0,    2*(J*y-K*x*z)-(A*C*p+B*G*q+E*L*r) == 0, 2*(K*y-M*x*z)-(A*D*p+B*H*q+E*N*r) == 0, 2*(L*y-N*x*z)-(A*E*p+B*I*q+E*O*r) == 0,    2*(C*x+D*y*z)-(B*A*p+F*B*q+I*E*r) == 0,    2*(G*x+H*y*z)-(B*B*p+F*F*q+I*I*r) == 0,    2*(J*x+K*y*z)-(B*C*p+F*G*q+I*L*r) == 0,    2*(K*x+M*y*z)-(B*D*p+F*H*q+I*N*r) == 0,    2*(L*x+N*y*z)-(B*E*p+F*I*q+I*O*r) == 0,    -(C*A*p+G*B*q+L*E*r) == 0,    -(C*B*p+G*F*q+L*I*r) == 0,     2*s-(C*C*p+G*G*q+L*L*r) == 0,    -(C*D*p+G*H*q+L*N*r) == 0,    -(C*E*p+G*I*q+L*O*r) == 0,    2*E-(D*A*p+H*B*q+N*E*r) == 0,    2*I-(D*B*p+H*F*q+N*I*r) == 0,    2*L-(D*C*p+H*G*q+N*L*r) == 0,    2*N-(D*D*p+H*H*q+N*N*r) == 0,    2*O-(D*E*p+H*I*q+N*O*r) == 0,    -(E*A*p+I*B*q+O*E*r) == 0,    -(E*B*p+I*F*q+O*I*r) == 0,    -(E*C*p+I*G*q+O*L*r) == 0,    -(E*D*p+I*H*q+O*N*r) == 0, 2*t-(E*E*p+I*I*q+O*O*r) == 0];
%eqns1 = [eqns(1) eqns(2)];



eqns1 = 2*(C*y-D*x*z)-(A^2*p+B^2*q+E^2*r);
eqns2 = (G*y-H*x*z)+(C*x+y*z*D)-(A*B*p+B*F*q+E*I*r);
eqns3 = (J*y-K*x*z)-(A*C*p+B*G*q+E*L*r);
eqns4 = (K*y-M*x*z)+E-(A*D*p+B*H*q+E*N*r);
eqns5 = (L*y-N*x*z)-(A*E*p+B*I*q+E*O*r);

eqns7 = 2*(G*x+H*y*z)-(B*B*p+F*F*q+I*I*r);
eqns8 = (J*x+K*y*z)-(B*C*p+F*G*q+I*L*r);
eqns9 = (K*x+M*y*z)+I-(B*D*p+F*H*q+I*N*r);
eqns10 = (L*x+N*y*z)-(B*E*p+F*I*q+I*O*r);

eqns13 = s-(C^2*p+G^2*q+L^2*r);
eqns14 = L-(C*D*p+G*H*q+L*N*r);
eqns15 = -(C*E*p+G*I*q+L*O*r);

eqns19 = 2*N-(D^2*p+H^2*q+N^2*r);
eqns20 = O-(D*E*p+H*I*q+N*O*r);

eqns25 = t-(E^2*p+I^2*q+O^2*r);

vars = [A, B, C, D, E, F, G, H, I, J, K, L, M, N, O];
E1 = [eqns1, eqns2, eqns3, eqns4, eqns5, eqns7, eqns8, eqns9, eqns10, eqns13, eqns14, eqns15, eqns19, eqns20, eqns25];
[solA, solB, solC, solD, solF, solE, solG, solH, solI, solJ, solK, solL, solM, solN, solO] = solve(E1, vars);
%S = solve(eqns1, eqns2, eqns3, eqns4, eqns5, eqns7, eqns8, eqns9, eqns10, eqns13, eqns14, eqns15, eqns19, eqns20, eqns25);


