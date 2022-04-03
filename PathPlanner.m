function [coef] = PathPlanner(state_x0, state_x1, d0, d1)

%starting position
%x =       state_x0(1);
%y =       state_x0(2);
%theta0 =  state_x0(3);
%theta1 =  state_x0(4)
%phi =     state_x0(5);

[yx_0,yxx_0,yxxx_0] = ComputeDerivatives(state_x0,d0,d1);
[yx_1,yxx_1,yxxx_1] = ComputeDerivatives(state_x1,d0,d1);

a0 = state_x0(2);
a1 = yx_0;
a2 = yxx_0/2;
a3 = yxxx_0/6;

x1 = state_x1(1);
A = [x1^7      x1^6       x1^5       x1^4;
     7*x1^6    6*x1^5     5*x1^4     4*x1^3;
     42*x1^5   30*x1^4    20*x1^3    12*x1^2;
     210*x1^4  120*x1^3   60*x1^2    24*x1^1;
     ];
 
B = [   state_x1(2) - ( a3*x1^3 + a2*x1^2 + a1*x1^1 + a0);
        yx_1 - ( a3*3*x1^2 + a2*2*x1^1 + a1);
        yxx_1 - (a3*6*x1^1 + a2*2);
        yxxx_1 - (6*a3);
    ];

coef = A\B;
a7 = coef(1);
a6 = coef(2);
a5 = coef(3);
a4 = coef(4);

coef = [a7 a6 a5 a4 a3 a2 a1 a0];




