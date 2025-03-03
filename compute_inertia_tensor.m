syms u v w;
syms x1 y1 z1 x2 y2 z2 x3 y3 z3 x4 y4 z4;
assume([u,v,x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4],'real');

% Define the vector function f(u, v)
f = [u*x1 + v*x2 + (1-u-v)*x3;
     u*y1 + v*y2 + (1-u-v)*y3;
     u*z1 + v*z2 + (1-u-v)*z3];

DfDu = diff(f,u);
DfDv = diff(f,v);
Ds = cross(DfDu,DfDv);
Ds = sqrt(Ds' * Ds);

A = [0 -f(3) f(2); f(3) 0 -f(1); -f(2) f(1) 0];
h = A'*A;
l = h * Ds;

I1 = int(int(l,u,[0 1-v]),v,[0 1]);
I1 = simplify(I1);
disp(I1);

% Define the vector function f(u, v)
f = [u*x1 + v*x2 + w*x3 + (1 - u - v - w) * x4;
     u*y1 + v*y2 + w*y3 + (1- u - v - w) * y4;
     u*z1 + v*z2 + w*z3 + (1 - u - v - w) * z4];

Jf = jacobian(f,[u,v,w]);

A = [0 -f(3) f(2); f(3) 0 -f(1); -f(2) f(1) 0];
h = A'*A;

l = h * Ds;

I2 = int(int(int(l,u,[0, 1 - v - w]),v,[0, 1 - w]),w,[0, 1]);
I2 = simplify(I2);
disp(I2);
