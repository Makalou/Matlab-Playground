syms qw qx qy qz;
syms I11 I12 I13 I22 I23 I33;
syms L1 L2 L3;
syms h;
assume([qw qx qy qz],'real');
assume([I11 I12 I13 I22 I23 I33],'real');
assume([L1 L2 L3],'real');
assume(h,'real');

R = [1 - 2*qy*qy - 2*qz*qz, 2*qx*qy - 2*qw*qz, 2*qx*qz + 2*qw*qy;
     2*qx*qy + 2*qw*qz    , 1 - 2*qx*qx - 2*qz*qz, 2*qy*qz - 2*qw*qx;
     2*qx*qz - 2*qw*qy,     2*qy*qz + 2*qw*qx,   1 - 2*qx*qx - 2*qy*qy];

I = [I11, I12, I13; 
     I12, I22, I23;
     I13, I23, I33];

L = [L1; L2; L3];

w = [R * I * R' * L; 0];

wx = w(1);
wy = w(2);
wz = w(3);
ww = w(4);

dqdt = 1/2 * [qw * wx + qx * ww + qy * wz - qz * wy;
              qw * wy + qy * ww + qz * wx - qx * wz;
              qw * wz + qz * ww + qx * wy - qy * wx;
              qw * ww - qx * wx - qy * wy - qz * wz];

f = -[qx; qy; qz; qw] + dqdt * h;
%flen = sqrt(f(1) * f(1) + f(2) * f(2) + f(3) * f(3) + f(4) * f(4));
%f = f/flen;

disp(f);

Jq = jacobian(f,[qx, qy, qz, qw]);

disp(Jq);
%disp(ccode(Jq));


