
c1 = [-1 -1 -1]; r1 = 2;
c2 = [1 1 1]; r2 = 2;
c3 = [0 -1 0]; r3 = 1;

f1 = @(x,y,z) (x-c1(1))^2 + (y-c1(2))^2 + (z-c1(3))^2 - r1^2;
f2 = @(x,y,z) (x-c2(1))^2 + (y-c2(2))^2 + (z-c2(3))^2 - r2^2;
f3 = @(x,y,z) (x-c3(1))^2 + (y-c3(2))^2 + (z-c3(3))^2 - r3^2;

F = @(X) [f1(X(1),X(2),X(3)); f2(X(1),X(2),X(3)); f3(X(1),X(2),X(3))];
df1dx = @(X) [2*(X(1) - c1(1)) 2*(X(2) - c1(2)) 2*(X(3) - c1(3))];
df2dx = @(X) [2*(X(1) - c2(1)) 2*(X(2) - c2(2)) 2*(X(3) - c2(3))];
df3dx = @(X) [2*(X(1) - c3(1)) 2*(X(2) - c3(2)) 2*(X(3) - c3(3))];
DFDX = @(X) [df1dx(X); df2dx(X); df3dx(X)];

y = 8 * rand(3,1) - 4.0;
N = 16;

%{
  x - y + DF(x)'*lambda = 0
  F(x) = 0  
%}
X1 = zeros(3,N)';
X1(1,:) = y;
lambda = [0 0 0]';
for i = 1 : N - 1
    x = X1(i,:)';
    rhs = -[x - y + DFDX(x)'*lambda;
            F(x)];
    H = zeros(3);%lambda(1)*diag([2 2 2]) + lambda(1) * diag([2 2 2]);
    JH = [eye(3) + H   DFDX(x)';
         DFDX(x)        zeros(3)];
    dz = JH \ rhs;
    z = [x; lambda] + dz;
    X1(i+1,:) = z(1:3);
    disp(F(x));
end

%{
    
%}
gamma = 1e3;
X2 = zeros(3,N)';
X2(1,:) = y;
for i = 1 : N - 1
    x = X2(i,:)';
    M = eye(3) + gamma*DFDX(x)'*DFDX(x);
    rhs = -(x - y) - gamma * DFDX(x)' * F(x);
    dx = M \ rhs;
    X2(i+1,:) = x + dx;
end

X3 = zeros(3,N)';
X3(1,:) = y;
for i = 1 : N -1
    x = X3(i,:)';
    x = x - f1(x(1),x(2),x(3))/(df1dx(x)*df1dx(x)')*df1dx(x)';
    x = x - f2(x(1),x(2),x(3))/(df2dx(x)*df2dx(x)')*df2dx(x)';
    x = x - f3(x(1),x(2),x(3))/(df3dx(x)*df3dx(x)')*df3dx(x)';
    X3(i+1,:) = x;
end

h1 = fimplicit3(f1);
hold on;
h2 = fimplicit3(f2);
hold on;
h3 = fimplicit3(f3);
h1.FaceAlpha = 0.2;
h1.EdgeColor = 'none';
h2.FaceAlpha = 0.2;
h2.EdgeColor = 'none';
h3.FaceAlpha = 0.2;
h3.EdgeColor = 'none';
plot3(X1(:,1),X1(:,2),X1(:,3),'*-',Color='r',LineWidth=1);
hold on;
plot3(X2(:,1),X2(:,2),X2(:,3),'*-',Color='b',LineWidth=1);
hold on;
plot3(X3(:,1),X3(:,2),X3(:,3),'*-',Color='g',LineWidth=1);
axis equal


