
c1 = [-1 -1 -1]; r1 = 2;
c2 = [1 1 1]; r2 = 2;

f1 = @(x,y,z) (x-c1(1))^2 + (y-c1(2))^2 + (z-c1(3))^2 - r1^2;
f2 = @(x,y,z) (x-c2(1))^2 + (y-c2(2))^2 + (z-c2(3))^2 - r2^2;

F = @(X) [f1(X(1),X(2),X(3)); f2(X(1),X(2),X(3))];
df1dx = @(X) [2*(X(1) - c1(1)) 2*(X(2) - c1(2)) 2*(X(3) - c1(3))];
df2dx = @(X) [2*(X(1) - c2(1)) 2*(X(2) - c2(2)) 2*(X(3) - c2(3))];
DFDX = @(X) [df1dx(X); df2dx(X)];

y = 8 * rand(3,1) - 4.0;
N = 10;

%{
  x - y + DF(x)'*lambda = 0
  F(x) = 0  
%}
X1 = zeros(3,N)';
X1(1,:) = y;
lambda = [0 0]';
for i = 1 : N - 1
    x = X1(i,:)';
    rhs = -[x - y + DFDX(x)'*lambda;
            F(x)];
    H = lambda(1)*diag([2 2 2]) + lambda(1) * diag([2 2 2]);
    JH = [eye(3) + H   DFDX(x)';
         DFDX(x)        zeros(2)];
    dz = JH \ rhs;
    z = [x; lambda] + dz;
    X1(i+1,:) = z(1:3);
    disp(F(x));
end

%{
  x - y + DF(x)'*lambda = 0
  F(x) = 0  
%}


h1 = fimplicit3(f1);
hold on;
h2 = fimplicit3(f2);
h1.FaceAlpha = 0.2;
h1.EdgeColor = 'none';
h2.FaceAlpha = 0.2;
h2.EdgeColor = 'none';
plot3(X1(:,1),X1(:,2),X1(:,3),'-*',Color='r');
axis equal


