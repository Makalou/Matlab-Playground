R = eye(3);

x10 = [1 0 0]';
x20 = [0 1 0]';
x30 = [0 0 1]';

m1 = 500; m2 = 200; m3 = 1;

center_mass = (m1 * x10 + m2 * x20 + m3*x30)/(m1 + m2 + m3);
x10 = x10 - center_mass;
x20 = x20 - center_mass;
x30 = x30 - center_mass;

% Compute the inertia tensor components

Ixx = m1*(x10(2)^2 + x10(3)^2) + m2*(x20(2)^2 + x20(3)^2) + m3*(x30(2)^2 + x30(3)^2);
Iyy = m1*(x10(1)^2 + x10(3)^2) + m2*(x20(1)^2 + x20(3)^2) + m3*(x30(1)^2 + x30(3)^2);
Izz = m1*(x10(1)^2 + x10(2)^2) + m2*(x20(1)^2 + x20(2)^2) + m3*(x30(1)^2 + x30(2)^2);

Ixy = -m1*x10(1)*x10(2) - m2*x20(1)*x20(2) - m3*x30(1)*x30(2);
Ixz = -m1*x10(1)*x10(3) - m2*x20(1)*x20(3) - m3*x30(1)*x30(3);
Iyz = -m1*x10(2)*x10(3) - m2*x20(2)*x20(3) - m3*x30(2)*x30(3);

% Construct the inertia tensor matrix
J = [Ixx, Ixy, Ixz;
     Ixy, Iyy, Iyz;
     Ixz, Iyz, Izz];

J_inv = inv(J);

torque = [0 0 0]';

dt = 0.01;
T = 50;
N = (T/dt);

X1 = zeros(3,N)';
X2 = zeros(3,N)';
X3 = zeros(3,N)';
L = zeros(3,N)';
W = zeros(3,N)';
match_errors = zeros(1,N)';

Y1 = zeros(3,N)';
Y2 = zeros(3,N)';
Y3 = zeros(3,N)';

y1_old = x10;
y2_old = x20;
y3_old = x30;

[V,D] = eig(R*J*R');

%w0 = 100 * V(:,1);
w0 =  [20 * rand(1)-10 20 * rand(1) - 10 20 * rand(1) - 10]';

w = w0;
w2 = w0;

for i = 1 : N
    L1 = R*J*R'*w;
    [w R] = rk4_step(w, R,J,J_inv,torque,dt);
    %w = w*0.999;
    fprintf("L : [%f %f %f]",L1(1),L1(2),L1(3));
    disp(norm(R*R'-eye(3),'fro'));

    X1(i,:) = R*x10;
    X2(i,:) = R*x20;
    X3(i,:) = R*x30;

    W(i,:) = w;
    L(i,:) = L1;

    mc_old = (m1 * y1_old + m2 * y2_old + m3 * y3_old)/(m1 + m2 + m3);
    y1_old_center = y1_old - mc_old;
    y2_old_center = y2_old - mc_old;
    y3_old_center = y3_old - mc_old;

    v1 = cross(w2,y1_old_center); % for now angular vel only
    v2 = cross(w2,y2_old_center);
    v3 = cross(w2,y3_old_center);

    y1_in = y1_old + v1 * dt;
    y2_in = y2_old + v2 * dt;
    y3_in = y3_old + v3 * dt;

    mc_in = (m1 * y1_in + m2 * y2_in + m3 * y3_in)/(m1 + m2 + m3);
   
    y1_in_center = y1_in - mc_in;
    y2_in_center = y2_in - mc_in;
    y3_in_center = y3_in - mc_in;

    % Compute covariance matrix
    % I don't know why it better to put mass here
    % I don't even sure it works for all scenarios :(
    H = (m1*y1_in_center * x10' + m2*y2_in_center * x20' + m3*y3_in_center * x30');
    [U,~,V] = svd(H);

    R_opt = U * V';

   if det(R_opt) < 0
        U(:,end) = -U(:,end);
        R_opt = U * V';
    end 

    disp(norm(R_opt * R_opt' - eye(3) ,'fro'));

    m_err = norm(R_opt * x10 - y1_in_center,2) + norm(R_opt * x20 - y2_in_center,2) + norm(R_opt * x30 - y3_in_center,2);

    match_errors(i) = m_err;

    y1_new = R_opt * x10 + mc_in;
    y2_new = R_opt * x20 + mc_in;
    y3_new = R_opt * x30 + mc_in;

    w2 = implicit_solve_w(R_opt,J,J_inv,w2,torque,dt);

    Y1(i,:) = y1_new;
    Y2(i,:) = y2_new;
    Y3(i,:) = y3_new;

    y1_old = y1_new;
    y2_old = y2_new;
    y3_old = y3_new;
end
figure;
L_err = L - L(1,:);
plot3(L_err(:,1),L_err(:,2),L_err(:,3));

fig2 = figure;
rotate3d on
axis equal;
axis([-2 2 -2 2 -2 2]);
plot3(Y1(:,1),Y1(:,2),Y1(:,3));
hold on;
plot3(Y2(:,1),Y2(:,2),Y2(:,3));
hold on;
plot3(Y3(:,1),Y3(:,2),Y3(:,3));

fig3 = figure;
rotate3d on
axis equal;
axis([-2 2 -2 2 -2 2]);
plot3(X1(:,1),X1(:,2),X1(:,3));
hold on;
plot3(X2(:,1),X2(:,2),X2(:,3));
hold on;
plot3(X3(:,1),X3(:,2),X3(:,3));

hold on;

% Initialize plot objects for the triangle edges
h1 = plot3([0 0], [0 0], [0 0], 'r-', 'LineWidth', 2); % Edge x1-x2
h2 = plot3([0 0], [0 0], [0 0], 'r-', 'LineWidth', 2); % Edge x2-x3
h3 = plot3([0 0], [0 0], [0 0], 'r-', 'LineWidth', 2); % Edge x3-x1

g1 = plot3([0 0], [0 0], [0 0], 'b-', 'LineWidth', 2); % Edge x1-x2
g2 = plot3([0 0], [0 0], [0 0], 'b-', 'LineWidth', 2); % Edge x2-x3
g3 = plot3([0 0], [0 0], [0 0], 'b-', 'LineWidth', 2); % Edge x3-x1

h4 = plot3([0 0], [0 0], [0 0], 'r-', 'LineWidth', 2); % rotation axis
h5 = plot3([0 0], [0 0], [0 0], 'b-', 'LineWidth', 2); % rotation axis

h6 = plot3([0 0], [0 0], [0 0], 'g-', 'LineWidth', 2); % angular momentumn axis

% Animation loop
for i = 1:N
    title('frame',i);
    if ~isvalid(fig3)
        break;
    end

    h1.XData = [X1(i, 1), X2(i, 1)];
    h1.YData = [X1(i, 2), X2(i, 2)];
    h1.ZData = [X1(i, 3), X2(i, 3)];
        
    h2.XData = [X2(i, 1), X3(i, 1)];
    h2.YData = [X2(i, 2), X3(i, 2)];
    h2.ZData = [X2(i, 3), X3(i, 3)];
        
    h3.XData = [X3(i, 1), X1(i, 1)];
    h3.YData = [X3(i, 2), X1(i, 2)];
    h3.ZData = [X3(i, 3), X1(i, 3)];

    h4.XData = [-abs(W(i, 1)), 0];
    h4.YData = [-abs(W(i, 2)), 0];
    h4.ZData = [-abs(W(i, 3)), 0];

    h5.XData = [0, abs(W(i, 1))];
    h5.YData = [0, abs(W(i, 2))];
    h5.ZData = [0, abs(W(i, 3))];

    h6.XData = [0, L(i, 1)];
    h6.YData = [0, L(i, 2)];
    h6.ZData = [0, L(i, 3)];

    
    g1.XData = [Y1(i, 1), Y2(i, 1)];
    g1.YData = [Y1(i, 2), Y2(i, 2)];
    g1.ZData = [Y1(i, 3), Y2(i, 3)];
        
    g2.XData = [Y2(i, 1), Y3(i, 1)];
    g2.YData = [Y2(i, 2), Y3(i, 2)];
    g2.ZData = [Y2(i, 3), Y3(i, 3)];
        
    g3.XData = [Y3(i, 1), Y1(i, 1)];
    g3.YData = [Y3(i, 2), Y1(i, 2)];
    g3.ZData = [Y3(i, 3), Y1(i, 3)];
    
    
    axis([-2 2 -2 2 -2 2]);
    drawnow;
    % Pause to create animation effect
    pause(dt);
end

function A = skew(v)
    A = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
end

function [w_new, R_new] = euler_step(w_old, R_old,J,J_inv,torque,dt)
    L = R_old*J*R_old'*w_old;
    alpha = R_old*J_inv*R_old'*(torque - cross(w_old,L));
    w_new = w_old + alpha * dt;
    R_new = expm(skew(w_new)*dt) * R_old;
end

function R = mymap(A)
    R = expm(A);
end

function [w_new, R_new] = rk4_step(w_old, R_old,J,J_inv,torque,dt)
    L = R_old*J*R_old'*w_old;
    alpha1 = R_old*J_inv*R_old'*(torque - cross(w_old,L));
    A1 = skew(w_old);

    w_in = w_old + dt/2 * alpha1;
    R_in = mymap(A1*dt/2) * R_old;

    L = R_in*J*R_in'*w_in;
    alpha2 = R_in*J_inv*R_in'*(torque - cross(w_in,L));
    A2 = skew(w_in);

    w_in = w_old + dt/2 * alpha2;
    R_in = mymap(A2*dt/2) * R_old;

    L = R_in*J*R_in'*w_in;
    alpha3 = R_in*J_inv*R_in'*(torque - cross(w_in,L));
    A3 = skew(w_in);

    w_in = w_old + dt * alpha3;
    R_in = mymap(A3*dt) * R_old;

    L = R_in*J*R_in'*w_in;
    alpha4 = R_in*J_inv*R_in'*(torque - cross(w_in,L));
    A4 = skew(w_in);

    w_new = w_old + (alpha1 + 2*alpha2 + 2*alpha3 + alpha4)*dt/6;

    R_new = mymap((A1 + 2*A2 + 2*A3 + A4)*dt/6)*R_old;
end

function w_new = explicit_solve_w(R,J,J_inv,w_old,torque,dt)
    alpha = R*J_inv*R'*(torque - cross(w_old,R*J*R'*w_old));
    w_new = w_old + alpha * dt;
end

function w_new = implicit_solve_w(R,J,J_inv,w_old,torque,dt)
    w_new = w_old;
    A = R*J*R';
    for i = 1 : 32
        f = A*(w_new - w_old) + cross(w_new,A*w_new) * dt - torque * dt;
        dfdw = A - skew(A*w_new) * dt;
        rhs = -f;
        dw = dfdw \ rhs;
        w_new = w_new + dw;
    end
end

