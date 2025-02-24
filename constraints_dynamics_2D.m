N = 10;  % Number of points
M = 10;  % Number of edges
X0 = 2 * rand(N,2) - 1;
E = zeros(N-1, 2);

perm = randperm(N);
for i = 2:N
    E(i-1, :) = [perm(i-1), perm(i)];
end

fixed_point = [perm(1)];

while size(E, 1) < M
    node1 = randi(N);
    node2 = randi(N);
    if node1 == node2
        continue; 
    end
    if any((E(:,1)==node1 & E(:,2)==node2) | (E(:,1)==node2 & E(:,2)==node1))
        continue; 
    end
    E(end+1, :) = [node1, node2];
end

% Construct constraints functions
Distance = zeros(size(E, 1),1);
for i = 1:size(E, 1)
    n1 = E(i, 1);
    n2 = E(i, 2);
    Distance(i) =  sqrt((X0(n1,1) - X0(n2,1))^2 + (X0(n1,2) - X0(n2,2))^2);
end

C = cell(size(E, 1),1);
for i = 1:size(E, 1)
    n1 = E(i, 1);
    n2 = E(i, 2);
    C{i} = @(X) (X(n1,1) - X(n2,1))^2 + (X(n1,2) - X(n2,2))^2 - Distance(i)^2;
end

DCDX = cell(size(E, 1),1);
for i = 1:size(E, 1)
    n1 = E(i, 1);
    n2 = E(i, 2);

    n1 = n1 * 2 - 1;
    m1 = n1 + 1;
    n2 = n2 * 2 - 1;
    m2 = n2 + 1;
    % X is flattened 2Nx1 array
    if n1 < n2
        DCDX{i} = @(X) [zeros(1, n1-1), 2*(X(n1) - X(n2)), 2*(X(m1) - X(m2)), zeros(1, n2-n1-2), -2*(X(n1) - X(n2)), -2*(X(m1) - X(m2)), zeros(1, length(X)-n2-1)];
    else
        DCDX{i} = @(X) [zeros(1, n2-1), -2*(X(n1) - X(n2)), -2*(X(m1) - X(m2)), zeros(1, n1-n2-2), 2*(X(n1) - X(n2)), 2*(X(m1) - X(m2)), zeros(1, length(X)-n1-1)];
    end
end

Nt = 1000;
dt = 0.01;

X = zeros(N,2,Nt);
X(:,:,1) = X0;

V = zeros(N, 2);
g = [0 -9.8]; % set value to large then you can see the difference between 2 methods
% newton solver tend to suffer from overshooting problem.

% Simulation loop
for t = 1 : Nt - 1
    V = V + g * dt;
    X_in = X(:,:,t) + V * dt;
    X(:,:,t+1) = gauss_newton_solver(X_in,X0,fixed_point,C,DCDX);
    V = (X(:,:,t+1) - X(:,:,t))/dt;
end

f1 = figure(1);
h1 = scatter(X0(:,1),X0(:,2),'filled','b');
axis([-5 5 -5 5]);
hold on;
scatter(X0(fixed_point,1),X0(fixed_point,2),'r');
hold on;
h_edges = gobjects(size(E, 1), 1);
for i = 1:size(E, 1)
    n1 = E(i, 1);
    n2 = E(i, 2);
    h_edges(i) = plot([X0(n1,1), X0(n2,1)], [X0(n1,2), X0(n2,2)], 'g-', 'LineWidth', 1.0);
    hold on;
end

%input("wait for any input");
for i = 1 : Nt
    if ~isvalid(f1)
        break;
    end
    title("frame",i);
    X_now = X(:,:,i);
    h1.XData = X_now(:,1);
    h1.YData = X_now(:,2);
    for j = 1:size(E, 1)
        n1 = E(j, 1);
        n2 = E(j, 2);
        h_edges(j).XData = [X_now(n1,1), X_now(n2,1)];
        h_edges(j).YData = [X_now(n1,2), X_now(n2,2)];
    end
    pause(dt);
end

function X = gauss_seidel_solver(X,X0,fixed_point,C,DCDX)
    N = size(X,1);
    for i = 1 : 5 % number of constraints solve step
        X(fixed_point,:) = X0(fixed_point,:);
        for k = 1 : size(C,1)
            Ck = C{k}(X);
            X_flatten = reshape(X', 1, []);
            DCDXk = DCDX{k}(X_flatten);
            X_flatten = X_flatten - Ck / (DCDXk * DCDXk') * DCDXk; 
            X = reshape(X_flatten, 2, N)';
        end
    end
    X(fixed_point,:) = X0(fixed_point,:);
end

function X = gauss_newton_solver(X,X0,fixed_point,C,DCDX)
    N = size(X,1);
    gamma = 1e3;
    Y_flatten = reshape(X', 1, []);
    for i = 1 : 5
        X(fixed_point,:) = X0(fixed_point,:);
        Cx = cellfun(@(Ck) Ck(X),C);
        X_flatten = reshape(X', 1, []);
        DCDx = cellfun(@(DCDXk) DCDXk(X_flatten),DCDX,UniformOutput=false);
        DCDx = cell2mat(DCDx);
        M = eye(2*N) + gamma * (DCDx' * DCDx);
        rhs = -(X_flatten' - Y_flatten') - gamma *DCDx'*Cx;
        dX = M \ rhs;
        X_flatten = X_flatten + dX';
        X = reshape(X_flatten, 2, N)';
    end
    X(fixed_point,:) = X0(fixed_point,:);
end