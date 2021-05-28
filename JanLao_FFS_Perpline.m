% Jan Lao
% May 2021
% ValeroLab - ValeroArm
% 2 Joint, 2 link planar, 3 muscle system
% No reiterating circles: Perpendicular line evaluation
clc; clear all; close all;
tic

%% Initialize your link parameters
q = [0.7854,0.7854]; % Radians
l = [1,1]; % length of link
num_joints = numel(q); % k
num_muscles = num_joints+1;
maxmotorforce = 1;
Rq = [-2,-3,1; -3,1,2]; % Optimal Moment arm matrix set
 
%% Limb Kinematics
% Endpoints
Gq = [l(1)*cos(q(1))+l(2)*cos(q(1)+q(2)); 
    l(1)*sin(q(1))+l(2)*sin(q(1)+q(2))];
Gq(3) = 0; % 2D for now: let z_center = 0
 
% Permutations of Jacobian
J = [-l(2)*sin(q(1)+q(2))-l(1)*sin(q(1)), -l(2)*sin(q(1)+q(2)); 
    l(2)*cos(q(1)+q(2))+l(1)*cos(q(1)), l(2)*cos(q(1)+q(2))];
J_inv = inv(J);
J_invT = transpose(J_inv);
 
%% Limb Mechanics
% f0(q,qdot)
f0diag = [maxmotorforce, maxmotorforce, maxmotorforce];
f0 = diag(f0diag);
 
% H Matrix
H = J_invT*Rq*f0;
 
% A possibilities of muscle activation - neural activation
a_poss = [1,1,1; 1,0,0; 1,0,1; 1,1,0; 0,1,1; 0,1,0; 0,0,1; 0,0,0];
a_T = transpose(a_poss);
 
% Wrench - Minkowski Sum
W = zeros(size(H,1),size(a_T,2)); % Initialize matrix dimensions
for n = 1:size(W,2) % Getting range of i from 1 to number of columns in W
   W(:,n)  = H*a_T(:,n);
   % update with H dotted with a-transpose at every i-th column of W 
end
W_T = transpose(W);

%% Plotting the arm
% Set shoulder base at (0,0) and creating first figure
x = 0;
y = 0;
q_n_k = 0;
figure(1)
tile = tiledlayout(1,1);
ax1 = axes(tile);
scatter(ax1,x,y, 'filled')
xlabel('Forces in X')
ylabel('Forces in Y')
xlim([-10 10])
ylim([-10 10])
axis square
hold on

ax2 = axes(tile);
for n = 1:num_joints
    q_n_k = q_n_k + q(n);
    [x_k,y_k] = sph2cart(q_n_k, 0, 1);
    x(n+1) = x(n)+x_k;
    y(n+1) = y(n)+y_k;

    plot(ax2, [x(n),x(n+1)],[y(n),y(n+1)])
    hold on
  
    scatter(x(n),y(n))
    scatter(x(n+1),y(n+1)) % end-effector location
end

ax2.XAxisLocation = 'top';
ax2.YAxisLocation = 'right';
xlabel('X Arm Position')
ylabel('Y Arm Position')
xlim([-10 10])
ylim([-10 10])
axis square
title('Feasible Force Set vs. Circle Center: Dist. Relation to End-Effector')

%% Convex Hull in FFS Space
% Creating convexhull of in field of W_T's data points
hull = convhull(W_T(:,1) + Gq(1), W_T(:,2) + Gq(2), 'simplify', true);

% Plot points of W and Minkowski Sum (polytope)
scatter(W_T(:,1) + Gq(1), W_T(:,2) + Gq(2),'*')
hold on
plot(W_T(hull,1) + Gq(1), W_T(hull,2) + Gq(2))
hold on

%% Creating All Potential Radii from Center to Edge of Polytope
% Init vertex arrays as set of points for projection function
vertex_x = W_T(hull,1);
vertex_y = W_T(hull,2);
vertex_z = zeros(size(W_T(hull))); % W_T(hull,3);
vector_v = []; vector_x = [];
proj_xv = zeros(numel(W_T(hull)),num_muscles);
D = [];
p_x = []; p_y = [];

% Last vertex in polygon gets compared to the first one (wrap around)
vertex_x(numel(W_T(hull))+1) = vertex_x(1);
vertex_y(numel(W_T(hull))+1) = vertex_y(1);
vertex_z(numel(W_T(hull))+1) = vertex_z(1);
proj_xv(numel(W_T(hull))+1) = proj_xv(1);

% Circle center init and graphing it
center = [Gq(1),Gq(2),Gq(3)];
figure(2)
plot(center(1),center(2),'*')
title('Distances From End-Effector to Polytope Edges')
xlabel('Forces in X')
ylabel('Forces in Y')
xlim([-10 10])
ylim([-10 10])
axis square
hold on

% Graphing hull
plot(W_T(hull,1) + Gq(1),W_T(hull,2) + Gq(2))

% Perpendicular line (D) and point (P) function
for n = 1:numel(W_T(hull))
    % Vertices get overwritten per iteration (sides of polytope)
    vertex1 = [vertex_x(n), vertex_y(n), vertex_z(n)];
    vertex2 = [vertex_x(n+1), vertex_y(n+1), vertex_z(n+1)];

    % Compute distance from center to line segment (polytope edgeline)
    vector_v(n,:) = vertex2 - vertex1;
    vector_x(n,:) = center - vertex1;
    proj_xv(n,:) = (dot(vector_v(n,:),vector_x(n,:))/(dot(vector_v(n,:),vector_v(n,:))))*vector_v(n,:);
    proj_xv(isnan(proj_xv)) = 0; % Setting any NaN to 0
    D(n,:) = vector_x(n,:) - proj_xv(n,:); % distance vectors from center
    D_mag(n,:) = sqrt(sum(D(n,:).^2)); % shortest dist from point to line

    % Finding the point of perpendendicular line and line segment span
    p_x = proj_xv(:,1) + vertex1(1);
    p_y = proj_xv(:,2) + vertex1(2); 
    
    % Graph each iteration of point P and vector D
    scatter(p_x(n), p_y(n), 'filled') % Point P
    hold on
    plot([Gq(1) p_x(n)],[Gq(2) p_y(n)]) % Vector D - potential radii
    hold on
end

%% Eliminate "bad" solutions
% Evaluating if P is outside of the polygon
inHull = inpolygon(p_x,p_y,W_T(hull,1),W_T(hull,2)); % returns boolean
perp_points = [inHull, p_x, p_y, D_mag]; 
% make array: boolean status, point (P), and associated distance (D_mag)

% If outside the polygon eliminate corresponding D_mag and P
perp_points = perp_points(~(inHull==0),:);

%% Finding the Largest Circle
% The smallest D_mag = the largest radii of the circle in the polytope.
% Also finding the min value in the last column as well as its row index
[radius, r_index] = min(perp_points(:,4)); 
space = linspace(0,2*pi);
circle = radius*[cos(space); sin(space)] + [Gq(1); Gq(2)];

% Graphing only the largest radii and the circle in the convhull
figure(3)
plot(W_T(hull,1) + Gq(1), W_T(hull,2) + Gq(2)) % Graph hull
hold on
plot(center(1),center(2),'*') % Graph center of circle
hold on
% Plotting the correct p_x and p_y
scatter(perp_points(r_index,2), perp_points(r_index,3), 'filled')
hold on
% Plotting largest radius
plot([Gq(1) perp_points(r_index,2)],[Gq(2) perp_points(r_index,3)])
% Plotting largest circle
plot(circle(1,:),circle(2,:))
title('Polytope vs. Circle and Its Largest Radii')
xlabel('Forces in X')
ylabel('Forces in Y')
xlim([-10 10])
ylim([-10 10])
axis square
hold off
toc