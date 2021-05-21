% Jan Lao
% May 2021
% ValeroArm - BBDL
% Koby Python script reworked (fast_ffs.py)
% 2 Joint, 2 link planar system
% Arm plot representation
 
clc; clear all; close all;
tic

%% Initialize your link parameters
q = [0.7854,0.7854]; % Radians
l = [1,1]; % length of link
k = numel(q); % Number of joints
maxmotorforce = 1;
Rq = [-1,-1,1; -1,1,1]; % Optimal Moment arm matrix set
 
%% Limb Kinematics
% for planar limbs, we can remove alpha and define forward kinematic model
% to only consider displacements
Gq = [l(1)*cos(q(1))+l(2)*cos(q(1)+q(2)); 
    l(1)*sin(q(1))+l(2)*sin(q(1)+q(2))]; %endpoints
 
% All permutations of Jacobian
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

%% Convex Hull in FFS Space
 
% Creating convexhull of in field of W_T's data points
hull = convhull(W_T(:,1), W_T(:,2), 'simplify', true);

% Plot points of W
plot(W_T(:,1),W_T(:,2),'*')
title('Feasible Force Set vs. Largest Circle Function')
xlabel('Forces in X')
ylabel('Forces in Y')
xlim([-5 5])
ylim([-5 5])
hold on

% Plot Minkowski Sum
plot(W_T(hull,1),W_T(hull,2))
hold on

%% Plotting the arm

% Set shoulder base at (0,0)
x = 0;
y = 0;
q_n_k = 0;
scatter(x,y, 'filled')

for n = 1:k
    q_n_k = q_n_k + q(n);
    [x_k,y_k] = sph2cart(q_n_k, 0, 1);
    x(n+1) = x(n)+x_k;
    y(n+1) = y(n)+y_k;

    plot([x(n),x(n+1)],[y(n),y(n+1)])
    hold on
  
    scatter(x(n),y(n))
    scatter(x(n+1),y(n+1)) % end-effector location
end

% Making sure end-effector location agrees with limb kinematics
if Gq(1)~=x(end)
    Gq(1)=x(end)
    if Gq(2)~=y(end)
        Gq(2)=y(end)
    end
end

%% Circle within the polytope

% Circle Init
space = linspace(0,2*pi);
circ = [cos(space); sin(space)];
radiusbounds = [0.0, 1.0]; % Set: Initializing bound radius
inHull = 1;
tolerance = 0.33;
num_iterations = 0;
hul = W_T;
x_center = Gq(1); % Endpoints of limb = x and y center
y_center = Gq(2);
 
% Manipulating the circle in the polygon to gain set of optimized solution
while inHull == 1;
    if radiusbounds(2) - radiusbounds(1) < tolerance
        % End the loop and
        % "Return" biggest radiusbound and the points in the circle
        inHull = 0;
    else
        % Radiusbounds to be updated after every iteration
        r = (radiusbounds(2) - radiusbounds(1))/2;
        % Updating points based on new radius
        points_t = r*circ + [x_center; y_center];
        points = transpose(points_t);
        
        % Examining all the vertices of W_T if they are nonzero(1=true)
        was_too_small = all([W_T(hull,1),W_T(hull,2)]);
        
        if was_too_small >= 0 %if any/all vertice is nonnegative
            radiusbounds = [r,radiusbounds(2)];
            % Move circle to the right
        else
            radiusbounds = [radiusbounds(1),r];
            % Move circle to the left
        end
        num_iterations = num_iterations+1;
        fprintf('Iteration #: %d \n', num_iterations);
    end
end

hold on
plot(points(:,1),points(:,2))
hold off
fprintf('Total # of iterations: %d \n', num_iterations);
toc