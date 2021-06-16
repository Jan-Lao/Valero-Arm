% Jan Lao
% June 2021
% ValeroLab - ValeroArm
% Perpendicular Line Unit Test
close all; clear all; clc;

%% Note to the Reader:
% Criterion for on-polygon as well as eliminating potentially redundant
% hull points is eliminate by the convhull() and inpolygon() commands

%% Trying Different Vertices
% Comment the 2 sets of vertex_x and vertex_y to test out the uncommented.
% Sanity Check (With redundant hull point)
vertex_x = [0; 8; 10; 10; 1; 0];
vertex_y = [0; 1.5; 10; 10; 7; 0];

% Random vertices (Note: as long as its perpendicular it is good!)
rng1 = randi(10, 3, 1);
rng2 = randi(10, 3, 1);
%vertex_x = [0; rng1; 0];
%vertex_y = [0; rng2; 0];

% Close to original code parameters (W_T(hull,1) and W_T(hull, 2))
%vertex_x = [0; 2; 3; 0; -2; -3; 0];
%vertex_y = [-5.7; -2.2; 4.4; 0; -3.4; -10; -5.7];

%% Setting up the rest of the points
% Test out different numbers here
center =  [0,-3]; %or [5, 5]; 
vertices = [vertex_x, vertex_y];
vertices(numel(vertex_x)+1,:) = vertices(1,:); %Wrap around for iterative

%% Single Point/Edge Analysis
% vector = head point - tail point
vertex_A = [vertices(1,1), vertices(1,2)];
vertex_B = [vertices(2,1), vertices(2,2)];

vector_A_to_B = vertex_B - vertex_A;
vector_A_to_Center = center - vertex_A;

% Projection of x onto v
mag_sq1 = vector_A_to_B(1)^2 + vector_A_to_B(2)^2; %intermediate calculation
proj_xv1 = ((dot(vector_A_to_B,vector_A_to_Center))/(mag_sq1))*vector_A_to_B; %something wrong here
w2 = vector_A_to_Center - proj_xv1; %Orthagonal component of x1
w2_mag = sqrt(w2(1)^2+w2(2)^2); %Potential radius
closest_point = [proj_xv1(1) + vertex_A(1), proj_xv1(2) + vertex_A(2)];

figure(1)
scatter(center(1), center(2), 'filled') %center
hold on
plot(vertex_x, vertex_y, 'k') %"hull"
hold on
scatter(vertex_A(1), vertex_A(2), '*') %First vertex
hold on
scatter(closest_point(1),closest_point(2)) %Point P
hold on
%plot([vertex_A(1) vertex_B(1)],[vertex_A(2) vertex_B(2)]) %vector v1
%hold on
plot([vertex_A(1) closest_point(1)],[vertex_A(2) closest_point(2)], '--r') %projection vector
hold on
plot([vertex_A(1) center(1)],[vertex_A(2) center(2)], 'g') %vector x1
hold on
plot([center(1) closest_point(1)],[center(2) closest_point(2)], ':b') %w2
hold off

axis square
xlim([-12 12])
ylim([-12 12])

%% Iterative Analysis
figure (2)
for n = 1:numel(vertex_x)
    vertex_1(n,:) = [vertices(n,1), vertices(n,2)];
    vertex_2(n,:) = [vertices(n+1,1), vertices(n+1,2)];
    
    vector_v(n,:) = vertex_2(n,:) - vertex_1(n,:);
    vector_x(n,:) = center - vertex_1(n,:);
    
    mag_sq =  vector_v(n,1)^2 + vector_v(n,2)^2; %intermediate calculation
    proj_xv(n,:) = ((dot(vector_v(n,:),vector_x(n,:)))/(mag_sq))*vector_v(n,:);
    proj_xv(isnan(proj_xv))= 0;
    
    D(n,:) = vector_x(n,:) - proj_xv(n,:);
    D_mag(n,:) = sqrt(D(n,1)^2 + D(n,2)^2);
    
    p(n,:) = [proj_xv(n,1) + vertex_1(n,1), proj_xv(n,2) + vertex_1(n,2)];
    
    scatter(vertex_1(n,1), vertex_1(n,2), '*') %First vertices
    hold on
    scatter(p(n,1),p(n,2)) %Point P
    hold on
    %plot([vertex_1(n,1) vertex_2(n,1)],[vertex_1(n,2) vertex_2(n,2)], 'g') %vector v
    %hold on
    plot([vertex_1(n,1) p(n,1)],[vertex_1(n,2) p(n,2)], '--r') %projection vector
    hold on
    plot([center(1) p(n,1)],[center(2) p(n,2)], ':b') %D
    hold on
end

D_info = [D, D_mag];
hold on
scatter(center(1), center(2), 'filled') %center
hold on
plot(vertex_x, vertex_y, 'k') %"hull"
hold off

axis square
xlim([-12 12])
ylim([-12 12])