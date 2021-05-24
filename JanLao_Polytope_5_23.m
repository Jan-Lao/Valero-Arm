% Jan Lao
% May 2021
% ValeroLab - ValeroArm
% % 2 Joint, 2 link planar system
% Creating Polytope using Francisco's textbook matlab script
clc; clear all; close all;
tic

%% Initialize your link parameters
q = [0.7854,0.7854]; % Radians
l = [1,1]; % length of link
num_joints = numel(q); % k
num_muscles = num_joints+1;
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
%H= rand(2,3)*2-1;

%% Polytope

vertices = [0 0;1 0; 0 1 ; 1 1];
new = [0 1];
num_muscles = num_muscles-1;
for n=2:num_muscles
    temp2 = [];
    for i=1:length(vertices')
        row =  vertices(i,:);
        temp1 = [0 row;
             1 row];
                 temp2 = [temp2;temp1];
    end
    vertices = temp2;
end
vertices = sortrows(vertices);
count = length(vertices');

H2= H(1:2,1:3);

% 2D case
% multiply each vertex by the matrix H
Y=[];
for i=1:count
    Y= [Y;(H2*vertices(i,:)')'];
end

% Find the convex hull
K = convhull(Y);

% plot all points and the convex hull.
figure(2)
hold on
plot(Y(K,1),Y(K,2),'r-')
toc