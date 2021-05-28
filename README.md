# Valero-Arm
Please refer to https://github.com/nkreid/ValeroArm for the tendon route optimization problem statement as well as Chapter 7 of Dr. Francisco Valero-Cuevas' BME 504 textbook: Fundamentals of Neuromechanics.

# Initial Task and Initial Hurdles
The challenge here was to take Koby's already working tendon route optimization code and to translate it completely into Matlab for the lab convenience and standardization purposes here at the ValeroLab. In order to do so, there were some issues with the iterative, largest circle criterion. As the while loop would take hours and even days to find the optimal solution within the feasible force set polytope: the largest circle possible without crossing the polytope (convex hull) boundaries.

# Additions to Koby's Tendon Route Optimization Code
Along with the intial parameters, limb mechanics, limb kinematics, and the creation of the Minkowski sum, and subsequently the convex hull, as derived from Dr. Valero-Cuevas' textbook and Koby's Tendon Route Optimization code, a graphical representation of the arm must be created in order to visualize the arm in relation to the shoulder joint @ (0,0) and the end-effector, the center of the optimization circle.

# Overcoming the Inherent Challenges of the Optimization Criteria
Instead of providing a while loop, I've opted in for the perpendicular line test. Such a test, takes draws a perpendicular line (the shortest distance) from the center of the circle to the edge of each polytope boundary. Then, the each line from the center the edge is evaluated, if the line crosses the polytope edge, due to the point of intersection being outside the polytope, that perpendicular line and intersection point is thrown out.

The line with the shortest distance from the center to the edge of the polytope that is still within the polytope is the radius of the largest circle possible within the polytope.
