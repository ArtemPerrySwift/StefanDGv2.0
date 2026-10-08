// Gmsh project created on Mon Apr 28 17:30:58 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 2.0};
//+
Point(2) = {1, 0, 0, 2.0};
//+
Point(3) = {0, 1, 0, 2.0};
//+
Point(4) = {0, 0, 1, 2.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {3, 2};
//+
Line(4) = {1, 4};
//+
Line(5) = {4, 2};
//+
Line(6) = {4, 3};
//+
Curve Loop(1) = {5, 1, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {4, 6, -2};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {5, -3, -6};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, 3, 1};
//+
Plane Surface(4) = {4};
//+
Surface Loop(1) = {3, 1, 4, 2};
//+
Volume(1) = {1};
//+
Physical Surface("D", 7) = {3, 2, 1};
//+
Physical Volume("M:L; 1.5;1.3;1.7;", 8) = {1};
