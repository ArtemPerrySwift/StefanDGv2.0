// Gmsh project created on Mon Apr 28 18:38:01 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.5};
//+
Point(2) = {1, 0, 0, 2.0};
//+
Point(3) = {0, 1, 0, 2.0};
//+
Point(4) = {0, 0, 1, 2.0};
//+
Point(5) = {-1, 0, 0, 2.0};
//+
Line(1) = {5, 1};
//+
Line(2) = {1, 4};
//+
Line(3) = {4, 5};
//+
Line(4) = {2, 1};
//+
Line(5) = {4, 2};
//+
Line(6) = {3, 2};
//+
Line(7) = {1, 3};
//+
Line(8) = {3, 5};
//+
Line(9) = {4, 3};
//+
Curve Loop(1) = {4, -1, -8, 6};
//+
Surface(1) = {1};
//+
Curve Loop(3) = {9, 6, -5};
//+
Surface(2) = {3};
//+
Curve Loop(5) = {8, -3, 9};
//+
Surface(3) = {5};
//+
Curve Loop(7) = {4, -1, -3, 5};
//+
Surface(4) = {7};
//+
Surface Loop(1) = {2, 4, 3, 1};
//+
Volume(1) = {1};
//+
Physical Surface("D", 10) = {2, 3, 4};
//+
Physical Volume("M:L; 1.0;0.0;0.0;", 11) = {1};
