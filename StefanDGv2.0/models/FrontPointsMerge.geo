// Gmsh project created on Tue Oct 20 21:37:24 2026
SetFactory("OpenCASCADE");
//+
Line(9) = {4, 12};
//+
Line(10) = {12, 11};
//+
Line(11) = {11, 3};
//+
Line(12) = {1, 9};
//+
Line(13) = {9, 10};
//+
Line(14) = {10, 2};
//+
Line(15) = {11, 10};
//+
Line(16) = {12, 9};
//+
Line(17) = {12, 8};
//+
Line(18) = {11, 7};
//+
Line(19) = {6, 10};
//+
Line(20) = {5, 9};
//+
Curve Loop(1) = {9, 16, -12, -1};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {10, 11, 4, 9};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {2, 3, 4, 1};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {2, -14, -13, -12};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {19, -13, -20, 6};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {7, 8, 5, 6};
//+
Plane Surface(6) = {6};
//+
Curve Loop(7) = {17, 5, 20, -16};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {18, -7, 19, -15};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {11, -3, -14, -15};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {10, 18, 8, -17};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {10, 15, -13, -16};
//+
Surface(11) = {11};
//+
Surface Loop(1) = {3, 4, 9, 2, 1, 11};
//+
Volume(1) = {1};
//+
Surface Loop(2) = {10, 8, 6, 7, 5, 11};
//+
Volume(2) = {2};
//+
Physical Volume("M:S; 580;1000;4187;", 21) = {1};
//+
Physical Volume("M:L; 580;1000;4187;", 22) = {2};
//+
Physical Surface("D:0.0", 23) = {3};
//+
Physical Surface("D:2.0", 24) = {6};
//+
//Physical Surface("D", 26) = {3};
//+
Physical Surface("S:1.0", 25) = {11};
