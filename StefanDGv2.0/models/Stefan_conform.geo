// Gmsh project created on Thu Oct 15 12:42:27 2026
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
Point(9) = {2, 0, 0, 1.0};
//+
Point(10) = {2, 1, 0, 1.0};
//+
Point(11) = {2, 1, 1, 1.0};
//+
Point(12) = {2, 0, 1, 1.0};
//+
Line(13) = {5, 12};
//+
Line(14) = {12, 9};
//+
Line(15) = {9, 6};
//+
Line(16) = {12, 11};
//+
Line(17) = {11, 7};
//+
Line(18) = {11, 10};
//+
Line(19) = {9, 10};
//+
Line(20) = {10, 8};
//+
Curve Loop(7) = {13, 16, 17, -6};
//+
Plane Surface(7) = {7};
//+
Curve Loop(8) = {13, 14, 15, 5};
//+
Plane Surface(8) = {8};
//+
Curve Loop(9) = {16, 18, -19, -14};
//+
Plane Surface(9) = {9};
//+
Curve Loop(10) = {15, 8, -20, -19};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {18, 20, 7, -17};
//+
Plane Surface(11) = {11};
//+
Surface Loop(2) = {7, 8, 9, 11, 10, 2};
//+
Volume(2) = {2};
//+
div = 0.5;
//+
MeshSize{ PointsOf{ Volume{1}; } } = div;
MeshSize{ PointsOf{ Volume{2}; } } = div;
//+
Physical Surface("D:0.0", 21) = {1};
//+
Physical Surface("D:2.0", 22) = {9};
//+
Physical Surface("S:1.0", 23) = {2};
//+
//Physical Surface("D", 23) = {1, 2, 9};
//+
Physical Volume("M:S; 580;1000;4187;", 24) = {1};
//+
Physical Volume("M:L; 580;1000;4187;", 25) = {2};
