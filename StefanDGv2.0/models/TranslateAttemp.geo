// Gmsh project created on Tue Oct 20 00:35:06 2026
SetFactory("Built-in");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {1, 0, 0, 1.0};
//+
Point(3) = {1, 1, 0, 1.0};
//+
Point(4) = {0, 1, 0, 1.0};
//+
Point(5) = {0, 1, 1, 1.0};
//+
Point(6) = {1, 1, 1, 1.0};
//+
Point(7) = {1, 0, 1, 1.0};
//+
Point(8) = {0, 0, 1, 1.0};
//+
Line(1) = {5, 4};
//+
Line(2) = {4, 3};
//+
Line(3) = {3, 6};
//+
Line(4) = {6, 5};
//+
Curve Loop(1) = {4, 1, 2, 3};
//+
Plane Surface(1) = {1};
//+
Line(5) = {5, 8};
//+
Line(6) = {8, 7};
//+
Line(7) = {7, 6};
//+
Curve Loop(2) = {6, 7, 4, 5};
//+
Plane Surface(2) = {2};
//+
Line(8) = {8, 1};
//+
Line(9) = {1, 4};
//+
Curve Loop(3) = {8, 9, -1, 5};
//+
Plane Surface(3) = {3};
//+
Line(10) = {1, 2};
//+
Line(11) = {2, 3};
//+
Curve Loop(4) = {10, 11, -2, -9};
//+
Plane Surface(4) = {4};
//+
Line(12) = {7, 2};
//+
Curve Loop(5) = {12, 11, 3, -7};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {6, 12, -10, -8};
//+
Plane Surface(6) = {6};
//+
Surface Loop(1) = {2, 6, 5, 4, 1, 3};
//+
Volume(1) = {1};

//+
Point(9) = {2, 0, 0, 1.0};
//+
Point(10) = {2, 1, 0, 1.0};
//+
Point(11) = {2, 1, 1, 1.0};
//+
Point(12) = {2, 0, 1, 1.0};
//+
Line(13) = {6, 11};
//+
Line(14) = {11, 10};
//+
Line(15) = {10, 3};
//+
Curve Loop(7) = {13, 14, 15, 3};
//+
Plane Surface(7) = {7};
//+
Line(16) = {10, 9};
//+
Line(17) = {9, 2};
//+
Curve Loop(8) = {15, -11, -17, -16};
//+
Plane Surface(8) = {8};
//+
Line(18) = {7, 12};
//+
Line(19) = {12, 9};
//+
Curve Loop(9) = {18, 19, 17, -12};
//+
Plane Surface(9) = {9};
//+
Line(20) = {12, 11};
//+
Curve Loop(10) = {19, -16, -14, -20};
//+
Plane Surface(10) = {10};
//+
Curve Loop(11) = {20, -13, -7, 18};
//+
Plane Surface(11) = {11};
//+
Surface Loop(2) = {11, 10, 9, 8, 7, 5};
//+
Volume(2) = {2};
//+
Physical Surface("D:0.0", 21) = {3};
//+
Physical Surface("D:2.0", 22) = {10};
//+
Physical Surface("S:1.0", 23) = {5};
//+
Physical Volume("M:S; 580;1000;4187;", 24) = {1};
//+
Physical Volume("M:L; 580;1000;4187;", 25) = {2};
