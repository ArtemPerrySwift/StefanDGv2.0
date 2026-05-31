//+
SetFactory("OpenCASCADE");
//+
front = 0.1;
far = 1.0;
//+
Box(1) = {-1, 0, 0, 0.75, 1, 1};
//+
MeshSize{ PointsOf{ Volume{1}; } } = far;
//+
Box(2) = {0.25, 0, 0, 0.75, 1, 1};
//+
MeshSize{ PointsOf{ Volume{2}; } } = far;
//+
Point(17) = {-0.25, 1, 2, front};
//+
Point(18) = {-0.25, 0, 2, front};
//+
Point(19) = {-0.25, 0, 3, front};
//+
Point(20) = {-0.25, 1, 3, front};
//+
Point(21) = {0.25, 1, 2, front};
//+
Point(22) = {0.25, 0, 2, front};
//+
Point(23) = {0.25, 0, 3, front};
//+
Point(24) = {0.25, 1, 3, front};
//+
Point(25) = {0, 1, 2, front};
//+
Point(26) = {0, 0, 2, front};
//+
Point(27) = {0, 0, 3, front};
//+
Point(28) = {0, 1, 3, front};
//+
Line(25) = {19, 18};
//+
Line(26) = {18, 17};
//+
Line(27) = {17, 20};
//+
Line(28) = {20, 19};
//+
Line(29) = {22, 21};
//+
Line(30) = {21, 24};
//+
Line(31) = {24, 23};
//+
Line(32) = {23, 22};
//+
Line(33) = {26, 25};
//+
Line(34) = {25, 28};
//+
Line(35) = {28, 27};
//+
Line(36) = {27, 26};
//+
Curve Loop(13) = {28, 25, 26, 27};
//+
Plane Surface(13) = {13};
//+
Curve Loop(14) = {35, 36, 33, 34};
//+
Plane Surface(14) = {14};
//+
Curve Loop(15) = {31, 32, 29, 30};
//+
Plane Surface(15) = {15};
//+
Line(37) = {17, 25};
//+
Line(38) = {28, 20};
//+
Line(39) = {19, 27};
//+
Line(40) = {18, 26};
//+
Line(41) = {26, 22};
//+
Line(42) = {27, 23};
//+
Line(43) = {28, 24};
//+
Line(44) = {25, 21};
//+
Curve Loop(16) = {27, -38, -34, -37};
//+
Plane Surface(16) = {16};
//+
Curve Loop(17) = {28, 39, -35, 38};
//+
Plane Surface(17) = {17};
//+
Curve Loop(18) = {25, 40, -36, -39};
//+
Plane Surface(18) = {18};
//+
Curve Loop(19) = {26, 37, -33, -40};
//+
Plane Surface(19) = {19};
//+
Curve Loop(20) = {44, -29, -41, 33};
//+
Plane Surface(20) = {20};
//+
Curve Loop(21) = {29, 30, 31, 32};
//+
Plane Surface(21) = {21};
//+
Curve Loop(22) = {31, -42, -35, 43};
//+
Plane Surface(22) = {22};
//+
Curve Loop(23) = {43, -30, -44, 34};
//+
Plane Surface(23) = {23};
//+
Curve Loop(24) = {42, 32, -41, -36};
//+
Plane Surface(24) = {24};
//+
Surface Loop(3) = {17, 13, 18, 19, 16, 14};
//+
Volume(3) = {3};
//+
Surface Loop(4) = {22, 24, 20, 23, 15, 14};
//+
Volume(4) = {4};
//+
Translate {0, 0, -1} {
  Volume{3}; Volume{4}; 
}
//+
Translate {0, 0, -1} {
  Volume{3}; Volume{4}; 
}
//+
Recursive Delete {
  Surface{21}; 
}
//+
Physical Volume("M:L; 580;1000;4187;", 85) = {1, 3};
//+
Physical Volume("M:S; 2230;900;2060;", 86) = {4, 2};
//+
Physical Surface("S: 273.15", 87) = {39};
//+
Physical Surface("D:293.15", 88) = {1};
//+
Physical Surface("D:271.15", 89) = {8};
//+
Physical Surface("C_0", 90) = {2, 41};
//+
Physical Surface("C_1", 91) = {7, 46};
