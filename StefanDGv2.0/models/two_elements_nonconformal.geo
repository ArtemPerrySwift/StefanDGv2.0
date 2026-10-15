// Gmsh project created on Sat Jul 05 03:01:13 2025
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 0.5};
//+
Point(2) = {-1, 0, 0, 2.0};
//+
Point(3) = {0, 1, 0, 2.0};
//+
Point(4) = {0, 0, 1, 2.0};
//+
Line(1) = {2, 1};
//+
Line(2) = {1, 3};
//+
Line(3) = {4, 1};
//+
Line(4) = {2, 3};
//+
Line(5) = {4, 2};
//+
Line(6) = {4, 3};

//+
Curve Loop(1) = {6, -4, -5};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {6, -2, -3};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {1, -3, 5};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {4, -2, -1};
//+
Plane Surface(4) = {4};
//+
Surface Loop(1) = {1, 2, 4, 3};
//+
Volume(1) = {1};
//+
Translate {1, 0, 0} {
  Duplicata { Volume{1}; }
}//+
Rotate {{0, 1, 0}, {0, 0, 0}, Pi/2} {
  Volume{2}; 
}
//+
Physical Surface("C", 13) = {2, 7};
//+
Physical Surface("D", 14) = {1, 5, 8, 3};
//Physical Surface("D", 14) = {1, 2, 5, 7, 8, 3};
//+
Translate {0, 0, 1} {
  Volume{2}; 
}
//+
Physical Volume("M:L; 1.5;1.7;1.3;", 15) = {1, 2};

