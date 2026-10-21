// Gmsh project created on Tue Oct 20 12:20:19 2026
SetFactory("OpenCASCADE");
//+
Point(1) = {0, 0, 0, 1.0};
//+
Point(2) = {0, 1, 0, 1.0};
//+
Point(3) = {0, 1, 1, 1.0};
//+
Point(4) = {0, 0, 1, 1.0};
//+
Point(5) = {1, 0, 0, 1.0};
//+
Point(6) = {1, 1, 0, 1.0};
//+
Point(7) = {1, 1, 1, 1.0};
//+
Point(8) = {1, 0, 1, 1.0};
//+
Line(1) = {3, 7};
//+
Line(2) = {7, 6};
//+
Line(3) = {6, 2};
//+
Line(4) = {2, 3};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Translate {0.1, 0, 0} {
  Point{7}; 
}
