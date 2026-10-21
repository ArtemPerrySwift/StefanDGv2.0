// Gmsh project created on Tue Oct 20 12:37:02 2026
SetFactory("OpenCASCADE");

Point(1) = {0, 0, 0};
Point(2) = {5, 0, 0};
Line(1) = {1, 2};
//+//+
Translate {1, 0, 0} {
  Curve{1}; 
}
