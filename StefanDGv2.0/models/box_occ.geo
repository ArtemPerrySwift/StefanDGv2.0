// Gmsh project created on Fri Jul 04 19:53:35 2025
SetFactory("OpenCASCADE");
//+
Box(1) = {0, 0, 0, 1, 1, 1};
//+
MeshSize{ PointsOf{ Volume{1}; } } = 1;
//+
Physical Volume("M:L; 1.7;1.5;1.8;", 14) = {1};
//+
Physical Surface("D", 13) = {1, 2};