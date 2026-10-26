// Gmsh project created on Tue Oct 20 21:32:33 2026
SetFactory("OpenCASCADE");
//+
GenMeshSize = 0.5;
//+
FrontMeshSize = 0.1;
//+
InfBorderRep = 10.0;
//+
Point(1) = {0, 0, 0, GenMeshSize };
//+
Point(2) = {0, 1, 0, GenMeshSize };
//+
Point(3) = {0, 1, 1, GenMeshSize };
//+
Point(4) = {0, 0, 1, GenMeshSize };
//+
Point(5) = {InfBorderRep, 0, 0, GenMeshSize };
//+
Point(6) = {InfBorderRep, 1, 0, GenMeshSize };
//+
Point(7) = {InfBorderRep, 1, 1, GenMeshSize };
//+
Point(8) = {InfBorderRep, 0, 1, GenMeshSize };
//+
Point(9) = {1, 0, 0, FrontMeshSize };
//+
Point(10) = {1, 1, 0, FrontMeshSize };
//+
Point(11) = {1, 1, 1, FrontMeshSize };
//+
Point(12) = {1, 0, 1, FrontMeshSize };
//+
Line(1) = {4, 1};
//+
Line(2) = {1, 2};
//+
Line(3) = {2, 3};
//+
Line(4) = {3, 4};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};
//+
//Point(13) = {1, 0.5, 0.5, GenMeshSize };
//+
//Point(14) = {1, 0.5, 0.7, GenMeshSize };
//+
Physical Point("S:1:0", 9) = {10, 11, 12, 9};
//+
