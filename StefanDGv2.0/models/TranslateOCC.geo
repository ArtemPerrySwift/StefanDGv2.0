// Gmsh project created on Tue Oct 20 11:06:43 2026

// Создаем точки плоскости
Point(1) = {0, 0, 0};
Point(2) = {1, 0, 0};
Point(3) = {1, 1, 0};
Point(4) = {0, 1, 0};

// Строим контур и саму плоскость
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};//+
Translate {1, 0, 0} {
  Point{3}; 
}
