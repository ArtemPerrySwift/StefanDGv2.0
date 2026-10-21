// 3. –азбиваем сферу плоскостью
// ѕервый аргумент в { } Ч что режем (объем 1), второй Ч чем режем (поверхность 2)
BooleanFragments{ Volume{1}; Delete; }{ Surface{2}; Delete; }
//+
Recursive Delete {
  Surface{4}; 
}
