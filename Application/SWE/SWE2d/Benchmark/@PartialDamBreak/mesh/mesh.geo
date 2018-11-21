cl = 5;
Point(1) = {0, 0, 0, cl};
Point(2) = {0, 200, 0, cl};
Point(3) = {200, 200, 0, cl};
Point(4) = {200, 0, 0, cl};
Point(5) = {95, 0, 0, cl};
Point(6) = {105, 0, 0, cl};
Point(7) = {95, 87.5, 0, cl};
Point(8) = {105, 87.5, 0, cl};
Point(9) = {95, 162.5, 0, cl};
Point(10) = {105, 162.5, 0, cl};
Point(11) = {95, 200, 0, cl};
Point(12) = {105, 200, 0, cl};
Line(1) = {2, 1};
Line(2) = {1, 5};
Line(3) = {5, 7};
Line(4) = {7, 8};
Line(5) = {8, 6};
Line(6) = {6, 4};
Line(7) = {4, 3};
Line(8) = {3, 12};
Line(9) = {12, 10};
Line(10) = {10, 9};
Line(11) = {9, 11};
Line(12) = {11, 2};
Line(13) = {10, 8};
Line Loop(15) = {1, 2, 3, 4, -13, 10, 11, 12};
Line Loop(17) = {13, 5, 6, 7, 8, 9};
Plane Surface(15) = {15};
Plane Surface(17) = {17};
Physical Line(2) = {1, 7, 2, 3, 4, 5, 6, 8, 9, 10, 11, 12};
Physical Surface(4) = {15};
Physical Surface(1) = {17};

// Recombine Surface{15};
// Recombine Surface{17};
Mesh.Smoothing = 10;
