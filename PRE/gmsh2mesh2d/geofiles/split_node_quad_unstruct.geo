// Gmsh project created on Fri Oct 22 09:32:47 2021
SetFactory("OpenCASCADE");
//+
// define some variables
rf=0.2;// mesh size on the fault
rb=3; // mesh size on boundaries
H=80;
theta=60;
L1=60;
PI=Acos(-1);
L2=L1 + H/Tan(theta*PI/180.0);

// define nodal points
Point(1) = {0, 0, 0, rb};
//+
Point(2) = {0, -H, 0, rb};
//+
Point(3) = {L2, -H, 0, rf};
//+
Point(4) = {L1, 0, 0, rf};
//+
Point(5) = {L2+L1, -H, 0, rb};
//+
Point(6) = {L2+L1, 0, 0, rb};
//+
Point(7) = {L1, 0, 0, rf};
//+
Point(8) = {L2, -H, 0, rf};

//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 1};
//+
Line(5) = {8, 5};
//+
Line(6) = {5, 6};
//+
Line(7) = {6, 7};
//+
Line(8) = {7, 8};

//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {5, 6, 7, 8};
//+
Plane Surface(2) = {2};

Physical Curve("bc1") = {2,5};
Physical Curve("bc2") = {6};
Physical Curve("bc3") = {4,7};
Physical Curve("bc4") = {1};
Physical Curve("bc5") = {3};
Physical Curve("bc6") = {8};

Physical Surface('mat1') = {1};
Physical Surface('mat2') = {2};

Recombine Surface{1, 2};
Mesh.Algorithm = 8;
Mesh.RecombinationAlgorithm = 2; // or 3
// RefineMesh;

