// Planar biased-probe plasma chamber.
// Regenerate with:
//   gmsh -2 -format msh2 -o biased_probe_2d.msh biased_probe_2d.geo

Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;
Mesh.Algorithm = 6;

chamber_length = 0.12;
chamber_half_height = 0.04;
probe_x = 0.075;
probe_radius = 0.012;
bulk_size = 0.006;
probe_size = 0.002;

Point(1) = {0, -chamber_half_height, 0, bulk_size};
Point(2) = {chamber_length, -chamber_half_height, 0, bulk_size};
Point(3) = {chamber_length, chamber_half_height, 0, bulk_size};
Point(4) = {0, chamber_half_height, 0, bulk_size};

Point(5) = {probe_x + probe_radius, 0, 0, probe_size};
Point(6) = {probe_x, probe_radius, 0, probe_size};
Point(7) = {probe_x - probe_radius, 0, 0, probe_size};
Point(8) = {probe_x, -probe_radius, 0, probe_size};
Point(9) = {probe_x, 0, 0, probe_size};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Circle(5) = {5, 9, 6};
Circle(6) = {6, 9, 7};
Circle(7) = {7, 9, 8};
Circle(8) = {8, 9, 5};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(1) = {1, 2};

Physical Curve("wall", 1) = {1, 3};
Physical Curve("outlet", 2) = {2};
Physical Curve("inlet", 3) = {4};
Physical Curve("probe", 4) = {5, 6, 7, 8};
Physical Surface("plasma", 10) = {1};

Field[1] = Distance;
Field[1].CurvesList = {5, 6, 7, 8};
Field[1].Sampling = 100;
Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = probe_size;
Field[2].SizeMax = bulk_size;
Field[2].DistMin = 0.006;
Field[2].DistMax = 0.025;
Background Field = 2;

Mesh.MeshSizeExtendFromBoundary = 0;
Mesh.MeshSizeFromPoints = 1;
Mesh.MeshSizeFromCurvature = 12;
