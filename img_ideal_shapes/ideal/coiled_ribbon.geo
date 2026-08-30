SetFactory("OpenCASCADE");

//Read STEP file and store volume(s)
v[] = ShapeFromFile("coiled_ribbon.step");

//Mesh only the surface
Physical Surface("surf") = CombinedBoundary{ Volume{v[]}; };
Delete { Volume{v[]};}
