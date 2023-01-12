SX0 = reshape(meshgrid([0 1 2 3 0 1 2 3 0 1 2 3 0 1 2 3], [0 0 0 0]), [64, 1]);
SY0 = reshape(meshgrid([0 1 2 3], [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0])', [64, 1]);
SZ0 = reshape(meshgrid([0 1 2 3], [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]), [64, 1]);

plot3(SX0, SY0, SZ0, 'o')