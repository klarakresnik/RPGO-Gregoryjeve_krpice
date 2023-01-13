% a0 = [sqrt(2)/ 2, 0, sqrt(2)/ 2];
% a3 = [sqrt(2)/ 2, 0, sqrt(2)/ 2];
% b0 = [sqrt(2)/ 2, 0, sqrt(2)/ 2];
% b3 = [sqrt(2)/ 2, 0, sqrt(2)/ 2];
% c0 = [0 1 0];
% c1 = [0 1 0];
% c2 = [0 1 0];
% 
% [a1, a2] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
% ----------------------------------------------------------

Bx = meshgrid([0 1 2 3], [0 0 0 0]);
By = meshgrid([0 1 2 3], [0 0 0 0])';
Bz = zeros(4, 4);

k = 1;
alpha = pi / 4;
d = 20;

a0 = k * [0, cos(alpha), sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [1 0 0];
c1 = [1 0 0];
c2 = [1 0 0];
[b11v0, b12v0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b11v0 = b11v0 + [1 0 0];
b12v0 = b12v0 + [2 0 0];

a0 = k * [-cos(alpha), 0, sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 1 0];
c1 = [0 1 0];
c2 = [0 1 0];
[b12u1, b22u1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b12u1 = b12u1 + [3 1 0];
b22u1 = b22u1 + [3 2 0];

a0 = k * [0, -cos(alpha), sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [-1 0 0];
c1 = [-1 0 0];
c2 = [-1 0 0];
[b21v1, b22v1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b21v1 = b21v1 + [2 3 0];
b22v1 = b22v1 + [1 3 0];

a0 = k * [cos(alpha), 0, sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 -1 0];
c1 = [0 -1 0];
c2 = [0 -1 0];
[b21u0, b11u0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b21u0 = b21u0 + [0 2 0];
b11u0 = b11u0 + [0 1 0];

b11 = @(u, v) (v .* b11u0 + u .* b11v0) ./ (u + v);
b21 = @(u, v) ((1 - v) .* b21u0 + u .* b21v1) ./ (1 - v + u);
b12 = @(u, v) (v .* b12u1 + (1 - u) .* b12v0) ./(v + 1 - u);
b22 = @(u, v) ((1 - v) .* b22u1 + (1 - u) .* b22v1) ./ (1 - u + 1 - v);

[u,v] = deal(linspace(0,1,d));
u = u';
v = v';

bx = nan(length(v), length(u));
by = nan(length(v), length(u));
bz = nan(length(v), length(u));

for i = (1 : length(u))
    for j = (1 : length(v))
        b11_tocka = b11(u(i), v(j));
        b21_tocka = b21(u(i), v(j));
        b12_tocka = b12(u(i), v(j));
        b22_tocka = b22(u(i), v(j));

        Bx(2, 2) = b11_tocka(1);
        By(2, 2) = b11_tocka(2);
        Bz(2, 2) = b11_tocka(3);

        Bx(3, 2) = b21_tocka(1);
        By(3, 2) = b21_tocka(2);
        Bz(3, 2) = b21_tocka(3);

        Bx(2, 3) = b12_tocka(1);
        By(2, 3) = b12_tocka(2);
        Bz(2, 3) = b12_tocka(3);

        Bx(3, 3) = b22_tocka(1);
        By(3, 3) = b22_tocka(2);
        Bz(3, 3) = b22_tocka(3);

        [x, y, z] = bezier2(Bx, By, Bz, u(i), v(j));

        bx(j, i) = x;
        by(j, i) = y;
        bz(j, i) = z;
    end
end

hold on
axis equal
surf(bx,by,bz)

bz1= -bz - 3;
hold on
surf(bx,by,bz1)



Bx = meshgrid([0 1 2 3], [0 0 0 0]);
By = zeros(4, 4);
Bz = meshgrid([0 -1 -2 -3], [0 0 0 0])';

beta = pi / 2 - alpha;

a0 = k * [0, -sin(beta), -cos(beta)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [1 0 0];
c1 = [1 0 0];
c2 = [1 0 0];
[b11v0, b12v0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);

b11v0 = b11v0 + [1 0 0];
b12v0 = b12v0 + [2 0 0];

a0 = k * [-cos(beta), -sin(beta), 0];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 0 -1];
c1 = [0 0 -1];
c2 = [0 0 -1];
[b12u1, b22u1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b12u1 = b12u1 + [3 0 -1];
b22u1 = b22u1 + [3 0 -2];

a0 = k * [0, -sin(beta), cos(beta)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [-1 0 0];
c1 = [-1 0 0];
c2 = [-1 0 0];
[b21v1, b22v1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b22v1 = b22v1 + [2 0 -3];
b21v1 = b21v1 + [1 0 -3];

a0 = k * [cos(beta), -sin(beta), 0];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 0 1];
c1 = [0 0 1];
c2 = [0 0 1];
[b21u0, b11u0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b21u0 = b21u0 + [0 0 -2];
b11u0 = b11u0 + [0 0 -1];

b11 = @(u, v) (v .* b11u0 + u .* b11v0) ./ (u + v);
b21 = @(u, v) ((1 - v) .* b21u0 + u .* b21v1) ./ (1 - v + u);
b12 = @(u, v) (v .* b12u1 + (1 - u) .* b12v0) ./(v + 1 - u);
b22 = @(u, v) ((1 - v) .* b22u1 + (1 - u) .* b22v1) ./ (1 - u + 1 - v);

[u,v] = deal(linspace(0,1,d));
u = u';
v = v';

bx = nan(length(v), length(u));
by = nan(length(v), length(u));
bz = nan(length(v), length(u));

for i = (1 : length(u))
    for j = (1 : length(v))
        b11_tocka = b11(u(i), v(j));
        b21_tocka = b21(u(i), v(j));
        b12_tocka = b12(u(i), v(j));
        b22_tocka = b22(u(i), v(j));

        Bx(2, 2) = b11_tocka(1);
        By(2, 2) = b11_tocka(2);
        Bz(2, 2) = b11_tocka(3);

        Bx(3, 2) = b21_tocka(1);
        By(3, 2) = b21_tocka(2);
        Bz(3, 2) = b21_tocka(3);

        Bx(2, 3) = b12_tocka(1);
        By(2, 3) = b12_tocka(2);
        Bz(2, 3) = b12_tocka(3);

        Bx(3, 3) = b22_tocka(1);
        By(3, 3) = b22_tocka(2);
        Bz(3, 3) = b22_tocka(3);

        [x, y, z] = bezier2(Bx, By, Bz, u(i), v(j));
        hold on
        plot3(Bx,By,Bz, 'o')

        bx(j, i) = x;
        by(j, i) = y;
        bz(j, i) = z;
    end
end

hold on
surf(bx,by,bz)

by1 = -by + 3;

surf(bx, by1, bz)










Bx = zeros(4, 4);
By = meshgrid([0 1 2 3], [0 0 0 0]);
Bz = meshgrid([0 -1 -2 -3], [0 0 0 0])';

beta = pi / 2 - alpha;

a0 = k * [-cos(beta), 0, -sin(beta)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 1 0];
c1 = [0 1 0];
c2 = [0 1 0];
[b11v0, b12v0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b11v0 = b11v0 + [0 1 0];
b12v0 = b12v0 + [0 2 0];


a0 = k * [-cos(alpha), -sin(alpha), 0];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 0 -1];
c1 = [0 0 -1];
c2 = [0 0 -1];
[b12u1, b22u1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b12u1 = b12u1 + [0 3 -1];
b22u1 = b22u1 + [0 3 -2];

a0 = k * [-cos(beta), 0, sin(beta)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 1 0];
c1 = [0 1 0];
c2 = [0 1 0];
[b21v1, b22v1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b21v1 = b21v1 + [0 2 -3];
b22v1 = b22v1 + [0 1 -3];

a0 = k * [-cos(alpha), sin(alpha), 0];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 0 1];
c1 = [0 0 1];
c2 = [0 0 1];
[b21u0, b11u0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b21u0 = b21u0 + [0 0 -2];
b11u0 = b11u0 + [0 0 -1];

b11 = @(u, v) (v .* b11u0 + u .* b11v0) ./ (u + v);
b21 = @(u, v) ((1 - v) .* b21u0 + u .* b21v1) ./ (1 - v + u);
b12 = @(u, v) (v .* b12u1 + (1 - u) .* b12v0) ./(v + 1 - u);
b22 = @(u, v) ((1 - v) .* b22u1 + (1 - u) .* b22v1) ./ (1 - u + 1 - v);

[u,v] = deal(linspace(0,1,d));
u = u';
v = v';

bx = nan(length(v), length(u));
by = nan(length(v), length(u));
bz = nan(length(v), length(u));

for i = (1 : length(u))
    for j = (1 : length(v))
        b11_tocka = b11(u(i), v(j));
        b21_tocka = b21(u(i), v(j));
        b12_tocka = b12(u(i), v(j));
        b22_tocka = b22(u(i), v(j));

        Bx(2, 2) = b11_tocka(1);
        By(2, 2) = b11_tocka(2);
        Bz(2, 2) = b11_tocka(3);

        Bx(3, 2) = b21_tocka(1);
        By(3, 2) = b21_tocka(2);
        Bz(3, 2) = b21_tocka(3);

        Bx(2, 3) = b12_tocka(1);
        By(2, 3) = b12_tocka(2);
        Bz(2, 3) = b12_tocka(3);

        Bx(3, 3) = b22_tocka(1);
        By(3, 3) = b22_tocka(2);
        Bz(3, 3) = b22_tocka(3);

        [x, y, z] = bezier2(Bx, By, Bz, u(i), v(j));

        bx(j, i) = x;
        by(j, i) = y;
        bz(j, i) = z;
    end
end

hold on
surf(bx,by,bz)

bx1 = -bx + 3;
surf(bx1,by,bz)

