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

Bx = meshgrid([0 5 10 15], [0 0 0 0]);
By = meshgrid([0 5 10 15], [0 0 0 0])';
Bz = zeros(4, 4);

% plot3(Bx, By, Bz, 'o')
k = 0.5;
alpha = pi / 4;

a0 = k * [cos(alpha), 0, sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 5 0];
c1 = [0 5 0];
c2 = [0 5 0];
[b11v0, b12v0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);

b11v0 = b11v0 + [0 5 0];
b12v0 = b12v0 + [0 10 0];

a0 = k * [0, -cos(alpha), sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [5 0 0];
c1 = [5 0 0];
c2 = [5 0 0];
[b12u1, b22u1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b12u1 = b12u1 + [5 15 0];
b22u1 = b22u1 + [10 15 0];

a0 = k * [-cos(alpha), 0, sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [0 -5 0];
c1 = [0 -5 0];
c2 = [0 -5 0];
[b21v1, b22v1] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b22v1 = b22v1 + [15 10 0];
b21v1 = b21v1 + [15 5 0];

a0 = k * [0, cos(alpha), sin(alpha)];
a3 = a0;
b0 = a0;
b3 = a0;
c0 = [-5 0 0];
c1 = [-5 0 0];
c2 = [-5 0 0];
[b21u0, b11u0] = Chiyokura_and_Kimura(a0, a3, b0, b3, c0, c1, c2);
b21u0 = b21u0 + [10 0 0];
b11u0 = b11u0 + [5 0 0];

b11 = @(u, v) (v .* b11u0 + u .* b11v0) ./ (u + v);
b21 = @(u, v) ((1 - v) .* b21u0 + u .* b21v1) ./ (1 - v + u);
b12 = @(u, v) (v .* b12u1 + (1 - u) .* b12v0) ./(v + 1 - u);
b22 = @(u, v) ((1 - v) .* b22u1 + (1 - u) .* b22v1) ./ (1 - u + 1 - v);

[u,v] = deal(linspace(0,1,10));
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
%         plot3(Bx,By,Bz, 'o')

%         plot3(b11_tocka(1), b11_tocka(2), b11_tocka(3), 'r*')
%         hold on
%         plot3(b21_tocka(1), b21_tocka(2), b21_tocka(3), 'r*')
%         hold on
%         plot3(b12_tocka(1), b12_tocka(2), b12_tocka(3), 'r*')
%         hold on
%         plot3(b22_tocka(1), b22_tocka(2), b22_tocka(3), 'r*')

        bx(j, i) = x;
        by(j, i) = y;
        bz(j, i) = z;
    end
end

% axis([0 15 0 15 0 1])
% pbaspect([1 1 1])
hold on
surf(bx,by,bz)

theta = pi / 2;
bx1 = cos(theta) * bx + sin(theta) * bz;
by1 = by;
bz1 = -sin(theta) * bx + cos(theta) * bz;

hold on 
surf(bx1,by1,bz1)
