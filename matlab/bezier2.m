function [bx,by,bz] = bezier2(Bx,By,Bz,u,v)
% Opis:
% bezier2 vrne točke na Bezierjevi ploskvi iz tenzorskega produkta
%
% Definicija:
% [bx,by,bz] = bezier2(Bx,By,Bz,u,v)
%
% Vhodni podatki:
% Bx, By, Bz matrike velikosti n+1 x m+1, ki predstavljajo
% koordinate kontrolnih točk,
% u, v vrstici dolžine M in N, ki predstavljata parametre
% v smereh u in v
%
% Izhodni podatki:
% bx, by, bz matrike velikosti N x M, ki predstavljajo točke na
% Bezierjevi ploskvi: [bx(J,I) by(J,I) bz(J,I)] je
% točka pri parametrih u(I) in v(J)

[n, m] = size(Bx);
M = length(u);
N = length(v);
bx = nan(N, M);
by = nan(N, M);
bz = nan(N, M);

vmesne_kontrolne_tocke_x = nan(N, m);
vmesne_kontrolne_tocke_y = nan(N, m);
vmesne_kontrolne_tocke_z = nan(N, m);

for i = (1 : m)
    b = [Bx(:, i), By(:, i), Bz(:, i)];
    tocke = deCasteljauTocke(b, v);
    vmesne_kontrolne_tocke_x(:, i) = tocke(:, 1); 
    vmesne_kontrolne_tocke_y(:, i) = tocke(:, 2); 
    vmesne_kontrolne_tocke_z(:, i) = tocke(:, 3);
end

for j = (1 : N)
    b = [vmesne_kontrolne_tocke_x(j, :)', vmesne_kontrolne_tocke_y(j, :)', vmesne_kontrolne_tocke_z(j, :)'];
    tocke = deCasteljauTocke(b, u);
    bx(j, :) = tocke(:, 1)';
    by(j, :) = tocke(:, 2)';
    bz(j, :) = tocke(:, 3)';
end



end