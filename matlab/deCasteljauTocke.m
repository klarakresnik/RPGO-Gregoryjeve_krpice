function tocke = deCasteljauTocke(b,t)
% DE_CASTELJAU izracuna tocko na Bezierovi krivulji pri parametru t
% s pomocjo de Casteljaujevega algoritma.
% Stolpci matrike b so kontrolne tocke Bezierove krivulje.

b = b';
tocke = zeros(length(t), length(b(:, 1)));
n = size(b, 2) - 1;
B = nan(1, 2);
for j = 1 : length(t)
    pomozen_vektor = b;
    for i = (1 : n)
        B = (1 - t(j)) * pomozen_vektor(:, 1 : end - 1) + t(j) * pomozen_vektor(:, 2 : end);
        pomozen_vektor = B;
    end
    tocke(j, :) = B;
end
end