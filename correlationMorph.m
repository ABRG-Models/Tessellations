%To show histograms of adjacent and randomly matched edges from a single
%director
load ('./logsMorph/correlate1.data');
load ('./logsMorph/random_correlate1.data');
figure(1)
h1 = histogram(random_correlate1,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0.0 Both morph 0')
xlabel('correlation of edges')
ylabel('frequency')
hold
h2 = histogram(correlate1,80)
h2.Normalization = 'probability'
h2.FaceColor = 'r'
h2.FaceAlpha = 0.4
print ( './logsMorph/corrMorph0Both.png','-dpng')

figure(2)
h1 = histogram(random_correlate1,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Rand morph 0')
xlabel('correlation of edges')
ylabel('frequency')
print ( './logsMorph/corrMorph0Rand.png','-dpng')

figure(3)
h1 = histogram(correlate1,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Adj morph 0')
xlabel('correlation of edges')
ylabel('frequency')
h1.FaceColor = 'r'
print ( './logsMorph/corrMorph0Adj.png','-dpng')


load ('./logsMorph/correlate2.data');
load ('./logsMorph/random_correlate2.data');
figure(4)
h1 = histogram(random_correlate2,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Both morph 1')
xlabel('correlation of edges')
ylabel('frequency')
hold
h2 = histogram(correlate2,80)
h2.Normalization = 'probability'
h2.FaceColor = 'r'
h2.FaceAlpha = 0.4
print ( './logsMorph/corrMorph1Both.png','-dpng')

figure(5)
h1 = histogram(random_correlate2,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Rand morph 1')
xlabel('correlation of edges')
ylabel('frequency')
print ( './logsMorph/corrMorph1Rand.png','-dpng')

figure(6)
h1 = histogram(correlate2,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Adj morph 1')
xlabel('correlation of edges')
ylabel('frequency')
h1.FaceColor = 'r'

print ( './logsMorph/corrMorph1Adj.png','-dpng')



load ('./logsMorph/correlate3.data');
load ('./logsMorph/random_correlate3.data');
figure(7)
h1 = histogram(random_correlate3,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Both morph 2')
xlabel('correlation of edges')
ylabel('frequency')
hold
h2 = histogram(correlate3,80)
h2.Normalization = 'probability'
h2.FaceColor = 'r'
h2.FaceAlpha = 0.4
print ( './logsMorph/corrMorph2Both.png','-dpng')

figure(8)
h1 = histogram(random_correlate3,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Rand morph 2')
xlabel('correlation of edges')
ylabel('frequency')
print ( './logsMorph/corrMorph2Rand.png','-dpng')

figure(9)
h1 = histogram(correlate3,80)
h1.Normalization = 'probability'
title('Correlations for Dn = 36.0 Chi = 0 Adj morph 2')
xlabel('correlation of edges')
ylabel('frequency')
h1.FaceColor = 'r'
print ( './logsMorph/corrMorph2Adj.png','-dpng')






