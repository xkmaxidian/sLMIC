Z = xlsread('D:\\���ڽ��й���\\work4\\kk\\feature.csv')

mappedX = tsne(Z);
gscatter(mappedX(:,1), mappedX(:,2),Clus)