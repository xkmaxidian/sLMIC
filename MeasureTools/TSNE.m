Z = xlsread('D:\\正在进行工作\\work4\\kk\\feature.csv')

mappedX = tsne(Z);
gscatter(mappedX(:,1), mappedX(:,2),Clus)