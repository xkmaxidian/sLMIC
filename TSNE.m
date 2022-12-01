%s=s'
%Z=xlsread('E:\\code data\\clusterGAN_OK\\biase\\ClusterGAN-maste0303\\ClusterGANmaster\\test_z.csv')

%m=xlsread('D:\\E≈Ã\\Compare methods\\dca-master\\dca\\dca_data.csv');
%m(any(isnan(m)'),:) = [];
%s = xlsread('D:\\scLDS2\\Code\\Human3 - ∏±±æ\\ClusterGAN-maste0222\\ClusterGANmaster\\test_z.csv')
% Z=xlsread('E:\\Compare methods\\compare methodss\\Dimension_reduction\\mESC_feature_NMF_RD.csv')
%E:\\code data\\clusterGAN_OK\\human2\\ClusterGAN-maste1225\\ClusterGANmaster\\test_z.csv
%E:\\Compare methods\\dca-master\\dca\\dca_data.csv
%Z=Z'
%V=V'
% s=s'
% Z(Z<0.01) = 0;
% Z = ( abs(A) + abs(A') ) / 2 
%Z=[C1 C2];


%Z = xlsread('D:\\E≈Ã\\code data\\ClusterGAN_Sc_data\\mouse2_shaixuan10\\GSM2230762_mouse2_10_tpm.csv')
% Z = rmmissing(Z) 
% M=xlsread('D:\\E≈Ã\\code data\\TSNE_DATA\\all\\human3_all\\tSNEM_human3_z.csv');
%s(any(isnan(s)'),:) = []
%  Z = rmmissing(Z); 
%  mappedX = tsne(Z);
% %color = lines(7);
%,color(1:14,:)
%Z = xlsread('D:\\E≈Ã\\code data\\ClusterGAN_Sc_data\\human2_shaixuan10\\GSM2230758_human2_10_tpm.csv');
mappedX = tsne(Z);
gscatter(mappedX(:,1), mappedX(:,2),gt)