
function NX=Normalize_TPM(X)
%%%=========Here, for scRNA-seq data, and normalize such that the total sum of counts, X,
%%%====is 10^6 in each cell j, which is essentially the Transcript per Million (TPM) normalization.
%         X=csvread('C:\\Users\\Administrator\\R\\data pre-processing\\Splate2_10.csv');
%         X=X.data
        %X=X';
        for i=1:size(X,1)
            x=X(i,:);
            n=find(x~=0);
            l=length(n);
            x=x./l;
           xX(i,:)=x;
        end
       ll=sum(xX);
       NX = xX./ll;
%        csvwrite('Splat2_tpm.csv',NX)
       NX=NX*1e6;
       size(NX);
%        csvwrite('1.csv',NX)
end
% function NX=Normalize_TPM(X)
% %%%=========Here, for scRNA-seq data, and normalize such that the total sum of counts, X,
% %%%====is 10^6 in each cell j, which is essentially the Transcript per Million (TPM) normalization.
%         X=importdata('C:\\Users\\Administrator\\Desktop\\gene_HUMAN3\\GSM2230759_human3_umifm_counts.csv');
%         X=X.data
%         for i=1:size(X,1)
%             x=X(i,:);
%             n=find(x~=0);
%             l=length(n);
%             x=x./l;
%            xX(i,:)=x;
%         end
%        ll=sum(xX);
%        NX = xX./ll;
%        csvwrite('human3_tpm.csv',NX)
%        NX=NX*1e6;
%        size(NX);
%        csvwrite('1.csv',NX)
% end












