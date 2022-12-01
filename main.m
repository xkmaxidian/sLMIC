%1024
clc;
clear;
close all;
addpath('MeasureTools');
addpath('function');
load('data\jDRdata.mat');
iter = 100;
M = 1;
nmix = zeros(M,1);
accr = zeros(M,1);
Fscore = zeros(M,1);
ar = zeros(M,1);
spa = zeros(M,1);
pur = zeros(M,1);
RandIdx = zeros(M,1);
Prec = zeros(M,1);
Rec = zeros(M,1);
views = 2;
gt = true_label;
X = {X1,X2};
k = size(unique(gt),1);
clusters = size(unique(gt),1);
N = size(gt,1);
err = zeros(iter,1);
% for j=1:size(X1,1)
%     g=X1(j,:);
%     [n1,v1]=find(g~=0);
%     shu(j)=length(v1);
% end
% [v2,n2]=find(shu<=30);
% X1(n2,:)=[];
% shu=[];
% for j=1:size(X2,1)
%     g=X2(j,:);
%     [n1,v1]=find(g~=0);
%     shu(j)=length(v1);
% end
% [v2,n2]=find(shu<=30);
% X2(n2,:)=[];

% X1=Normalize_TPM(X1);
% X2=Normalize_TPM(X2);
normData = 1;
    if normData == 1
        for i = 1:2
            for  j = 1:N
                normItem = std(X{i}(:,j));
                if (0 == normItem)
                    normItem = eps;
                end;
                X{i}(:,j) = (X{i}(:,j)-mean(X{i}(:,j)))/(normItem);
            end;
        end;
    end;
    
wwhh = 1;
lambda_1 =0.1 %1;%0.1 1
lambda_c = 0.1;%0.1 1
lambda_d = 1;    %0.0001 1
lambda_2 = 1;
lambda_3 = 1;
% X{1}=X{1}';
% X{2}=X{2}'
% X{3}=X{3}';
% X{4}=X{4}';
% plot(err(1:30),'r-o','MarkerFaceColor','r')
Data = X;

for why = 1:wwhh
    M=1;
    kk = size(unique(gt),1);
    nmix = zeros(M,1);
    accr = zeros(M,1);
    Fscore = zeros(M,1);
    ar = zeros(M,1);

    for m = 1:M
        tol=10^-7;
        clusters = size(unique(gt),1);
        N = size(gt,1);
        mu = 10^-6;
        rho=1.1;
        maxmu=10^6;
        I_1 = ones(N,1);
        S = rands(N,N);
        Z = zeros(N,N);
        J = zeros(N,N);
        Y2 = zeros(N,N);
        w = ones(views,1);
        E_21 = zeros(views,1);

    %% initialize S0: Constructing the SIG matrices
        pn = 5;
        S0 = cell(1,views);
        WW = cell(1,views);
        HH = cell(1,views);
        options = [];
        option.Metric = 'Cosine';
        options.NeighborMode = 'KNN';%KNN
        options.k =5;%5 nearest neighbors
        options.WeightMode = 'Cosine';%Weights are 0 or 1, it can eplace with 'HeatKernel', 'Euclidean' 
        for i = 1:views
            %% NMF
            %[WW{i}, HH{i}] = NMF(Data{i},300,100);
            %[S0{i}, ~] = InitializeSIGs(HH{i}, pn, 0);
            %[S0{i}, ~] = InitializeSIGs(Data{i}, pn, 0);
            %S0{i} = constructW(HH{i}',options)
            S0{i} = constructW(Data{i}',options);
        end;
        S0_initial = S0;
        %% constructing PMI matrices
        P0 = cell(1,views);
        kp = 1;
        for i = 1:views
            P0{i} = PMI(S0{i},kp);
        end;    
        P0_initial = P0;
            for v = 1:views
                X{v} = Data{v}./repmat(sqrt(sum(Data{v}.^2,1)),size(Data{v},1),1);
                %X{v} = S0_initial{v};
                %X{v}(all(X{v}==0,2),:) = [];
                %C{v} = S0_initial{v};
                C{v} = rands(N,N);
                %C{v} = 0.00001.*S0_initial{v};
                E{v} = zeros(size(X{v}));
                Y1{v} = zeros(size(X{v}));
            end
            x = zeros(views, 1);
            for i = 1:views
                x(i) = norm(Data{i},'fro')^2;
            end
            xx=sum(x);
            idxx = cell(1,m);
            ed = cell(1,m);
            for v = 1:m
                ed{v} = L2_distance_1(X{v}, X{v});
                [~, idxx{v}] = sort(ed{v}, 2); % sort each row
            end;
            for s = 1:iter
                tic;
                %% Update C(v)
                for i =1:views
                    sum_C_v = zeros(N,N);
                    for j=1:views
                        if j~=i
                            sum_C_v = sum_C_v + w(j).*C{j};
                        end    
                    end    
                    Ac = zeros(N,N);
                    Cc = zeros(N,N);
                    Ac = mu.*X{i}'*X{i} +2*lambda_3*w(i)*w(i).*eye(N,N);
                    Cc = mu.*X{i}'*(X{i}-X{i}*S-E{i})+X{i}'*Y1{i}-lambda_2*w(i).*sum_C_v;
                    C{i} = max(Ac\Cc,0);  
                    %C{i}(C{i}<0) = 0;
                    %% ---E(v)  二一范数优化
                    tempP = X{i}-X{i}*S-X{i}*C{i} + Y1{i}/mu;
                    E{i} = solve_l1l2(tempP,1/mu);
                    %% ---update w
                    for nn = 1:views
                        w(v) = 0.5/norm(S - C{v},'fro');
                    end    
                end    
                w=w./sum(w);
                w(w<0)=0.00001;         
                %% Update S
                A = zeros(N,N);
                B = zeros(N,N);
                for i = 1:views
                    A = A + X{i}'*X{i};
                    B = B + X{i}'*(X{i}+Y1{i}/mu-X{i}*C{i}-E{i});
                end    
                S = max((A + eye(N,N))\(B + Z -Y2/mu + lambda_d.*J),0);

                %% Update Z
                [U,SS,V] = svd(S + Y2./mu);
                a = diag(SS)-lambda_1/mu;
                a(a<0)=0;
                T = diag(a);
                Z = U*T*V';

                %% ---W(ibdr)
                L_C = (abs(J) + abs(J')).*0.5; 
                [W_c,~] = eig(L_C);
                
                len = min(length(W_c(1,:)),kk);
                W_c = W_c(:,len)*W_c(:,len)';
                %% ---J
                tmp_W = diag(W_c)*I_1'-W_c;
                tmp = abs(S)-0.5*lambda_c*lambda_d.*(tmp_W+tmp_W');
                J = max(0,tmp);
                J = sign(S).*J;    
                %% converge
                for i = 1:views
                    p(i) = norm(X{i}-X{i}*(S+C{i}),'fro')^2;
                    eq1{i} = X{i}-X{i}*(S+C{i})-E{i};
                end
                P = sum(p);
                eq2=S-Z;  
                err(s)=abs(xx-P);
                if s==iter || err(s)<tol
            %             break;
                else
                    xx = P;
                    for i=1:views
                        Y1{i}=Y1{i}+mu*eq1{i};
                    end
                    Y2 = Y2 + mu*eq2;
                    mu = min(maxmu,mu*rho);
                end
                s = s+1;
            end        
             %% clustering
            Zt=zeros(N,N);
            for i = 1:views
                Zt=Zt+w(i)*(abs(C{i})+abs(C{i})');
            end
            Z = (abs(S)+abs(S'))/2+Zt/(2*2);
            Clus = SpectralClustering(Z,clusters);
            time = toc;
            ll=gt;%%%  the label originally identified by the authors
            l=Clus;%%% Labels obtained by DRjCC
            %[ACC MIhat Purity F P R RI];
%             result = ClusteringMeasure_new(ll, l);  
%             nmix(m) = result(2);
%             accr(m)  = result(1);
%             Fscore(m)  = result(4);
%             ar(m)  = result(7);
            [nmi,ACC,f,RI]=clustering(Z, clusters, gt);
            nmix(m) = nmi;
            accr(m)  = ACC;
            Fscore(m)  = f;
            ar(m)  = RI;
            resultt = [mean(accr),mean(nmix),mean(ar),mean(Fscore)]
            resultts = [std(accr),std(nmix),std(ar),std(Fscore)]
    end 
end    
        
% Z = (abs(S)+abs(S'))/2;
Z = cnormalize(Z, Inf);
% A = abs(Z_n) + abs(Z_n)';
% imshow(A)
% A(A<0.15) = 0;
% spy(A)
% Z = (abs(S)+abs(S'))/2;
% Z_n = cnormalize(S, Inf);
% A = abs(Z_n) + abs(Z_n)';
% imshow(A)
% Z(Z<0.03) = 0;
% spy(Z)
% spy(A)
% %plot(err(1:30),'r-o','MarkerFaceColor','r')
% %save('ORL_err.mat','err')
% S0_initial{1}(logical(eye(size(S0_initial{1}))))=0
% ss = (S0_initial{1}+S0_initial{1}')/2;
% graph(ss(180:210,180:210));
% plot(ans);