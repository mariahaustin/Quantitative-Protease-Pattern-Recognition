clc
clear
set(0,'DefaultFigureWindowStyle','normal')
close all


%% Collagenase Curve Fitting
%Each protease and substrate are fit individually using this same repeated
%process
CollagenaseTrain = readtable('CollagenaseTrainingData.xlsx'); 
time = table2array(CollagenaseTrain(:,1));
Collagenase_Peptide = table2array(CollagenaseTrain(:,2:25));%Separate each substrate into individual arrays
Collagenase_P3NPro = table2array(CollagenaseTrain(:,27:50));
Collagenase_P1NAsn = table2array(CollagenaseTrain(:,52:75));
Collagenase_P3pNAla = table2array(CollagenaseTrain(:,77:100));

figure('Name','Collagenase Training')
subplot(2,2,1)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Collagenase_Peptide(:,i));%Fit first replicate of highest concentration curve
A(i,1)= fitresult.A;%Record plateau value
train_k(i,1) = fitresult.k;%Record exponential constant
end
Plat1 = A(1,1);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Collagenase_Peptide(:,i),Plat1);%Fit remaining curves with plateau value from first fit
A(i,1)= fitresult.A;%Record plateau value for verification
train_k(i,1) = fitresult.k;%Record exponential constants (will become array)
end
title('Peptide')

subplot(2,2,2)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Collagenase_P3NPro(:,i));
A(i,2)= fitresult.A;
train_k(i,2) = fitresult.k;
end
Plat2 = A(1,2);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Collagenase_P3NPro(:,i),Plat2);
A(i,2)= fitresult.A;
train_k(i,2) = fitresult.k;
end
title('P3 NPro')

subplot(2,2,3)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Collagenase_P1NAsn(:,i));
A(i,3)= fitresult.A;
train_k(i,3) = fitresult.k;
end
Plat3 = A(1,3);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Collagenase_P1NAsn(:,i),Plat3);
A(i,3)= fitresult.A;
train_k(i,3) = fitresult.k;
end
title('P1 NAsn')

subplot(2,2,4)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Collagenase_P3pNAla(:,i));
A(i,4)= fitresult.A;
train_k(i,4) = fitresult.k;
end
Plat4 = A(1,4);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Collagenase_P3pNAla(:,i),Plat4);
A(i,4)= fitresult.A;
train_k(i,4) = fitresult.k;
end
title('P3p NAla')
hold off

%% Proteinase K Curve Fitting
ProtKTrain = readtable('ProteinaseKTrainingData.xlsx');
time = table2array(ProtKTrain(:,1));
ProtK_Peptide = table2array(ProtKTrain(:,2:25));
ProtK_P3NPro = table2array(ProtKTrain(:,27:50));
ProtK_P1NAsn = table2array(ProtKTrain(:,52:75));
ProtK_P3pNAla = table2array(ProtKTrain(:,77:100));

figure('Name','ProtK Training')
subplot(2,2,1)
hold on
for i=[1]
[fitresult,gof] = createFit(time,ProtK_Peptide(:,i));
A(i+24,1)= fitresult.A;
train_k(i+24,1) = fitresult.k;
end
Plat1 = A(25,1);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,ProtK_Peptide(:,i),Plat1);
A(i+24,1)= fitresult.A;
train_k(i+24,1) = fitresult.k;
end
title('Peptide')

subplot(2,2,2)
hold on
for i=[1]
[fitresult,gof] = createFit(time,ProtK_P3NPro(:,i));
A(i+24,2)= fitresult.A;
train_k(i+24,2) = fitresult.k;
end
Plat2 = A(25,2);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,ProtK_P3NPro(:,i),Plat2);
A(i+24,2)= fitresult.A;
train_k(i+24,2) = fitresult.k;
end
title('P3 NPro')

subplot(2,2,3)
hold on
for i=[1]
[fitresult,gof] = createFit(time,ProtK_P1NAsn(:,i));
A(i+24,3)= fitresult.A;
train_k(i+24,3) = fitresult.k;
end
Plat3 = A(25,3);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,ProtK_P1NAsn(:,i),Plat3);
A(i+24,3)= fitresult.A;
train_k(i+24,3) = fitresult.k;
end
title('P1 NAsn')

subplot(2,2,4)
hold on
for i=[1]
[fitresult,gof] = createFit(time,ProtK_P3pNAla(:,i));
A(i+24,4)= fitresult.A;
train_k(i+24,4) = fitresult.k;
end
Plat4 = A(25,4);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,ProtK_P3pNAla(:,i),Plat4);
A(i+24,4)= fitresult.A;
train_k(i+24,4) = fitresult.k;
end
title('P3p NAla')
hold off

%% Elastase Curve Fitting
ElastaseTrain = readtable('ElastaseTrainingData.xlsx');
time = table2array(ElastaseTrain(:,1));
Elastase_Peptide = table2array(ElastaseTrain(:,2:25));
Elastase_P3NPro = table2array(ElastaseTrain(:,27:50));
Elastase_P1NAsn = table2array(ElastaseTrain(:,52:75));
Elastase_P3pNAla = table2array(ElastaseTrain(:,77:100));

figure('Name','Elastase Training')
subplot(2,2,1)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Elastase_Peptide(:,i));
A(i+48,1)= fitresult.A;
train_k(i+48,1) = fitresult.k;
end
Plat1 = A(49,1);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Elastase_Peptide(:,i),Plat1);
A(i+48,1)= fitresult.A;
train_k(i+48,1) = fitresult.k;
end
title('Peptide')

subplot(2,2,2)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Elastase_P3NPro(:,i));
A(i+48,2)= fitresult.A;
train_k(i+48,2) = fitresult.k;
end
Plat2 = A(49,2);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Elastase_P3NPro(:,i),Plat2);
A(i+48,2)= fitresult.A;
train_k(i+48,2) = fitresult.k;
end
title('P3 NPro')

subplot(2,2,3)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Elastase_P1NAsn(:,i));
A(i+48,3)= fitresult.A;
train_k(i+48,3) = fitresult.k;
end
Plat3 = A(49,3);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Elastase_P1NAsn(:,i),Plat3);
A(i+48,3)= fitresult.A;
train_k(i+48,3) = fitresult.k;
end
title('P1 NAsn')

subplot(2,2,4)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Elastase_P3pNAla(:,i));
A(i+48,4)= fitresult.A;
train_k(i+48,4) = fitresult.k;
end
Plat4 = A(49,4);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Elastase_P3pNAla(:,i),Plat4);
A(i+48,4)= fitresult.A;
train_k(i+48,4) = fitresult.k;
end
title('P3p NAla')
hold off



%% Thermolysin Curve Fitting
ThermolysinTraining = readtable('ThermolysinTrainingData.xlsx');
time = table2array(ThermolysinTraining(:,1));
Thermolysin_Peptide = table2array(ThermolysinTraining(:,2:25));
Thermolysin_P3NPro = table2array(ThermolysinTraining(:,27:50));
Thermolysin_P1NAsn = table2array(ThermolysinTraining(:,52:75));
Thermolysin_P3pNAla = table2array(ThermolysinTraining(:,77:100));

figure('Name','Thermolysin Training')
subplot(2,2,1)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Thermolysin_Peptide(:,i));
A(i+72,1)= fitresult.A;
train_k(i+72,1) = fitresult.k;
end
Plat1 = A(73,1);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Thermolysin_Peptide(:,i),Plat1);
A(i+72,1)= fitresult.A;
train_k(i+72,1) = fitresult.k;
end
title('Peptide')

subplot(2,2,2)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Thermolysin_P3NPro(:,i));
A(i+72,2)= fitresult.A;
train_k(i+72,2) = fitresult.k;
end
Plat2 = A(73,2);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Thermolysin_P3NPro(:,i),Plat2);
A(i+72,2)= fitresult.A;
train_k(i+72,2) = fitresult.k;
end
title('P3 NPro')

subplot(2,2,3)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Thermolysin_P1NAsn(:,i));
A(i+72,3)= fitresult.A;
train_k(i+72,3) = fitresult.k;
end
Plat3 = A(73,3);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Thermolysin_P1NAsn(:,i),Plat3);
A(i+72,3)= fitresult.A;
train_k(i+72,3) = fitresult.k;
end
title('P1 NAsn')

subplot(2,2,4)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Thermolysin_P3pNAla(:,i));
A(i+72,4)= fitresult.A;
train_k(i+72,4) = fitresult.k;
end
Plat4 = A(73,4);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Thermolysin_P3pNAla(:,i),Plat4);
A(i+72,4)= fitresult.A;
train_k(i+72,4) = fitresult.k;
end
title('P3p NAla')
hold off


%% Chymotrypsin Curve Fitting
ChymotrypsinTraining = readtable('ChymotrypsinTrainingData.xlsx');
Chymotrypsin_Peptide = table2array(ChymotrypsinTraining(:,2:25));
Chymotrypsin_P3NPro = table2array(ChymotrypsinTraining(:,27:50));
Chymotrypsin_P1NAsn = table2array(ChymotrypsinTraining(:,52:75));
Chymotrypsin_P3pNAla = table2array(ChymotrypsinTraining(:,77:100));

figure('Name','Chymotrypsin Training')
subplot(2,2,1)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Chymotrypsin_Peptide(:,i));
A(i+96,1)= fitresult.A;
train_k(i+96,1) = fitresult.k;
end
Plat1 = A(97,1);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Chymotrypsin_Peptide(:,i),Plat1);
A(i+96,1)= fitresult.A;
train_k(i+96,1) = fitresult.k;
end
title('Peptide')

subplot(2,2,2)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Chymotrypsin_P3NPro(:,i));
A(i+96,2)= fitresult.A;
train_k(i+96,2) = fitresult.k;
end
Plat2 = A(97,2);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Chymotrypsin_P3NPro(:,i),Plat2);
A(i+96,2)= fitresult.A;
train_k(i+96,2) = fitresult.k;
end
title('P3 NPro')

subplot(2,2,3)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Chymotrypsin_P1NAsn(:,i));
A(i+96,3)= fitresult.A;
train_k(i+96,3) = fitresult.k;
end
Plat3 = A(97,3);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Chymotrypsin_P1NAsn(:,i),Plat3);
A(i+96,3)= fitresult.A;
train_k(i+96,3) = fitresult.k;
end
title('P1 NAsn')

subplot(2,2,4)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Chymotrypsin_P3pNAla(:,i));
A(i+96,4)= fitresult.A;
train_k(i+96,4) = fitresult.k;
end
Plat4 = A(97,4);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Chymotrypsin_P3pNAla(:,i),Plat4);
A(i+96,4)= fitresult.A;
train_k(i+96,4) = fitresult.k;
end
title('P3p NAla')
hold off


%% Papain Curve Fitting
PapainTraining = readtable('PapainTrainingData.xlsx');
Papain_Peptide = table2array(PapainTraining(:,2:25));
Papain_P3NPro = table2array(PapainTraining(:,27:50));
Papain_P1NAsn = table2array(PapainTraining(:,52:75));
Papain_P3pNAla = table2array(PapainTraining(:,77:100));

figure('Name','Papain Training')
subplot(2,2,1)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Papain_Peptide(:,i));
A(i+120,1)= fitresult.A;
train_k(i+120,1) = fitresult.k;
end
Plat1 = A(121,1);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Papain_Peptide(:,i),Plat1);
A(i+120,1)= fitresult.A;
train_k(i+120,1) = fitresult.k;
end
title('Peptide')

subplot(2,2,2)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Papain_P3NPro(:,i));
A(i+120,2)= fitresult.A;
train_k(i+120,2) = fitresult.k;
end
Plat2 = A(121,2);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Papain_P3NPro(:,i),Plat2);
A(i+120,2)= fitresult.A;
train_k(i+120,2) = fitresult.k;
end
title('P3 NPro')

subplot(2,2,3)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Papain_P1NAsn(:,i));
A(i+120,3)= fitresult.A;
train_k(i+120,3) = fitresult.k;
end
Plat3 = A(121,3);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Papain_P1NAsn(:,i),Plat3);
A(i+120,3)= fitresult.A;
train_k(i+120,3) = fitresult.k;
end
title('P1 NAsn')

subplot(2,2,4)
hold on
for i=[1]
[fitresult,gof] = createFit(time,Papain_P3pNAla(:,i));
A(i+120,4)= fitresult.A;
train_k(i+120,4) = fitresult.k;
end
Plat4 = A(121,4);
for i=[2:24]
[fitresult,gof] = constrainedFit(time,Papain_P3pNAla(:,i),Plat4);
A(i+120,4)= fitresult.A;
train_k(i+120,4) = fitresult.k;
end
title('P3p NAla')
hold off

%% Compile k values

log_train_k = log10(train_k);%Transform exponential constants by taking log10 value
log_train_k_avg = squeeze(mean(reshape(log_train_k,[3,48,4]),1));%Average the three replicates for each value
labels = readtable('AvgLabels.xlsx');%Input labels corresponding to the array
Proteases = table2array(labels(:,1));%Protease names only
Samples = table2array(labels(:,2));%Protease name and concentration
Concentrations = table2array(labels(:,3));%Concentrations only
Substrates = ["Peptide","A NPro","C NAsn","F NAla"];
figure('Name','Training Data')
h = heatmap(log_train_k_avg)%Generate labeled heatmap with each array value
h.YDisplayLabels = Samples
h.XDisplayLabels = Substrates

%% PCA with all training data
[coeff,score,latent,tsquared,explained] = pca(log_train_k_avg);%Built in MATLAB PCA function

figure('Name','Scree Plot w CV')
y= scree(latent,0.000001)%Evaluate variance captured per PC

figure('Name','PCA with Proteases')
PC1=explained(1);
PC2=explained(2);
PC3=explained(3);

p =gscatter(score(:, 1), score(:, 2), Proteases,hsv(6),[],20);%Generate 2D PCA plot colored by protease
xlabel(['PC 1 (' num2str(PC1,3) '%)']);ylabel(['PC 2 (' num2str(PC2,3) '%)']);
xlim([-6 8]);ylim([-4 3])
hold on

figure('Name','PCA 3D')
col=zeros(48,1);
col(1:8)=1;col(9:16)=2;col(17:24)=3;col(25:32)=4;col(33:40)=5;col(41:48)=6;%Assign colors (have to do it manually in 3D)
scatter3(score(:,1),score(:,2),score(:,3),40,col,'filled')%Generate 3D PCA plot
hold on
substrates = {'Peptide', 'A NPro', 'C NAsn', 'F NAla'};
bp = biplot(coeff(:,1:3),'VarLabels', substrates);%Add loading plot
colormap(gca,"hsv(6)")
xlabel(['PC 1 (' num2str(PC1,3) '%)']);ylabel(['PC 2 (' num2str(PC2,3) '%)']);zlabel(['PC 3 (' num2str(PC3,3) '%)']);

%% LDA with all training data
ltk_max = max(log_train_k_avg);
ltk_min = min(log_train_k_avg);
diff = ltk_max-ltk_min;
log_train_k_norm = (log_train_k_avg-ltk_min)./diff;%Normalize training array 
LDAMdl = fitcdiscr(log_train_k_norm,Proteases)%Generate LDA model using all array values with samples labeled by protease identity
[W, LAMBDA] = eig(LDAMdl.BetweenSigma, LDAMdl.Sigma);%Extract eigenvalues for projection 
lambda = diag(LAMBDA);
[lambda, SortOrder] = sort(lambda, 'descend');
F = W(:, SortOrder);

figure('Name','LDA Training Data')
Z = (log_train_k_norm)*F;%Multiple coefficients
p2 =gscatter(Z(:, 1), Z(:, 2), Proteases,hsv(6),[],20);%Generate 2D LDA plot

%Divide scores by protease
score1=Z(1:8,1:2);
x1=score1(:,1); y1=score1(:,2);
score2=Z(8:16,1:2);
x2=score2(:,1); y2=score2(:,2);
score3=Z(17:24,1:2);
x3=score3(:,1); y3=score3(:,2);
score4=Z(25:32,1:2);
x4=score4(:,1); y4=score4(:,2);
score5=Z(33:40,1:2);
x5=score5(:,1); y5=score5(:,2);
score6=Z(41:48,1:2);
x6=score6(:,1); y6=score6(:,2);

%find centers for data from each cluster
meanx1 = mean(x1);
meany1 = mean(y1);
meanx2 = mean(x2);
meany2 = mean(y2);
meanx3 = mean(x3);
meany3 = mean(y3);
meanx4 = mean(x4);
meany4 = mean(y4);
meanx5 = mean(x5);
meany5 = mean(y5);
meanx6 = mean(x6);
meany6 = mean(y6);
%confidence intervals


p=95; 
%find upper and lower bounds for 95% confidence in x and y directions
% CIFcn_x1 = @(x1,p)prctile(x1,abs([0,100]-(100-p)/2));
CIFcn_x1 = prctile(x1,abs([0,100]-(100-p)/2));
CIFcn_y1 = prctile(y1,abs([0,100]-(100-p)/2));
CIFcn_x2 = prctile(x2,abs([0,100]-(100-p)/2));
CIFcn_y2 = prctile(y2,abs([0,100]-(100-p)/2));
CIFcn_x3 = prctile(x3,abs([0,100]-(100-p)/2));
CIFcn_y3 = prctile(y3,abs([0,100]-(100-p)/2));
CIFcn_x4 = prctile(x4,abs([0,100]-(100-p)/2));
CIFcn_y4 = prctile(y4,abs([0,100]-(100-p)/2));
CIFcn_x5 = prctile(x5,abs([0,100]-(100-p)/2));
CIFcn_y5 = prctile(y5,abs([0,100]-(100-p)/2));
CIFcn_x6 = prctile(x6,abs([0,100]-(100-p)/2));
CIFcn_y6 = prctile(y6,abs([0,100]-(100-p)/2));

CI_x1 = abs((CIFcn_x1(1)-CIFcn_x1(2))/2); 
CI_y1 = abs((CIFcn_y1(1)-CIFcn_y1(2))/2); 
CI_x2 = abs((CIFcn_x2(1)-CIFcn_x2(2))/2); 
CI_y2 = abs((CIFcn_y2(1)-CIFcn_y2(2))/2); 
CI_x3 = abs((CIFcn_x3(1)-CIFcn_x3(2))/2); 
CI_y3 = abs((CIFcn_y3(1)-CIFcn_y3(2))/2); 
CI_x4 = abs((CIFcn_x4(1)-CIFcn_x4(2))/2); 
CI_y4 = abs((CIFcn_y4(1)-CIFcn_y4(2))/2); 
CI_x5 = abs((CIFcn_x5(1)-CIFcn_x5(2))/2); 
CI_y5 = abs((CIFcn_y5(1)-CIFcn_y5(2))/2); 
CI_x6 = abs((CIFcn_x6(1)-CIFcn_x6(2))/2); 
CI_y6 = abs((CIFcn_y6(1)-CIFcn_y6(2))/2); 

mean1 = [meanx1,meany1]
std1 = [CI_x1,CI_y1]
mean2 = [meanx2,meany2]
std2 = [CI_x2,CI_y2]
mean3 = [meanx3,meany3]
std3 = [CI_x3,CI_y3]
mean4 = [meanx4,meany4]
std4 = [CI_x4,CI_y4]
mean5 = [meanx5,meany5]
std5 = [CI_x5,CI_y5]
mean6 = [meanx6,meany6]
std6 = [CI_x6,CI_y6]

ellipse_colors = hsv(6);
%add ellipses 
E = plotEllipses(mean1,std1); 
E.FaceColor = [ellipse_colors(1,1:3) .2]; 
E.EdgeColor = ellipse_colors(1,1:3); 
E.LineWidth = 0.5; 
 
E = plotEllipses(mean2,std2); 
E.FaceColor = [ellipse_colors(2,1:3) .2]; 
E.EdgeColor = ellipse_colors(2,1:3); 
E.LineWidth = 0.5; 

E = plotEllipses(mean3,std3); 
E.FaceColor = [ellipse_colors(3,1:3) .2]; 
E.EdgeColor = ellipse_colors(3,1:3); 
E.LineWidth = 0.5; 

E = plotEllipses(mean4,std4); 
E.FaceColor = [ellipse_colors(4,1:3) .2]; 
E.EdgeColor =ellipse_colors(4,1:3); 
E.LineWidth = 0.5; 
 
E = plotEllipses(mean5,std5); 
E.FaceColor = [ellipse_colors(5,1:3) .2]; 
E.EdgeColor = ellipse_colors(5,1:3); 
E.LineWidth = 0.5; 

E = plotEllipses(mean6,std6); 
E.FaceColor = [ellipse_colors(6,1:3) .2]; 
E.EdgeColor = ellipse_colors(6,1:3); 
E.LineWidth = 0.5; 

%% Cross-validated classicfication
cv = cvpartition(Proteases, 'Holdout',0.375)%Partition the data using 0.375 held for testing (3/8 samples). Holdout is stratified by protease identity.
idxTrain = training(cv);
tblTrain = log_train_k_norm(idxTrain,:);%Apply holdout indices to training data
idxNew = test(cv);
tblNew = log_train_k_norm(idxNew,:);%Apply holdout indices to testing data
trainProteases = Proteases(idxTrain,:);%Apply holdout indices to training protease labels
testProteases = Proteases(idxNew,:);%Apply holdout indices to training protease labels
Mdl = fitcdiscr(tblTrain,trainProteases)%Generate a new LDA model using only training datapoints
[W, LAMBDA] = eig(Mdl.BetweenSigma, Mdl.Sigma); 
lambda = diag(LAMBDA);
[lambda, SortOrder] = sort(lambda, 'descend');
G = W(:, SortOrder);
Y = (tblNew)*G;


[LDA_Test,LDA_Test_score,LDA_Test_cost] = predict(Mdl,tblNew);%Use trained model to predict test samples
C = confusionmat(testProteases,LDA_Test);%Construct confusion matrix to visually evaluate classification accuracy
figure('Name','LDA Testing Data Confusion Matrix')
names = {'Collagense','Proteinase K','Elastase','Thermolysin','Chymotrypsin','Papain'};
names = categorical(names);
confusionchart(C,names);
figure('Name','LDA Testing Data Posterior Probabilities')
hmlabels = Mdl.ClassNames
post = heatmap(LDA_Test_score)%Visually evaluate posterior probabilities for each prediction
post.YDisplayLabels = testProteases
post.XDisplayLabels = hmlabels


%% Multiple Regression 
conc = table2array(labels(:,4));%One array of concentrations used
q = length(LDA_Test);
logconc = log10(conc);%Convert concentrataions to log10 values as to fit properly with log10 k constants
%Data remains partitioned for multiple regression. Each of the testing
%samples is fed to the following for loop, which will estimate
%concentration based on the fit of the training samples for the protease it
%was classified as.
for i=[1:q]
if strcmp(LDA_Test(i),'Collagenase')%If protease was identified as collagenase, use this fit
    Collk(:,2:5)= tblTrain(1:5,:);%Select collagenase samples from training table (only 5/8 samples)
    Collk(:,1) = 1;%Constant for fitting
    trainConc = logconc(idxTrain(1:8),:);%Select training concentrations using training indices
    b = regress(trainConc,Collk);%Carry out multiple regression using all four substrates
    g = transpose(tblNew(i,:));
    logfit(i) = b(1)+b(2)*g(1)+b(3)*g(2)+b(4)*g(3)+b(5)*g(4);%Estimate and collect new concentrations based on regression coefficients
end
if strcmp(LDA_Test(i),'ProteinaseK')
    ProtKk(:,2:5)= tblTrain(6:10,:);
    ProtKk(:,1) = 1;
    trainConc = logconc(idxTrain(9:16),:);
    b = regress(trainConc,ProtKk);
    g = transpose(tblNew(i,:));
    logfit(i) = b(1)+b(2)*g(1)+b(3)*g(2)+b(4)*g(3)+b(5)*g(4);
end
if strcmp(LDA_Test(i),'Elastase')
    Elk(:,2:5)= tblTrain(11:15,:);
    Elk(:,1) = 1;
    trainConc = logconc(idxTrain(17:24),:);
    b = regress(trainConc,Elk);
    g = transpose(tblNew(i,:));
    logfit(i) = b(1)+b(2)*g(1)+b(3)*g(2)+b(4)*g(3)+b(5)*g(4);
end
if strcmp(LDA_Test(i),'Thermolysin')
    Thk(:,2:5)= tblTrain(16:20,:);
    Thk(:,1) = 1;
    trainConc = logconc(idxTrain(25:32),:);
    b = regress(trainConc,Thk);
    g = transpose(tblNew(i,:));
    logfit(i) = b(1)+b(2)*g(1)+b(3)*g(2)+b(4)*g(3)+b(5)*g(4);
end
if strcmp(LDA_Test(i),'Chymotrypsin')
    Chk(:,2:5)= tblTrain(21:25,:);
    Chk(:,1) = 1;
    trainConc = logconc(idxTrain(33:40),:);
    b = regress(trainConc,Chk);
    g = transpose(tblNew(i,:));
    logfit(i) = b(1)+b(2)*g(1)+b(3)*g(2)+b(4)*g(3)+b(5)*g(4);
end
if strcmp(LDA_Test(i),'Papain')
    Papk(:,2:5)= tblTrain(26:30,:);
    Papk(:,1) = 1;
    trainConc = logconc(idxTrain(41:48),:);
    b = regress(trainConc,Papk);
    g = transpose(tblNew(i,:));
    logfit(i) = b(1)+b(2)*g(1)+b(3)*g(2)+b(4)*g(3)+b(5)*g(4);
end
end
TestConc = transpose(10.^(logfit));%Convert predicted concentrations from log10 values
SampConc = Concentrations(idxNew,:);%Apply testing indicies to list of concentrations
Concerr = TestConc-SampConc;%Calculate error in sample concentrations
Concpcterr = (abs(Concerr)./SampConc)*100;%Convert error to percentage
Summary = table;%Compile summary table
Summary.Protease = testProteases;
Summary.PredictedProtease = LDA_Test;
Summary.PredictedConcentration = TestConc;
Summary.Error = Concerr;
Summary.PctError = Concpcterr;
Summary
figure('Name','Concentrations with Multiple Regression')
colorarray = hsv(length(unique(LDA_Test)));%Assign each sample a color based on predicted protease identity
for i=[1:q]
if strcmp(LDA_Test(i),'Collagenase')
    clr(i,:) = colorarray(1,:);
elseif strcmp(LDA_Test(i),'ProteinaseK')
    clr(i,:) = colorarray(2,:);
elseif strcmp(LDA_Test(i),'Elastase')
    clr(i,:) = colorarray(3,:);
elseif strcmp(LDA_Test(i),'Thermolysin')
    clr(i,:) = colorarray(4,:);
elseif strcmp(LDA_Test(i),'Chymotrypsin')
    clr(i,:) = colorarray(5,:);
elseif strcmp(LDA_Test(i),'Papain')
    clr(i,:) = colorarray(6,:);
end
end
%Collagenase
subplot(2,3,1)
scatter(SampConc(1:3),TestConc(1:3),[],clr(1:3,:),'filled')%Make plot of actual vs predicted concentration
hold on
title('Collagenase','Color',colorarray(1,:))
xlim([0 80])
xlabel('Actual Concentration (ug/mL)')
ylabel('Predicted Concentration (ug/mL)')
hline = refline(1,0);%Add reference line for perfect concentration determination
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
legend('off')
hold off
%Proteinase K
subplot(2,3,2)
scatter(SampConc(4:6),TestConc(4:6),[],clr(4:6,:),'filled')
title('Proteinase K','Color',colorarray(2,:))
xlim([0 80])
xlabel('Actual Concentration (ug/mL)')
ylabel('Predicted Concentration (ug/mL)')
hold on
hline = refline(1,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
legend('off')
hold off
%Elastase
subplot(2,3,3)
scatter(SampConc(7:9),TestConc(7:9),[],clr(7:9,:),'filled')
title('Elastase','Color',colorarray(3,:))
xlim([0 80])
xlabel('Actual Concentration (ug/mL)')
ylabel('Predicted Concentration (ug/mL)')
hold on
hline = refline(1,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
legend('off')
hold off
%Thermolysin
subplot(2,3,4)
scatter(SampConc(10:12),TestConc(10:12),[],clr(10:12,:),'filled')
title('Thermolysin','Color',colorarray(4,:))
xlim([0 80])
xlabel('Actual Concentration (ug/mL)')
ylabel('Predicted Concentration (ug/mL)')
hold on
hline = refline(1,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
legend('off')
hold off
%Chymotrypsin
subplot(2,3,5)
scatter(SampConc(13:15),TestConc(13:15),[],clr(13:15,:),'filled')
title('Chymotrypsin','Color',colorarray(5,:))
xlim([0 80])
xlabel('Actual Concentration (ug/mL)')
ylabel('Predicted Concentration (ug/mL)')
hold on
hline = refline(1,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
legend('off')
hold off
%Papain
subplot(2,3,6)
scatter(SampConc(16:18),TestConc(16:18),[],clr(16:18,:),'filled')
title('Papain','Color',colorarray(6,:))
xlim([0 80])
xlabel('Actual Concentration (ug/mL)')
ylabel('Predicted Concentration (ug/mL)')
hold on
hline = refline(1,0);
hline.Color = 'k';
hline.LineStyle = '--';
hline.LineWidth = 1;
legend('off')
hold off


%% Exponential Plateau Fit
function [fitresult, gof] = createFit(time, y)
%CREATEFIT(TIME,Y)
%  Create a fit.
%
%  Data for 'Exponential Plateau Fit' fit:
%      X Input : time
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 31-Jan-2022 14:53:57


%% Fit: 'Exponential Plateau Fit'.
[xData, yData] = prepareCurveData( time, y );

% Set up fittype and options.
ft = fittype( 'A-(A)*exp(-k*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [200000 1e-05];
opts.StartPoint = [500000 0.1];
opts.Upper = [700000 Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
%figure( 'Name', 'Exponential Plateau Fit' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. time', 'Exponential Plateau Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Time', 'Interpreter', 'none' );
ylabel( 'Fluorescence (RFU)', 'Interpreter', 'none' );
grid on
end


%% Constrained Plateau Fit
function [fitresult, gof] = constrainedFit(time, y,Plat)
%CREATEFIT(TIME,Y)
%  Create a fit.
%
%  Data for 'Exponential Plateau Fit' fit:
%      X Input : time
%      Y Output: y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 31-Jan-2022 14:53:57


%% Fit: 'Constrained Plateau Fit'.
[xData, yData] = prepareCurveData( time, y);

% Set up fittype and options.
ft = fittype( 'A-(A)*exp(-k*x)', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [Plat 1e-05];
opts.StartPoint = [Plat 0.1];
opts.Upper = [Plat Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
%figure( 'Name', 'Exponential Plateau Fit' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. time', 'Exponential Plateau Fit', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Time', 'Interpreter', 'none' );
ylabel( 'Fluorescence (RFU)', 'Interpreter', 'none' );
grid on
end


%% Scree Plot Function
function y = scree( latent, alpha );
newplot

p = 100*latent/sum(latent);
pe = cumsum(p);

if ( nargin < 2 )
    alpha = 0.05;
end;

i = find(pe > 100*(1 - alpha),1);
if ( isempty(i) ), i = length(latent); end;

if nargout > 0
    y = [pe(1:i) p(1:i)];
end

line(1:i, pe(1:i),'marker', 'o', 'color', 'b', 'markerfacecolor', 'b' );
line(1:i, p(1:i),'marker', 'o', 'color', 'k', 'markerfacecolor', 'k' ); 


xlabel('Principle Component Number','FontWeight', 'bold' );
ylabel('Percent Variance Captured','FontWeight', 'bold' );
legend( {'Cumulative', 'Individual'}, 'location', 'northwest' );

title( 'Scree Plot with Percent Variance');
end

function h = plotEllipses(cnt,rads,axh)
% cnt is the [x,y] coordinate of the center (row or column index).
% rads is the [horizontal, vertical] "radius" of the ellipses (row or column index).
% axh is the axis handle (if missing or empty, gca will be used)
% h is the object handle to the plotted rectangle.
% The easiest approach IMO is to plot a rectangle with rounded edges. 
% EXAMPLE
%    center = [1, 2];         %[x,y] center (mean)
%    stdev = [1.2, 0.5];      %[x,y] standard dev.
%    h = plotEllipses(center, stdev)
%    axis equal
% get axis handle if needed
if nargin < 3 || isempty(axh)
   axh = gca();  
end
% Compute the lower, left corner of the rectangle.
llc = cnt(:)-rads(:);
% Compute the width and height
wh = rads(:)*2; 
% Draw rectangle 
h = rectangle(axh,'Position',[llc(:).',wh(:).'],'Curvature',[1,1]); 
end
function y = scree( latent, alpha )
newplot

p = 100*latent/sum(latent);
pe = cumsum(p);

if ( nargin < 2 )
    alpha = 0.05;
end;

i = find(pe > 100*(1 - alpha),1);
if ( isempty(i) ), i = length(latent); end;

if nargout > 0
    y = [pe(1:i) p(1:i)];
end

line(1:i, pe(1:i),'marker', 'o', 'color', 'b', 'markerfacecolor', 'b' );
line(1:i, p(1:i),'marker', 'o', 'color', 'k', 'markerfacecolor', 'k' ); 


xlabel('Principle Component Number','FontWeight', 'bold' );
ylabel('Percent Variance Captured','FontWeight', 'bold' );
legend( {'Cumulative', 'Individual'}, 'location', 'northwest' );

title( 'Scree Plot with Percent Variance');