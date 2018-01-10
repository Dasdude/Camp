close all
addpath('../Dataset');
dataset  = 'same_density_0to70.csv';
if ~exist('data','var')==1
    data  = readtable(dataset);
end
densities = 0:5:50
for j= 2:length(densities)
row_index = (data.Av_Density<densities(j));
row_index2 = data.Av_Density>=densities(j-1);
row_index = logical(row_index.*row_index2);
if ~any(row_index)
    continue
end
variable_names = {'Av_Density','RSS','Range'};
data_density_sep = data(row_index,variable_names);
dd = data_density_sep;



%% 
range_quantization_par = 1;
range_iter = 0:range_quantization_par:400;
hist_range = -100:-30;
dd.Range = floor(dd.Range/range_quantization_par)*range_quantization_par;
[range_mesh,bin_mesh] = meshgrid(range_iter,hist_range);
pdf_concat=[]
cmap_concat = []

for i=range_iter
    index = dd.Range==i;
    ddr = dd(index,variable_names);
    a = size(ddr);

        [pdf,edge]=hist(ddr.RSS,-100:-30);
        pdf = pdf/sum(pdf);
        cmap_temp  = pdf/max(pdf(:));
        pdf_concat = [pdf_concat,pdf'];
        cmap_concat = [cmap_concat,cmap_temp'];

%         title(['Range:' ,num2str(i)])
%         saveas(gcf,['../Matlab/Plots/',dataset,'Range:' ,num2str(i),'Fading Distribution.png'])
%         pause

end

figure
s=surf(range_mesh,bin_mesh,pdf_concat,cmap_concat,'FaceAlpha',1);
% s=surf(range_mesh,bin_mesh,pdf_concat,'FaceAlpha',1);
s.EdgeColor='none'
colorbar
title(['Density',num2str(densities(j-1)),':',num2str(densities(j))])
end
% colormap jet