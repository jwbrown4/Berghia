%%Read in csv file and partition data per time index
clear all
name=input('Please type the name of the FIJI file to be read in, including the .csv suffix: \n','s');
M=csvread(name,1,0);
M=[M(:,4) M(:,2) M(:,3)];
time_ind=unique(M(:,1),'stable');

for i=1:length(time_ind)
    q=find(M(:,1)==time_ind(i));
    points{i}=M(q,2:3);
end

% figure; scatter(points{1}(:,1),points{1}(:,2))
% hold on; scatter(points{2}(:,1),points{2}(:,2),'r')

for j=1:size(points,2)
    clear temp_points D mins D_half D_all
temp_points=points{j};
D=squareform(pdist(temp_points));
D_half=triu(D);
D_all=nonzeros(D_half);
D(D==0)=NaN;

% for i=1:size(D,2)
%     D_sort{i}=sort(D(i+1:size(D,1),i));
% end

D_sort=sort(D,2);
D_3min=D_sort(:,1:3);


for i=1:size(D,1)
    mins(i)=min(D(:,i));
end

avg_min(j)=mean(mins);%mean(unique(mins));
avg_3min(j)=mean(unique(D_3min(:)));
avg_tot(j)=mean(D_all);
med_min(j)=median(unique(mins));
med_3min(j)=median(unique(D_3min(:)));
med_tot(j)=median(D_all);

end

avg_min=reshape(avg_min,2,[])';
avg_3min=reshape(avg_3min,2,[])';
avg_tot=reshape(avg_tot,2,[])';
med_min=reshape(med_min,2,[])';
med_3min=reshape(med_3min,2,[])';
med_tot=reshape(med_tot,2,[])';

output=[avg_min avg_3min avg_tot med_min med_3min med_tot];
% 
name2=input('Please type the desired name of your processed csv output file, including the .csv suffix: \n','s');
csvwrite(name2,output)