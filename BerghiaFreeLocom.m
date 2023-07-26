%Import trajectories and save in cell
%Import inf and relin columns of spreadsheet

for i=1:size(data,2)
    clear coor init maxshort inf relin P* V* lP* DP
    coor=data{i};
%     figure;plot(coor(:,1),coor(:,2))
%     saveas(gcf,sprintf('Traj%d.tif',i));

   init=idxs(i,1);
   maxshort=idxs(i,3);
   inf=idxs(i,5);
   relin=idxs(i,6);


P1=coor(relin-100,:);
P2=coor(relin,:);
if size(coor,1) <= relin+100
    P3=coor(end,:);
else
P3=coor(relin+100,:);
end

V12(1)=P2(1)-P1(1);
V12(2)=P2(2)-P1(2);
V23(1)=P3(1)-P2(1);
V23(2)=P3(2)-P2(2);

DP=-V12(1)*V23(1)-V12(2)*V23(2);

lP12=sqrt((P2(2)-P1(2))^2 + (P2(1)-P1(1))^2);
lP23=sqrt((P3(2)-P2(2))^2 + (P3(1)-P2(1))^2);

angle(i)=acosd((DP)/(lP12*lP23));
end

angle=180-angle;
angle=angle';

%Calculate and Bin Speeds in mm/s
tstep=1/30;
pixeltomm=0.110966057;
for i=1:size(data,2)
maxsizes(i)=size(data{i},1);
end
maxsize=max(maxsizes);

dec_fac=15;
spd=zeros(size(data,2),maxsize);
spd_dec=zeros(size(data,2),floor(maxsize/dec_fac));


for i=1:size(data,2)
    clear coor dist coor_dec dist_dec
    coor=data{i};
    for t=1:size(coor,1)-1
        dist(t)=pixeltomm*sqrt((coor(t+1,1)-coor(t,1))^2 + (coor(t+1,2)-coor(t,2))^2);
    end
    spd(i,1:length(dist))=dist/tstep;

    coor_dec(:,:)=coor(1:dec_fac:end,:);
    for t=1:size(coor_dec,1)-1
        dist_dec(t)=pixeltomm*sqrt((coor_dec(t+1,1)-coor_dec(t,1))^2 + (coor_dec(t+1,2)-coor_dec(t,2))^2);
    end
    spd_dec(i,1:length(dist_dec))=dist_dec/(tstep*dec_fac);
end

clearvars -except angle data idxs pixeltomm spd spd_dec tstep
save 070921_An13.mat

%Plot Speed Histograms
for i=1:size(data,2)
    figure;bar(spd(i,:))
end

for i=1:size(hspd_dec_1A,1)
    figure;bar(hspd_dec_1A(i,:))
end



% %Plot XYT trajectory
% figure; plot3(coor(:,1),coor(:,2),0:tstep:tstep*length(dist)); %Plots positions (time on z-axis)
% % hold on
% % [X,Y]=meshgrid(min(coor(:,1)):max(coor(:,1)),min(coor(:,2)):max(coor(:,2)));
% % Z=zeros(size(X,1),size(X,2));
% % Z(:)=20;
% % surf(X,Y,Z,'FaceAlpha',0.01);
% xlabel('x coordinate')
% ylabel('y coordinate')
% zlabel('Time (s)')
% set(gca,'FontSize',10)
% % title('XYT Trajectory of Animal in 50 µM #2 5-HT (1-3-20), t_{tot}=1 hr.')

%Center all trajectories with one another (-5 s before stim, +20 s after
%stim)

tspd_dec_1A(isnan(tspd_dec_1A))=0;
tail_1A_cen=zeros(size(tspd_dec_1A,1),51);
for i=1:size(tspd_dec_1A,1)
    if tidxs_1A(i,7)-10 < 1
        q=11-tidxs_1A(i,7);
        tail_1A_cen(i,1:q)=0;
        tail_1A_cen(i,(q+1):end)=tspd_dec_1A(i,((tidxs_1A(i,7)-10)+q):(tidxs_1A(i,7)+40));
      
    else
        tail_1A_cen(i,:)=tspd_dec_1A(i,(tidxs_1A(i,7)-10):(tidxs_1A(i,7)+40));
    end
end

%Attenuate any (artifactual) speed bars > median + 4*MAD
for i=1:size(tail_1A_cen,1)
    median=median(tail_1A_cen(i,:));
    mad=mad(tail_1A_cen(i,:),1);
tail_1A_cen(i,tail_1A_cen(i,:)>median+4*mad)=median+4*mad;
clear median mad
end
  
%Calcuate mean and STE for every animal (tail and tail) across all time
%steps
for i=1:size(tail_1A_cen,2)
tail_1A_mean(i)=mean(nonzeros(tail_1A_cen(:,i)));
tail_1A_ste(i)=std(nonzeros(tail_1A_cen(:,i)))/sqrt(size(nonzeros(tail_1A_cen(:,i)),1));
end

% 
% for i=1:size(tail_1A_cen,1)
%     figure;bar(tail_1A_cen(i,:))
% end

figure;bar(tail_1A_mean)