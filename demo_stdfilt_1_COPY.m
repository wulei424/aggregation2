
clear

frames=2:100:400;
stats={};
for i=1:length(frames);
    frame=frames(i);
    I = imread(sprintf('1204-3-4_t%d.TIF',frame));
    Ib = imread(sprintf('1204-3-4_background.tif'));
    I3=I-mean(I(:));
    I2=(stdfilt(I,ones(7,7)));
    
    Inorm_b=( (Ib-min(Ib(:)))/(max(Ib(:))-min(Ib(:))));
    th_b=graythresh( Inorm_b );
    bw_b=im2bw(Inorm_b,th_b/4);
    
    Inorm_2=( (I2-min(I2(:)))/(max(I2(:))-min(I2(:))));
    th_2=graythresh( Inorm_2 );
    bw_2=im2bw(Inorm_2,th_2/4);
    
    Inorm_3=( (I3-min(I3(:)))/(max(I3(:))-min(I3(:))));
    th_3=graythresh( Inorm_3 );
    bw_3=im2bw(Inorm_3,th_3/4);
    
%     maybe try to undo the broadening due to doing filter over a window    
%     se=strel('disk',3);
%     bw2=imerode(bw,se);

%     to visualize the individual areas
    % figure(1);imshowpair(I,bw);
    % title(['im = ',num2str(i)])
     % pause;
    
    %A2(i)=mean(bw(:));
    
    % bootstrapping error estimate
    % let's split the image up into blocks to get a sense for the error
 
    A_b(i)=mean(bw_b(:));
    
    ncells=32;
    Asplit_2 = mat2cell( bw_2, ncells*ones(floor(size(bw_2,1)/ncells),1), ncells*ones(floor(size(bw_2,2)/ncells),1));
    meanValues_2=cellfun(@(x) mean(x(:)), Asplit_2);
    sigma_A_2(i)=2*std(meanValues_2(:))./sqrt(length(meanValues_2));
    A_2(i)=mean(meanValues_2(:));
    
    Asplit_3 = mat2cell( bw_3, ncells*ones(floor(size(bw_3,1)/ncells),1), ncells*ones(floor(size(bw_3,2)/ncells),1));
    meanValues_3=cellfun(@(x) mean(x(:)), Asplit_3);
    sigma_A_3(i)=2*std(meanValues_3(:))./sqrt(length(meanValues_3));
    A_3(i)=mean(meanValues_3(:));
    
     [Afromcorr_2, radial_average_2, radii_list_2, Cnorm_2]=meanareafromautocorr(I2,250);
     [Afromcorr_3, radial_average_3, radii_list_3, Cnorm_3]=meanareafromautocorr(I3,250);
     [Afromcorr_b, radial_average_b, radii_list_b, Cnorm_b]=meanareafromautocorr(Ib,250);
    %figure(1);plot(radii_list* .180^2,radial_average/Cnorm,'-','LineWidth',3,'Color',mycolors(i,:));xlabel('Radii [\mum]');ylabel('normalized autocorrelation');ylim([0,1]);hold all;
    
    stats{i}.surfacecoverage_b=A_b(i);
    %stats{i}.sigma_surfacecoverage_b=sigma_A_b(i);
    stats{i}.Acorr_b=Afromcorr_b;
    stats{i}.image_autocorr_vals_b=radial_average_b;
    stats{i}.image_autocorr_rad_b=radii_list_b;
    stats{i}.image_autocorr_norm_b=Cnorm_b;
    
    stats{i}.surfacecoverage_2=A_2(i);
    stats{i}.sigma_surfacecoverage_2=sigma_A_2(i);
    stats{i}.Acorr_2=Afromcorr_2;
    stats{i}.image_autocorr_vals_2=radial_average_2;
    stats{i}.image_autocorr_rad_2=radii_list_2;
    stats{i}.image_autocorr_norm_2=Cnorm_2;
    
    stats{i}.surfacecoverage_3=A_3(i);
    stats{i}.sigma_surfacecoverage_3=sigma_A_3(i);
    stats{i}.Acorr_3=Afromcorr_3;
    stats{i}.image_autocorr_vals_3=radial_average_3;
    stats{i}.image_autocorr_rad_3=radii_list_3;
    stats{i}.image_autocorr_norm_3=Cnorm_3;
    
    
end

save('rad_stats.mat','stats')
%figure(1);
%subplot (2,1,1);
%errorbar(frames/10,A_2,sigma_A_2,'b');xlabel('Time [min]');ylabel('surfacecoverage' );hold on;
%subplot (2,1,2);
%errorbar(frames/10,A_2,sigma_A_2,'r');xlabel('Time [min]');ylabel('surfacecoverage');hold on;
%figure(2);plot(A_2); hold on;
%figure (3);plot (A_3);hold on;


%%
mycolors=jet(length(frames));
for i=1:length(frames);
    Afromcorr_b= stats{i}.Acorr_b;
    radial_average_b=stats{i}.image_autocorr_vals_b;
    radii_list_b=stats{i}.image_autocorr_rad_b;
    Cnorm_b=stats{i}.image_autocorr_norm_b;
    
    Afromcorr_2= stats{i}.Acorr_2;
    radial_average_2=stats{i}.image_autocorr_vals_2;
    radii_list_2=stats{i}.image_autocorr_rad_2;
    Cnorm_2=stats{i}.image_autocorr_norm_2;
    
    Afromcorr_3= stats{i}.Acorr_3;
    radial_average_3=stats{i}.image_autocorr_vals_3;
    radii_list_3=stats{i}.image_autocorr_rad_3;
    Cnorm_3=stats{i}.image_autocorr_norm_3;

    Cth=0.2;
    r_2=find(radial_average_2./Cnorm_2 < Cth, 1, 'First');
    r_3=find(radial_average_3./Cnorm_3 < Cth, 1, 'First');
    r_b=find(radial_average_b./Cnorm_b < Cth, 1, 'First');
    
      if isempty(r_b)
        warning('correlation did not drop below threshold within sampling window. make maxlag larger or get larger images');
        A_b= pi* radii_list_b(end)^2;
    else
        A_b= pi* radii_list_b(r_b)^2;
    end
    
    if isempty(r_2)
        warning('correlation did not drop below threshold within sampling window. make maxlag larger or get larger images');
        A_2= pi* radii_list_2(end)^2;
    else
        A_2= pi* radii_list_2(r_2)^2;
    end
    
    if isempty(r_3)
        warning('correlation did not drop below threshold within sampling window. make maxlag larger or get larger images');
        A_3= pi* radii_list_3(end)^2;
    else
        A_3= pi* radii_list_3(r_3)^2;
    end
    
   % Acorr2(i)=sum(radial_average./Cnorm);
   Acorr2_b(i)=A_b;
   A_blob_b=pi.*radii_list_b.^2 * .180^2;
   
   Acorr2_2(i)=A_2;
   A_blob_2=pi.*radii_list_2.^2 * .180^2;
   
   Acorr2_3(i)=A_3;
   A_blob_3=pi.*radii_list_3.^2 * .180^2;
   A_blob_2_f=A_blob_2-A_blob_b*(512*512)/(70*70);
   A_blob_3_f=A_blob_3-A_blob_b*(512*512)/(70*70);
   radial_average_2_f=radial_average_2-radial_average_b*(512*512)/(70*70);
   radial_average_3_f=radial_average_3-radial_average_b*(512*512)/(70*70);
   Cnorm_2_f=Cnorm_2-Cnorm_b*(512*512)/(70*70);
   Cnorm_3_f=Cnorm_3-Cnorm_b*(512*512)/(70*70);
  figure(1);semilogx(A_blob_2,radial_average_2_f/Cnorm_2_f,'-',A_blob_3,radial_average_3_f/Cnorm_3_f,':','LineWidth',2,'Color',mycolors(i,:)); xlabel('Blob area [\mum^2]'); ylabel('Normalized autocorrelation'); ylim([0,1]);hold all;
  hold on;

    
end


%% 
   