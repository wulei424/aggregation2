
frames=2:1:400;
stats={};
for i=1:length(frames);
    frame=frames(i);
    I = imread(sprintf('1204-3-4_t%d.TIF',frame));
    I2=(stdfilt(I,ones(7,7)));
    Inorm=( (I2-min(I2(:)))/(max(I2(:))-min(I2(:))));
    th=graythresh( Inorm );
    bw=im2bw(Inorm,th/4);
    
%     maybe try to undo the broadening due to doing filter over a window    
%     se=strel('disk',3);
%     bw2=imerode(bw,se);

%     to visualize the individual areas
     figure(1);imshowpair(I,bw);
     title(['im = ',num2str(i)])
     pause;
    
    %A2(i)=mean(bw(:));
    
    % bootstrapping error estimate
    % let's split the image up into blocks to get a sense for the error
    ncells=32;
    Asplit = mat2cell( bw, ncells*ones(floor(size(bw,1)/ncells),1), ncells*ones(floor(size(bw,2)/ncells),1));
    meanValues=cellfun(@(x) mean(x(:)), Asplit);
    sigma_A(i)=2*std(meanValues(:))./sqrt(length(meanValues));
    A(i)=mean(meanValues(:));
    %[Afromcorr, radial_average, radii_list, Cnorm]=meanareafromautocorr(I2,250);
    %figure(1);plot(radii_list* .180^2,radial_average/Cnorm,'-','LineWidth',3,'Color',mycolors(i,:));xlabel('Radii [\mum]');ylabel('normalized autocorrelation');ylim([0,1]);hold all;
    
    stats{i}.surfacecoverage=A(i);
    stats{i}.sigma_surfacecoverage=sigma_A(i);
    stats{i}.Acorr=Afromcorr;
    stats{i}.image_autocorr_vals=radial_average;
    stats{i}.image_autocorr_rad=radii_list;
    stats{i}.image_autocorr_norm=Cnorm;
    
    
end

save('rad_stats.mat','stats')

figure(1);errorbar(frames/10,A,sigma_A);
xlabel('Time [min]');ylabel('surfacecoverage [\mum^2]');hold on;
figure(2);plot(A);


%%
mycolors=jet(length(frames));
for i=1:length(frames);
    Afromcorr= stats{i}.Acorr;
    radial_average=stats{i}.image_autocorr_vals;
    radii_list=stats{i}.image_autocorr_rad;
    Cnorm=stats{i}.image_autocorr_norm;
    Cth=0.2;
    r=find(radial_average./Cnorm < Cth, 1, 'First');
    
    if isempty(r)
        warning('correlation did not drop below threshold within sampling window. make maxlag larger or get larger images');
        A= pi* radii_list(end)^2;
    else
        A= pi* radii_list(r)^2;
    end
    
   % Acorr2(i)=sum(radial_average./Cnorm);
   Acorr2(i)=A;
   A_blob=pi.*radii_list.^2 * .180^2;
   
    figure(1);semilogx(A_blob,radial_average/Cnorm,'-','LineWidth',3,'Color',mycolors(i,:)); xlabel('Blob area [\mum^2]'); ylabel('Normalized autocorrelation'); ylim([0,1]);hold all;
    hold on;
    
    
end
   figure(2);histfit(A_blob,30);xlabel('Blob area [\mum^2]');ylabel('PDF');
   hold on;
   figure(3);plot(frames,Acorr2*(0.180)^2);hold on;
   
