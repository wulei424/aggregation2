
frames=2:400;
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
%     figure(1);imshowpair(I,bw);
%     title(['im = ',num2str(i)])
%     pause;
    
    %A2(i)=mean(bw(:));
    
    % bootstrapping error estimate
    % let's split the image up into blocks to get a sense for the error
    ncells=32;
    Asplit = mat2cell( bw, ncells*ones(floor(size(bw,1)/ncells),1), ncells*ones(floor(size(bw,2)/ncells),1));
    meanValues=cellfun(@(x) mean(x(:)), Asplit);
    sigma_A(i)=2*std(meanValues(:))./sqrt(length(meanValues));
    A(i)=mean(meanValues(:));
end

figure;errorbar(frames,A,sigma_A)


%figure;plot(A);