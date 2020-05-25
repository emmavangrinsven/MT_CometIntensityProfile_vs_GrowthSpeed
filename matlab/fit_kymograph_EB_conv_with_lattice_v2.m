clear all


%% parameters
%fitting window 
%halfwidth=30;
halfwidthright=15;
halfwidthleft=35;
bMakeMovie=false;

nUmPerPixel=0.065;
%nMinPerFrame=0.5/60;
nMinPerFrame=0.1/60;


%filebase='SLAIN2 100NM A_K002';
[ch1_filenames,folder_name] = load_all_filenames();

%going through all profiles
nTotFiles=numel(ch1_filenames);
nSumOfFit=[];
SumOfProfiles={};
nProfilesTotal=0;
for nFile=1:nTotFiles
    filebase = strcat(folder_name,'\',ch1_filenames{nFile});
    sWork=strcat('working on (',num2str(nFile),'/',num2str(nTotFiles),') ',ch1_filenames{nFile});
    disp(sWork);

    wb = waitbar(0,sWork);


    %read tiff file
    filenameG=strcat(filebase,'.tif');

    I_G=double(imread(filenameG)');

    szG=size(I_G);

    imwidth=szG(1);
    imheight=szG(2);

    %read Kymo_ROI
    sROI = ReadImageJROI(strcat(filebase,'.roi'));
    if(strcmp(sROI.strType,'Line'))
       sROI.mnCoordinates=zeros(2,2);
       sROI.mnCoordinates(1,:)=sROI.vnLinePoints(1:2);
       sROI.mnCoordinates(2,:)=sROI.vnLinePoints(3:4);
    end
    sROI.mnCoordinates=round(sROI.mnCoordinates);
    nPoints=size(sROI.mnCoordinates);
    nPoints=nPoints(1);

    nTotProfiles=sROI.mnCoordinates(end,2)-sROI.mnCoordinates(1,2);

    %fitting functions: exponential convolved with gaussian of specified width
    erfampG = @(bkg, amp1, amp2, mu, lambda, x) ...
            bkg + 0.5*amp1*exp((1.72*1.72*lambda*lambda*0.5)+lambda.*x-lambda*mu).*(1-erf((lambda*1.72*1.72+x-mu)/(sqrt(2)*1.72)))+...
            0.5*amp2*(erfc((x-mu)/(1.72*sqrt(2))));


    clear 'MovieProfG';
    MovieProfG(nTotProfiles) = struct('cdata',[],'colormap',[]);

    %arrays for storage of fitting results
    nFitParamNG=5;

    expfitresultG = zeros(nTotProfiles,nFitParamNG+4);
    xinterval=zeros(nTotProfiles,2);
    iterNG=zeros(nTotProfiles,1);
    fitoutG=cell(nTotProfiles,1);

    i=1;

    for(nPoint=2:nPoints)

        %line along comet
        x1=sROI.mnCoordinates(nPoint-1,1)+1;
        y1=sROI.mnCoordinates(nPoint-1,2)+1;
        x2=sROI.mnCoordinates(nPoint,1)+1;
        y2=sROI.mnCoordinates(nPoint,2)+1;

        invslope=(x2-x1)/(y2-y1);



        %figure

        for(ycurr=(y1+1):y2)

            %update waitbar
            waitbar(i/ nTotProfiles);


            xbeg=round(invslope*(ycurr-y1)+x1-halfwidthleft);
            xend=round(invslope*(ycurr-y1)+x1+halfwidthright);
            if(xend>imwidth)
                xend=imwidth;
            end
            if(xbeg<1)
                xbeg=1;
            end

            %store fitting range
            xinterval(i,:)=[xbeg xend];
            %% GREEN CHANNEL FITTING
            fitx = (xbeg:xend)';    
            fity = I_G(xbeg:xend,ycurr);
            [maxy,maxind]=max(fity);
            miny=min(fity);
            maxamprange=maxy-miny;
            %taking care of initial conditions
            ini=zeros(nFitParamNG,1);
            ini(1) = miny;
            ini(2) = 2*maxamprange;
            ini(3) = 0.2*maxamprange;           
            %ini(4)=xbeg+0.5*(xend-xbeg);
            ini(4)=round(invslope*(ycurr-y1)+x1);
            %ini(4)=fitx(maxind);
            ini(5) = 1.0;

            %debug stuff
            if(ycurr==280)
                ycurr=ycurr+1;
                ycurr=ycurr-1;
            end
    % {
                options = fitoptions('Method','NonlinearLeastSquares','StartPoint', ini,...
                    'Lower',[miny 0 0 xbeg 0],...
                    'Upper',[maxy 10*maxamprange maxamprange xend 20],...
                    'MaxFunEvals',1000,'MaxIter',1000);
            % }
            bFitSuccess=true;
            try
                [fitt,gof,output] = fit(fitx,fity,erfampG, options);
            catch
                disp('NaN warning');
                bFitSuccess=false;
            end

            iterNG(i)=output.iterations;
            fitoutG{i,1}=output.message;
            if(bFitSuccess)
                coeff = coeffvalues(fitt);
                expfitresultG(i,:) =horzcat(coeff,ycurr,0,0,nFile);
                fittedy = erfampG(coeff(1),coeff(2),coeff(3),coeff(4),coeff(5),fitx);
                %store profile
                fittedMax=max(fittedy);
                %normalize intensity
                normfity=(fity-coeff(1))/(fittedMax-coeff(1));
                shiftedfitx=fitx-coeff(4);
                nProfilesTotal=nProfilesTotal+1;
                SumOfProfiles{nProfilesTotal,1}=horzcat(shiftedfitx,normfity);
                
            else
                coeff = ini';
                expfitresultG(i,:) =horzcat(NaN(1,nFitParamNG),ycurr,0,0,nFile);
            end
            if(bMakeMovie)                
                h=figure;
                plot(fitx,fity);
                hold on
                plot(fitx,fittedy,'--');
                %plot(fitx,fittedgauss,'--g');
                %line([coeff(4) coeff(4)],[coeff(1) coeff(2)+coeff(1)],'Color','b');
                %line([coeff(4)+coeff(6) coeff(4)+coeff(6)],[coeff(1) coeff(2)+coeff(1)],'Color','g');
                hold off
                MovieProfG(i) = getframe(gcf);
                close(h)
            end



            i=i+1;
        end
        
end
    close(wb) 
    speedgrowth=diff(expfitresultG(:,4));
    expfitresultG(1:(end-1),nFitParamNG+3)=speedgrowth;
    expfitresultG(:,nFitParamNG+2)=nUmPerPixel*(1./expfitresultG(:,5));
    
    %outliers removal
    cometposition = expfitresultG(:,4);
    cometsliding = smooth(cometposition,5);
    stdval=median(abs(cometposition-cometsliding));
    diff_x=abs(cometposition-cometsliding);
    

    %% MTA analysis 
	xvals=expfitresultG(:,6)*nMinPerFrame;
	yvals=smooth(expfitresultG(:,4),10)*nUmPerPixel;
    %yvals=expfitresultG(:,4)*nUmPerPixel;
    %yvals=expfitresultG(:,4);
	xyArr=horzcat(xvals,yvals);
	[optimal_epoches, slopes, xyApprox] = mta_analysis(xyArr,2);
 
	%xvals=expfitresultG(:,6);
	yvals=expfitresultG(:,4)*nUmPerPixel;

    %%fill average speed values
    MTASpeed=zeros(length(xvals),1);
    nEpN=length(slopes);
    for(ind=1:nEpN)
        MTASpeed(optimal_epoches(ind):optimal_epoches(ind+1),1)=slopes(ind);
    end
    %average slope number in the transition (epoch) points
    if(nEpN>2)
        for(ind=2:nEpN)
            MTASpeed(optimal_epoches(ind),1)=mean([slopes(ind-1) slopes(ind)]);
        end
    end
    figure
    plot(xvals,yvals,'o')    
    hold on
    plot(xyApprox(:,1),xyApprox(:,2),'-r')  
    expfitresultG(:,nFitParamNG+3)=MTASpeed;
    nSumOfFit=vertcat(nSumOfFit,expfitresultG);

end

nComLength_Speed=horzcat(nSumOfFit(:,7),nSumOfFit(:,8));
nComLengthMSDN=horzcat(mean(nComLength_Speed(:,1)),std(nComLength_Speed(:,1)),numel(nComLength_Speed(:,1)));
nSpeedMSDN=horzcat(mean(nComLength_Speed(:,2)),std(nComLength_Speed(:,2)),numel(nComLength_Speed(:,2)));
nSpeed=zeros(nTotFiles,2);
nCom=zeros(nTotFiles,2);
for i=1:nTotFiles
    filt=nSumOfFit(:,9)==i;
    nCom(i,1)=mean(nComLength_Speed(filt,1));
    nCom(i,2)=std(nComLength_Speed(filt,1))/sqrt(numel(nComLength_Speed(filt,1)));
    nSpeed(i,1)=mean(nComLength_Speed(filt,2));
    nSpeed(i,2)=std(nComLength_Speed(filt,2))/sqrt(numel(nComLength_Speed(filt,2)));

end

% 
% expfitresultG(:,nFitParamNG+2)=(expfitresultG(:,3)*(1.5*(2*pi)^0.5))./expfitresultG(:,2);
% speedgrowth=diff(expfitresultG(:,4));
% expfitresultG(1:(end-1),nFitParamNG+3)=speedgrowth;
% 
% %expfitresultG(:,nFitParamN+2)=(expfitresultG(:,3)*(1.5*(2*pi)^0.5))./expfitresultG(:,2);
% speedgrowth=diff(expfitresultR(:,3));
% expfitresultR(1:(end-1),nFitParamNR+2)=speedgrowth;
% 
% %close waitbar
% close(wb) 
% %{
%  R = corrcoef(expfitresultG(:,(nFitParamN+2):(nFitParamN+3)));
%  Rminus=R(1,2);
%  R = corrcoef(horzcat(expfitresultG(2:end,(nFitParamN+2)),expfitresult(1:(end-1),(nFitParamN+3))));
%  Rplus=R(1,2);
%  %}
% 
% if(bShowPlot)
%     maxIm=max(max(I_G));
%     minIm=min(min(I_G));
%     I2=(I_G-minIm)/(maxIm-minIm);
%     handlefig=figure;
%     imshow(I2);
%     hold on
%     
%     %position of error function
%     plot(expfitresultG(:,7),expfitresultG(:,4),'Color','b');
% 
%     t= 0:pi/10:2*pi;
% 
%     patchR=0.5;
%     alphavals=expfitresultG(:,8);
%     maxAlpha = quantile(alphavals,0.80);
%     alphavals=alphavals/maxAlpha;
%     morex=alphavals>=1.0;
%     alphavals(morex)=1;
%     for i=1:size(expfitresultG(:,7))
%         %ddy1=expfitresultR(i,3)+expfitresultR(i,4);
%         %ddy2=expfitresultR(i,3)-expfitresultR(i,4);
%         % {
%         ddy1=expfitresultG(i,4);
%         %if(bFitGaussianExp)
% %            ddy2=expfitresultR(i,3)-1/expfitresultR(i,nFitParamNR);
%  %       else
%             ddy2=expfitresultG(i,4)+expfitresultG(i,5);
%   %      end
% 
%         line([expfitresultG(i,7) expfitresultG(i,7)],[ddy1 ddy2],'Color','b');
% % }
%         %accumulated intensity ratio
%         x=expfitresultG(i,nFitParamNG+1);
%         y=expfitresultG(i,4)+expfitresultG(i,6); 
%         pb=patch((patchR*sin(t)+ x),(patchR*cos(t)+y),'g','edgecolor','none','FaceAlpha',alphavals(i));
%     end
% 
%     
%     %plot(expfitresultG(:,nFitParamNG+1),expfitresultG(:,4),'Color','g');
%     plot(expfitresultR(:,7),expfitresultR(:,4),'Color','r');
%     saveas(handlefig,strcat(filebase,'_G.fig'),'fig');
%     
%      %% EB channel
%     maxIm=max(max(I_R));
%     minIm=min(min(I_R));
%     I2=(I_R-minIm)/(maxIm-minIm);
%     handlefig=figure;
%     imshow(I2);
%     hold on
%     plot(expfitresultR(:,7),expfitresultR(:,4),'Color','r');
% 
%     t= 0:pi/10:2*pi;
% 
%     patchR=0.5;
%     alphavals=expfitresultG(:,8);
%     maxAlpha = quantile(alphavals,0.80);
%     alphavals=alphavals/maxAlpha;
%     morex=alphavals>=1.0;
%     alphavals(morex)=1;
%     for i=1:size(expfitresultG(:,7))
%         %ddy1=expfitresultR(i,3)+expfitresultR(i,4);
%         %ddy2=expfitresultR(i,3)-expfitresultR(i,4);
%         ddy1=expfitresultR(i,4);
%        % if(bFitGaussianExp)
%             ddy2=expfitresultR(i,4)-1/expfitresultR(i,6);
%         %else
%           %  ddy2=expfitresultR(i,3)+expfitresultR(i,nFitParamNR);
%         %end
% 
%         line([expfitresultR(i,nFitParamNR+1) expfitresultR(i,nFitParamNR+1)],[ddy1 ddy2],'Color','r');
% 
%         %{
%         x=expfitresultG(i,nFitParamNG+1);
%         y=expfitresultG(i,4)+expfitresultG(i,6);
%         
%         %accumulated intensity
%         pb=patch((patchR*sin(t)+ x),(patchR*cos(t)+y),'g','edgecolor','none','FaceAlpha',alphavals(i));
%         %}
%     end
% 
%     %position of error function
%    % plot(expfitresultG(:,nFitParamNG+1),expfitresultG(:,4),'Color','g');
%     %plot(expfitresultR(:,nFitParamNR+1),expfitresultR(:,3),'Color','r');
%      
% end
% 
% %% MTA analysis
% 
%  xvals=expfitresultR(:,7);
%  yvals=smooth(expfitresultR(:,4),3);
%  xyArr=horzcat(xvals,yvals);
% [optimal_epoches, slopes, xyApprox] = mta_analysis(xyArr,2);
% 
%  xvals=expfitresultR(:,7);
%  yvals=expfitresultR(:,4);
% 
% %approximation
% %  nSlopesFound = length(slopes);
% %  for i=1:nSlopesFound
% %      rangeApprox = ((optimal_epoches(i)+i-1):(optimal_epoches(i+1)+i-1))';      
%       %plot(xyArrNoNoise(rangeApprox,1),xyArrNoNoise(rangeApprox,2),':r','LineWidth',2);
%        
%    %   plot(xyApprox(rangeApprox,1),xyApprox(rangeApprox,2),'Color',[0 0.5 0],'LineWidth',2);
%     %  hold on
% %  end
% 
%   
% %%fill average speed values
% MTASpeed=zeros(length(xvals),1);
% nEpN=length(slopes);
% for(ind=1:nEpN)
%     MTASpeed(optimal_epoches(ind):optimal_epoches(ind+1),1)=slopes(ind);
% end
% %average slope number in the transiotion (epoch) points
% if(nEpN>2)
%     for(ind=2:nEpN)
%         MTASpeed(optimal_epoches(ind),1)=mean([slopes(ind-1) slopes(ind)]);
%     end
% end
% 
% if(bShowPlot)
%     %{
%     maxIm=max(max(I_R));
%     minIm=min(min(I_R));
%     I2R=(I_R-minIm)/(maxIm-minIm);
%     handlefig=figure;
%     imshow(I2R); 
%     hold on
%     plot(xvals,yvals,'Color','r');
%     hold on
%     %}
%     plot(xyApprox(:,1),xyApprox(:,2),'Color','g');
%     hold off
%      saveas(handlefig,strcat(filebase,'_R.fig'),'fig');
% end
% 
% expfitresultR=horzcat(expfitresultR,MTASpeed);