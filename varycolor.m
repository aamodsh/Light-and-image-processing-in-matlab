function ColorSet=varycolor(NumberOfPlots, opts)
% VARYCOLOR Produces colors with maximum variation on plots with multiple
% lines.
%
%     VARYCOLOR(X) returns a matrix of dimension X by 3.  The matrix may be
%     used in conjunction with the plot command option 'color' to vary the
%     color of lines.  
%
%     Yellow and White colors were not used because of their poor
%     translation to presentations.
% 
%     Example Usage:
%         NumberOfPlots=50;
%
%         ColorSet=varycolor(NumberOfPlots);
% 
%         figure
%         hold on;
% 
%         for m=1:NumberOfPlots
%             plot(ones(20,1)*m,'Color',ColorSet(m,:))
%         end

%Created by Daniel Helmick 8/12/2008

%error(nargchk(1,1,nargin))%correct number of input arguements??
error(nargoutchk(0, 1, nargout))%correct number of output arguements??

if exist('opts', 'var')
    switch opts
        case 'r2g'
            N = NumberOfPlots / 2 ;
            greenColorMap = [zeros(1, floor(NumberOfPlots/2)), linspace(0, 1, ceil(NumberOfPlots/2))];
            redColorMap = [linspace(1, 0, ceil(NumberOfPlots/2)), zeros(1, floor(NumberOfPlots/2))];
            ColorSet = [redColorMap; greenColorMap; zeros(1, length(greenColorMap))]';
            return ;
    end
end

%Take care of the anomolies
if NumberOfPlots<1
    ColorSet=[];
elseif NumberOfPlots==1
    ColorSet=[0 0 0];
elseif NumberOfPlots==2
    ColorSet=[1 0 0; 0 0 1];
elseif NumberOfPlots==3
    ColorSet=[31, 118, 180; 254, 127, 14; 44, 160, 44] / 255;
elseif NumberOfPlots==4
    ColorSet=[31, 118, 180; 254, 127, 14; 44, 160, 44; 213, 39, 40] / 255;
elseif NumberOfPlots==5
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0];
elseif NumberOfPlots==6
    ColorSet=[0 1 0; 0 1 1; 0 0 1; 1 0 1; 1 0 0; 0 0 0];

else %default and where this function has an actual advantage

    %we have 5 segments to distribute the plots
    EachSec=floor(NumberOfPlots/5); 
    
    %how many extra lines are there? 
    ExtraPlots=mod(NumberOfPlots,5); 
    
    %initialize our vector
    ColorSet=zeros(NumberOfPlots,3);
    
    %This is to deal with the extra plots that don't fit nicely into the
    %segments
    Adjust=zeros(1,5);
    for m=1:ExtraPlots
        Adjust(m)=1;
    end
    
    SecOne   =EachSec+Adjust(1);
    SecTwo   =EachSec+Adjust(2);
    SecThree =EachSec+Adjust(3);
    SecFour  =EachSec+Adjust(4);
    SecFive  =EachSec;

    for m=1:SecOne
        ColorSet(m,:)=[0 1 (m-1)/(SecOne-1)];
    end

    for m=1:SecTwo
        ColorSet(m+SecOne,:)=[0 (SecTwo-m)/(SecTwo) 1];
    end
    
    for m=1:SecThree
        ColorSet(m+SecOne+SecTwo,:)=[(m)/(SecThree) 0 1];
    end
    
    for m=1:SecFour
        ColorSet(m+SecOne+SecTwo+SecThree,:)=[1 0 (SecFour-m)/(SecFour)];
    end

    for m=1:SecFive
        ColorSet(m+SecOne+SecTwo+SecThree+SecFour,:)=[(SecFive-m)/(SecFive) 0 0];
    end
    
end