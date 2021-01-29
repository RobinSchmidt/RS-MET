clear all;

height = 64;   % height of the rectangle in pixels
width  = 128;  % width of the rectangle in pixels

filename = strcat('bilinearGradienBlue', num2str(width), 'x', ...
                  num2str(height), '.bmp');

topLeftRed       = 1.0;  % red-value of the top-left corner
topRightRed      = 0.25;  % red-value of the top-left corner
bottomLeftRed    = 0.25;
bottomRightRed   = 1.0;
topLeftGreen     = 1.0;  
topRightGreen    = 0.25; 
bottomLeftGreen  = 0.25;
bottomRightGreen = 1.0;
topLeftBlue      = 1.0;  
topRightBlue     = 0.5; 
bottomLeftBlue   = 0.5;
bottomRightBlue  = 1.0;

% allocate gray-matrix:
colorMatrix = zeros(height, width, 3);

colorMatrix(1, 1, 1)          = topLeftRed;
colorMatrix(1, width, 1)      = topRightRed;
colorMatrix(height, 1, 1)     = bottomLeftRed;
colorMatrix(height, width, 1) = bottomRightRed;
colorMatrix(1, 1, 2)          = topLeftGreen;
colorMatrix(1, width, 2)      = topRightGreen;
colorMatrix(height, 1, 2)     = bottomLeftGreen;
colorMatrix(height, width, 2) = bottomRightGreen;
colorMatrix(1, 1, 3)          = topLeftBlue;
colorMatrix(1, width, 3)      = topRightBlue;
colorMatrix(height, 1, 3)     = bottomLeftBlue;
colorMatrix(height, width, 3) = bottomRightBlue;

for x=1:width
 if( width == 1 )
  p = 1;
 else
  p = (x-1)/(width-1); % relative x-postion (from 0...1)
 end
 colorMatrix(1, x, 1)      = (1-p)*topLeftRed      + p*topRightRed; 
 colorMatrix(1, x, 2)      = (1-p)*topLeftGreen    + p*topRightGreen; 
 colorMatrix(1, x, 3)      = (1-p)*topLeftBlue     + p*topRightBlue; 
 colorMatrix(height, x, 1) = (1-p)*bottomLeftRed   + p*bottomRightRed; 
 colorMatrix(height, x, 2) = (1-p)*bottomLeftGreen + p*bottomRightGreen; 
 colorMatrix(height, x, 3) = (1-p)*bottomLeftBlue  + p*bottomRightBlue; 
end

for x=1:width
 for y=2:height-1
  p = (y-1)/(height-1); % relative y-postion (from 0...1)  
  colorMatrix(y,x,1) = (1-p)*colorMatrix(1,x,1) + p*colorMatrix(height,x,1);  
  colorMatrix(y,x,2) = (1-p)*colorMatrix(1,x,2) + p*colorMatrix(height,x,2);    
  colorMatrix(y,x,3) = (1-p)*colorMatrix(1,x,3) + p*colorMatrix(height,x,3);    
 end 
end

image(colorMatrix);
imwrite(colorMatrix, filename, 'bmp');
%colormap(bone);


