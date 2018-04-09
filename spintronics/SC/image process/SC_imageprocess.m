clear;

picture=imread('1.bmp');
[m,n]=size(picture);
newpicture=picture;
robertsNum=0;

for i=1:1:m-1
    for j=1:1:n-1
        robertsNum = abs(picture(i,j)-picture(i+1,j+1)) + abs(picture(i+1,j)-picture(i,j+1));
        newpicture(i,j)=255-robertsNum;
    end
end

imshow(newpicture);
