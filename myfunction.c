#include <stdbool.h>
#include "stdio.h"
#include "time.h"
typedef struct {
   unsigned char red;
   unsigned char green;
   unsigned char blue;
} pixel;

typedef struct {
    int red;
    int green;
    int blue;
    // int num;
} pixel_sum;


/* Compute min and max of two integers, respectively */
int min(int a, int b) { return (a < b ? a : b); }
int max(int a, int b) { return (a > b ? a : b); }

int calcIndex(int i, int j, int n) {
	return ((i*n)+(j));
}

/*
 * initialize_pixel_sum - Initializes all fields of sum to 0
 */
void initialize_pixel_sum(pixel_sum *sum) {
	sum->red = sum->green = sum->blue = 0;
	// sum->num = 0;
	return;
}

/*
 * assign_sum_to_pixel - Truncates pixel's new value to match the range [0,255]
 */
static void assign_sum_to_pixel(pixel *current_pixel, pixel_sum sum, int kernelScale) {

	// divide by kernel's weight
	sum.red = sum.red / kernelScale;
	sum.green = sum.green / kernelScale;
	sum.blue = sum.blue / kernelScale;

	// truncate each pixel's color values to match the range [0,255]
	current_pixel->red = (unsigned char) (min(max(sum.red, 0), 255));
	current_pixel->green = (unsigned char) (min(max(sum.green, 0), 255));
	current_pixel->blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return;
}

/*
* sum_pixels_by_weight - Sums pixel values, scaled by given weight
*/
static void sum_pixels_by_weight(pixel_sum *sum, pixel p, int weight) {
	sum->red += (p.red) * weight;
	sum->green += (p.green) * weight;
	sum->blue += ( p.blue) * weight;
	// sum->num++;
	return;
}

/*
 *  Applies kernel for pixel at (i,j)
 */
static pixel applyKernel(int dim, int i, int j, pixel *src, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int ii, jj;
	pixel_sum sum;
	pixel current_pixel;
	int min_intensity = 766; // arbitrary value that is higher than maximum possible intensity, which is 255*3=765
	int max_intensity = -1; // arbitrary value that is lower than minimum possible intensity, which is 0
	int min_row, min_col, max_row, max_col;
    int k=0;
    int kRow, kCol;
	pixel loop_pixel;
	sum.red = sum.green = sum.blue = 0;
	// sum->num = 0;
    int minidim=min(i+1,dim-1);
    int mindjm=min(j+1,dim-1);
    int weight;
    int x=0;
    int x1=0;
    //if i make optimize here it will reduce runtime in 10  ms
	for(ii = i-1; ii <= minidim; ii++,x++) {
        x1=0;
		for(jj = j-1; jj <= mindjm; jj++,x1++) {
			// compute row index in kernel
            k=ii*dim+jj;
			// apply kernel on pixel at [ii,jj]
             weight=kernel[x][x1];
            sum.red += src[k].red * weight ;
            sum.green += src[k].green * weight;
            sum.blue += src[k].blue * weight;
		}
	}
	int colorsum=0;

	if (filter) {
		// find min and max coordinates
        ii = i-1;
		for(; ii <= minidim; ii++) {
            jj = j-1;
			for(; jj <= mindjm; jj++) {
				// check if smaller than min or higher than max and update
				loop_pixel = src[(ii*dim)+jj];
				colorsum= loop_pixel.red +  loop_pixel.green + loop_pixel.blue;
				if (colorsum <= min_intensity) {
					min_intensity = colorsum;
					min_row = ii;
					min_col = jj;
				}
				if (colorsum > max_intensity) {
					max_intensity = colorsum;
					max_row = ii;
					max_col = jj;
				}
			}
		}
		// filter out min and max
        k=min_row*dim+min_col;
        sum.red -= src[k].red  ;
        sum.green -= src[k].green ;
        sum.blue -= src[k].blue ;
        k=max_row*dim+max_col;
        sum.red -= src[k].red  ;
        sum.green -= src[k].green ;
        sum.blue -= src[k].blue ;
	}

	// assign kernel's result to pixel at [i,j]
	//assign_sum_to_pixel(&current_pixel, sum, kernelScale);
    sum.red = sum.red / kernelScale;
    sum.green = sum.green / kernelScale;
    sum.blue = sum.blue / kernelScale;
    // truncate each pixel's color values to match the range [0,255]
    current_pixel.red = (unsigned char) (min(max(sum.red, 0), 255));
    current_pixel.green = (unsigned char) (min(max(sum.green, 0), 255));
    current_pixel.blue = (unsigned char) (min(max(sum.blue, 0), 255));
	return current_pixel;
}

/*
* Apply the kernel over each pixel.
* Ignore pixels where the kernel exceeds bounds. These are pixels with row index smaller than kernelSize/2 and/or
* column index smaller than kernelSize/2
*/
void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	int i, j;
	for (i = kernelSize / 2 ; i < dim - kernelSize / 2; i++) {
		for (j =  kernelSize / 2 ; j < dim - kernelSize / 2 ; j++) {
			dst[i*dim +j] = applyKernel(dim, i, j, src, kernelSize, kernel, kernelScale, filter);
		}
	}
}

void charsToPixels(Image *charsImg, pixel* pixels) {

	
	int i=0;
	while(i<n*m){
	               pixels[i].red = image->data[3*i];
			pixels[i].green = image->data[3*i+ 1];
			pixels[i].blue = image->data[3*i + 2];
			++i;
			if(i<n*m){
			pixels[i].red = image->data[3*i];
			pixels[i].green = image->data[3*i+ 1];
			pixels[i].blue = image->data[3*i + 2];
			++i;
			}
			if(i<n*m){
			pixels[i].red = image->data[3*i];
			pixels[i].green = image->data[3*i+ 1];
			pixels[i].blue = image->data[3*i + 2];
			++i;
			}
			if(i<n*m){
			pixels[i].red = image->data[3*i];
			pixels[i].green = image->data[3*i+ 1];
			pixels[i].blue = image->data[3*i + 2];
			++i;
			}
			
	}
}

void pixelsToChars(pixel* pixels, Image *charsImg) {
int row, col;
	int ni = 0;
	for (row =0 ; row < m ; row++) {

		for ( col=0; col < n ; col++) {
		int z = 3*ni + 3*col;
			image->data[z] = pixels[ni + col].red;
			image->data[z + 1] = pixels[ni + col].green;
			image->data[z + 2] = pixels[ni + col].blue;
		}
		ni+=m;
	}
 
}
void copyPixels(pixel* src, pixel* dst) {
	int i=0;
	while(i<n*m){
	
	               dst[i].red = src[i].red;
			dst[i].green = src[i].green;
			dst[i].blue = src[i].blue;
			++i;
			if(i<n*m){
			 dst[i].red = src[i].red;
			dst[i].green = src[i].green;
			dst[i].blue = src[i].blue;
			++i;
			}
			if(i<n*m){
			 dst[i].red = src[i].red;
			dst[i].green = src[i].green;
			dst[i].blue = src[i].blue;
			++i;
			}
			if(i<n*m){
			 dst[i].red = src[i].red;
			dst[i].green = src[i].green;
			dst[i].blue = src[i].blue;
			++i;
			}
	}



}

void doConvolution(Image *image, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter) {

	pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	pixel* backupOrg = malloc(m*n*sizeof(pixel));

	charsToPixels(image, pixelsImg);
	copyPixels(pixelsImg, backupOrg);

	smooth(m, backupOrg, pixelsImg, kernelSize, kernel, kernelScale, filter);

	pixelsToChars(pixelsImg, image);

	free(pixelsImg);
	free(backupOrg);
}

void myfunction(Image *image, char* srcImgpName, char* blurRsltImgName, char* sharpRsltImgName, char* filteredBlurRsltImgName, char* filteredSharpRsltImgName, char flag) {

	/*
	* [1, 1, 1]
	* [1, 1, 1]
	* [1, 1, 1]
	*/
	int blurKernel[3][3] = {{1, 1, 1}, {1, 1, 1}, {1, 1, 1}};

	/*
	* [-1, -1, -1]
	* [-1, 9, -1]
	* [-1, -1, -1]
	*/
    clock_t start, end;
    start=clock();
	int sharpKernel[3][3] = {{-1,-1,-1},{-1,9,-1},{-1,-1,-1}};
    register int i=0;
    register int j=0;
   register int nm=n*m;
    int k=nm/10;
   register int c=0;
    int k2=nm%10;
    int im=0;
    int msubone=m-1;
	if (flag == '1') {	
		// blur image  doconvlution
		//doConvolution(image, 3, blurKernel, 9, false);
		pixel* pixelsImg = malloc(m*n*sizeof(pixel));
	        pixel* backupOrg = malloc(m*n*sizeof(pixel));
		//charsToPixels(image, pixelsImg);
         i=0;
         j=0;
         c=k;
        while(c>0){
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }//good
        c=k2;
        while(c>0){
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }
    //copy pixels
	//copyPixels(pixelsImg, backupOrg);
	 i=0;
       c=k;
	while(c>0){
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        --c;
	}//good
    c=k2;
    while(c>0){
        backupOrg[i].red = pixelsImg[i].red;
        backupOrg[i].green = pixelsImg[i].green;
        backupOrg[i].blue = pixelsImg[i].blue;
        ++i;
        --c;
    }
    //smooth
	//void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);
         im=m;
    start=clock();
        int msubone=m-1;
        for (i = 1 ; i < msubone; i++) {
            for (j =  1 ; j < msubone ; j++) {
                ++im;
                pixelsImg[im] = applyKernel(m, i, j, backupOrg, 3, blurKernel,9, false);
            }
            im-=j-1;
            im+=m;
        }//good
        end=clock()-start;
        printf("%f",(double )end/CLOCKS_PER_SEC);
    //pixels to chars
        i=0;
        j=0;
        c=k;
        while(c>0){

                image->data[i] = pixelsImg[j].red;
                ++i;
                image->data[i] = pixelsImg[j].green;
                ++i;
                image->data[i] = pixelsImg[j].blue;
                ++i;
                ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
            }
        c=k2;
        while(c>0){
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
        }


        free(pixelsImg);
        free(backupOrg);
//all good untill here

		// write result image to file
		writeBMP(image, srcImgpName, blurRsltImgName);	

		// sharpen the resulting image
		//doConvolution(image, 3, sharpKernel, 1, false);

         pixelsImg = malloc(m*n*sizeof(pixel));
         backupOrg = malloc(m*n*sizeof(pixel));
        //charsToPixels(image, pixelsImg);
         i=0;
         j=0;
         c=k;
        while(c>0){
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }//good
        c=k2;
        while(c>0){
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }
        //copy pixels
        //copyPixels(pixelsImg, backupOrg);
        i=0;
        c=k;
        while(c>0){
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            --c;
        }//good
        c=k2;
        while(c>0){
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            --c;
        }
        //smooth
        //void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);
         im=m;

        start=clock();
        for (i = 1 ; i < msubone; i++) {
            for (j =  1 ; j < msubone ; j++) {
                ++im;
                pixelsImg[im] = applyKernel(m, i, j, backupOrg, 3, sharpKernel,1, false);
            }
            im-=j-1;
            im+=m;
        }
        end=clock()-start;
        printf("%f",(double )end/CLOCKS_PER_SEC);
        //pixels to chars

        i=0;
        j=0;
        c=k;
        while(c>0){

            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
        }
        c=k2;
        while(c>0){
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
        }
        free(pixelsImg);
        free(backupOrg);
		
		// write result image to file
		writeBMP(image, srcImgpName, sharpRsltImgName);

	} else {
		// apply extermum filtered kernel to blur image
		//doConvolution(image, 3, blurKernel, 7, true);
        // blur image  doconvlution

        pixel* pixelsImg = malloc(m*n*sizeof(pixel));
        pixel* backupOrg = malloc(m*n*sizeof(pixel));
        //charsToPixels(image, pixelsImg);
         i=0;
         j=0;

        c=k;
        while(c>0){
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }//good
        c=k2;
        while(c>0){
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }
        //copy pixels
      //  copyPixels(pixelsImg, backupOrg);
        i=0;
        c=k;
        while(c>0){
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            --c;
        }//good
        c=k2;
        while(c>0){
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            --c;
        }
        //smooth
        //void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);

         im=m;
        for (i = 1 ; i < msubone; i++) {
            for (j =  1 ; j < msubone ; j++) {
                ++im;
                pixelsImg[im] =  applyKernel(m, i, j, backupOrg, 3, blurKernel,7, true);
            }
            im-=j-1;
            im+=m;
        }//good
        //pixels to chars

        i=0;
        j=0;
        c=k;
        while(c>0){

            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
        }
        c=k2;
        while(c>0){
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
        }
        free(pixelsImg);
        free(backupOrg);
//all good untill here


        // write result image to file
		writeBMP(image, srcImgpName, filteredBlurRsltImgName);

		// sharpen the resulting image
		//doConvolution(image, 3, sharpKernel, 1, false);
        //doConvolution(image, 3, sharpKernel, 1, false);

        pixelsImg = malloc(m*n*sizeof(pixel));
        backupOrg = malloc(m*n*sizeof(pixel));
        //charsToPixels(image, pixelsImg);
         i=0;
         j=0;
         c=k;
        while(c>0){
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            j=3*i;
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }//good
        c=k2;
        while(c>0){
            pixelsImg[i].red = image->data[j];
            ++j;
            pixelsImg[i].green = image->data[j];
            ++j;
            pixelsImg[i].blue = image->data[j];
            ++j;
            ++i;
            --c;
        }
        //copy pixels
        //copyPixels(pixelsImg, backupOrg);
        i=0;
        c=k;
        while(c>0){
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            --c;
        }//good
        c=k2;
        while(c>0){
            backupOrg[i].red = pixelsImg[i].red;
            backupOrg[i].green = pixelsImg[i].green;
            backupOrg[i].blue = pixelsImg[i].blue;
            ++i;
            --c;
        }
        //smooth
        //void smooth(int dim, pixel *src, pixel *dst, int kernelSize, int kernel[kernelSize][kernelSize], int kernelScale, bool filter);

        im=m;

        for (i = 1 ; i < msubone; i++) {
            for (j =  1 ; j < msubone ; j++) {
                ++im;
                pixelsImg[im] = applyKernel(m, i, j, backupOrg, 3, sharpKernel,1, false);
            }
            im-=j-1;
            im+=m;
        }//good
        //pixels to chars
        i=0;
        j=0;
        c=k;
        while(c>0){

            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
        }
        c=k2;
        while(c>0){
            image->data[i] = pixelsImg[j].red;
            ++i;
            image->data[i] = pixelsImg[j].green;
            ++i;
            image->data[i] = pixelsImg[j].blue;
            ++i;
            ++j;
            --c;
        }

        free(pixelsImg);
        free(backupOrg);
        // write result image to file
        end=clock()-start;
        printf("lolo%f\n",(double )end/CLOCKS_PER_SEC);
		writeBMP(image, srcImgpName, filteredSharpRsltImgName);

	}
}

