/* Decompress
 Project
 SF2568 Parallel Computations for Large-Scale Problems
 Josephine Thuvander
 josthu@kth.se
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <jpeglib.h>


/* Function declarations */
void writeFile(int s, unsigned char *im, char file[]);
void readFile(int s, char file[], int *im);
void printImage(int h, int w, unsigned char *im);
int decompress(char file[], unsigned char *im);
int compress(int w, int h, char in_file[], char out_file[]);
void prewitt(int h, int w, unsigned char im_in[h][w], unsigned char out_im[h][w]);
void ptrToMap(int h, int w, unsigned char *im_ptr, unsigned char im_map[h][w]);
void mapToPtr(int h, int w, unsigned char im_map[h][w], unsigned char *im_ptr);
void sobel_x(int h, int w, unsigned char im_in[h][w], unsigned char im_out[h][w]);
void sobel_y(int h, int w, unsigned char im_in[h][w], unsigned char im_out[h][w]);
/* Functions definitions */ 

void prewitt(int h, int w, unsigned char im_in[h][w], unsigned char out_im[h][w]) {
	int i, j;
	/* For now, skip the outermost pixels */
	int t1 , t2;
	unsigned char min, max;
	min = 255;
	max = 0;

	/* Set edges to 0 */
	for (i=0;i<h; i++) {
		out_im[i][0] = 0;
		out_im[i][w-1] = 0;
	}

	for (j=1;j<w; j++) {
		out_im[0][j] = 0;
		out_im[h-1][j] = 0;
	}

	for (i=1; i<h-1; i++) {
		for (j=1; j<w-1; j++) {
			t1 = im_in[i+1][j-1]-im_in[i-1][j-1]+
				im_in[i+1][j]-im_in[i-1][j]+
				im_in[i+1][j+1]-im_in[i-1][j+1];
			t2 = im_in[i-1][j+1]-im_in[i-1][j-1]+
				im_in[i][j+1]-im_in[i][j-1]+
				im_in[i+1][j+1]-im_in[i+1][j-1];
			//out_im[i][j] = abs(t1) + abs(t2);
			out_im[i][j] = sqrt((t1*t1)+(t2*t2));

			if (out_im[i][j] < min) {
				min = out_im[i][j];
			}
			if (out_im[i][j] > max) {
				max = out_im[i][j];
			}
		}
	}

	/* Scale the values */
	for (i=1; i<h-1; i++) {
		for (j=1; j<w-1; j++) {
			out_im[i][j] = 255*(out_im[i][j]-min)/(max-min);
		}
	}
}


void sobel_x(int h, int w, unsigned char im_in[h][w], unsigned char im_out[h][w]) {
	int i,j;

	/* Set edges to 0 */
	
	for (i=0; i<h; i++) {
		im_out[i][0] = 0;
		im_out[i][w-1] = 0;
	}

	for (i=0; i<h; i++) {
		for (j=0; j<w; j++) {
			
		}
		
	}

}


void sobel_y(int h, int w, unsigned char im_in[h][w], unsigned char im_out[h][w]) {
	int i,j;

	for (j=0; j<w; j++) {

	}

	for (j=1;j<w; j++) {
		im_out[0][j] = 0;
		im_out[h-1][j] = 0;
	}


}



void writeFile(int s, unsigned char *im ,char file[]) {
	printf("Writing to file...\n");
	printf("Size is %d \n", s);
	FILE *fp;

	fp = fopen(file,"w");

   	if(fp == NULL) {
   		printf("Not able to open file \n");
		exit(EXIT_FAILURE);
   	}

	fwrite(im,sizeof(unsigned char),s,fp);
	fclose(fp);
}

void readFile(int s, char file[], int *im) {
	FILE *fp;
	int i;

	fp = fopen(file,"r");

	if(fp == NULL) {
		printf("No such file... \n");
		exit(EXIT_FAILURE);
   	}

   	fread(im,sizeof(int),s,fp);
   	fclose(fp);
}

void printImage(int h, int w, unsigned char *im) {
	int i, j;
	printf("Printing image...\n");

	for (i=0; i<h; i++) {
		for (j=0; j<w; j++) {
			printf("%u ",im[i*w+j]);
		}
		printf("\n");
	}

 	printf("\n");
}

int decompress(char in_file[], unsigned char *im) {
	unsigned char *jdata;
	struct jpeg_decompress_struct cinfo;
    struct jpeg_error_mgr jerr;

    unsigned long x, y, data_size, bmp_size;
    int i, pixel_size, row_stride,channels;

    FILE *file = fopen(in_file,"r");

    cinfo.err = jpeg_std_error(&jerr);
    jpeg_create_decompress(&cinfo);
    
    if (!file) {
    	fprintf(stderr,"Error reading %s \n", in_file);
     	return 0;
    } 

    jpeg_stdio_src(&cinfo, file);
    jpeg_read_header(&cinfo, TRUE); // read jpeg file header
    jpeg_start_decompress(&cinfo); // decompress the file
    
    cinfo.out_color_space = JCS_GRAYSCALE;
    x = cinfo.output_width;
    y = cinfo.output_height;
    pixel_size = cinfo.output_components;

    //printf("pixel_size %u \n", pixel_size);

    bmp_size = x * y * pixel_size;

    row_stride = x * pixel_size;
    
    data_size = x*y;

    //jdata = (unsigned char *)malloc(data_size);
    jdata = (unsigned char *) malloc(bmp_size);

    while (cinfo.output_scanline < cinfo.output_height) {
      unsigned char *buffer_array[1];
      buffer_array[0] = jdata + \
               (cinfo.output_scanline) * row_stride;
      jpeg_read_scanlines(&cinfo, buffer_array, 1);
    }

    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);

    // unsigned char c[data_size*pixel_size];

    for (i=0; i<data_size*pixel_size; i++) {
      //c[i] = jdata[i];
    	im[i] = jdata[i];
      	//fprintf(resultfiletxt,"%u ", c[i]);
      	//printf("%u \n",im[i]);
    }
    fclose(file);
    return 0;
}

int compress(int w, int h, char in_file[], char out_file[]) {
	
	int quality = 50;

  	//FILE *buffer = fopen("decompressed_file.bin","r");
  	//FILE *outfile = fopen("compressed_image.jpg", "wb");

	FILE *buffer = fopen(in_file,"r");
  	FILE *outfile = fopen(out_file, "wb");	

  	unsigned char *jdata;
  	// Jpeg compression object
  	struct jpeg_compress_struct cinfo;

  	// For error handling
	struct jpeg_error_mgr jerr;
 
	cinfo.err = jpeg_std_error(&jerr);
	
	jpeg_create_compress(&cinfo);

	// Data destination module
	jpeg_stdio_dest(&cinfo, outfile);

 	// Parameters of compression
	cinfo.image_width = w;
	cinfo.image_height = h;
	cinfo.input_components = 1;
	cinfo.in_color_space = JCS_GRAYSCALE;

	jpeg_set_defaults(&cinfo);
	/*set the quality [0..100]  */
	jpeg_set_quality (&cinfo, 75, 1);

	// Begin a compression cycle
	jpeg_start_compress(&cinfo, TRUE);

	JSAMPROW row_pointer[1];          /* pointer to a single row */
 	int i = 0; 

 	int row_stride;			/* physical row width in buffer */
 	row_stride = cinfo.image_width * cinfo.input_components;	/* JSAMPLEs per row in image_buffer */

 	printf("Row stride %d \n", row_stride);
 	// jdata = (unsigned char *)malloc(400*500);

	printf("Before while: Scan line %d \n",cinfo.next_scanline);
	printf("Before while: cinfo height  %d \n",cinfo.image_height);

	unsigned char *b = (unsigned char*) malloc(row_stride);

	while (cinfo.next_scanline < cinfo.image_height) {
		//row_pointer[0] = (JSAMPROW) &buffer[cinfo.next_scanline*row_stride];
		fseek(buffer,cinfo.next_scanline*row_stride,SEEK_SET);
		fread(b,sizeof(unsigned char),row_stride,buffer);
		//row_pointer[0] = (JSAMPROW) &b;
		row_pointer[0] = b;

		jpeg_write_scanlines(&cinfo, row_pointer, 1);
		   	    
	    printf("Scan line %d \n",cinfo.next_scanline);
	    printf("next*row_stride %d \n", cinfo.next_scanline*row_stride);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);
	
	fclose(buffer);
	fclose(outfile);
  	return 0;
}

void ptrToMap(int h, int w, unsigned char *im_ptr, unsigned char im_map[h][w]) {
	int i, j;
	for (i=0; i<h; i++) {
		for (j=0; j<w; j++) {
			im_map[i][j] = im_ptr[i*w+j];
		}
	}
	printf("Image is rearranged to 2d map \n");
}

void mapToPtr(int h, int w, unsigned char im_map[h][w], unsigned char *im_ptr) {
	int i, j;
	for (i=0; i<h; i++) {
		for (j=0; j<w; j++) {
			im_ptr[i*w+j] = im_map[i][j];
		}
	}
	printf("Image is rearranged to ptr \n");
}


/* Main program */
int main(int argc, char **argv) {

	int i, width, height, pixel_size, size;
	unsigned char *decompressed_image, *compressed_image;
	
	char filename[] = "grayimage.jpg";
	char decompressed_file[] = "decompressed_file.bin";
	char compressed_file[] = "compressed_file.jpg";

	//"decompressed_file.bin","r"

	width = 512;
	height = 512;
	pixel_size = 1;
	size = width*height*pixel_size;

	// Allocate memory
	decompressed_image = (unsigned char *) malloc(size*sizeof(unsigned char));
	compressed_image = (unsigned char *) malloc(size*sizeof(unsigned char));

	int done;
	unsigned char image_map[width][height];
	unsigned char detected_map[width][height];

	// Decompress file
	done = decompress(filename,decompressed_image);
	//printImage(width,height,decompressed_image);
	ptrToMap(height,width,decompressed_image,image_map);
	prewitt(height,width,image_map,detected_map);
	mapToPtr(height,width,detected_map,decompressed_image);
	writeFile(size,decompressed_image,decompressed_file);

	//done  = compress(compressed_file, decompressed_image);
	done = compress(width, height, decompressed_file, compressed_file);

	/*
	printImage(size, input_image);
	
	writeFile(size,input_file,input_image);
	readFile(size, input_file, output_image);
	printImage(size, output_image);

	*/

	return 0;
}