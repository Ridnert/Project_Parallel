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
#include "mpi.h"
#include <time.h>


/* Function declarations */
void writeFile(int s, unsigned char *im, char file[]);
void readFile(int s, char file[], int *im);
void printImage(int h, int w, unsigned char *im);
int decompress(char file[], unsigned char *im);
int compress(int w, int h, char in_file[], char out_file[]);
void prewitt(int h, int w, unsigned char **im_in, unsigned char **out_im);
void ptrToMap(int h, int w, unsigned char *im_ptr, unsigned char **im_map);
void mapToPtr(int h, int w, unsigned char **im_map, unsigned char *im_ptr);
void sobel(int h, int w, unsigned char** im_in, unsigned char** im_out);
void convertToGray(int s, unsigned char* im_in, unsigned char* im_out);
unsigned char** make2DArray(int w, int h);

/* Functions definitions */
void writeFile(int s, unsigned char *im ,char file[]) {
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
    jpeg_read_header(&cinfo, TRUE);
    jpeg_start_decompress(&cinfo);
    cinfo.out_color_space = JCS_RGB;
    x = cinfo.output_width;
    y = cinfo.output_height;
    pixel_size = cinfo.output_components;
    bmp_size = x * y * pixel_size;
    row_stride = x * pixel_size;
    data_size = x*y;
    jdata = (unsigned char *) malloc(bmp_size);
    while (cinfo.output_scanline < cinfo.output_height) {
			//printf("Scanline %d \n", cinfo.output_scanline);
      unsigned char *buffer_array[1];
      buffer_array[0] = jdata + \
               (cinfo.output_scanline) * row_stride;
      jpeg_read_scanlines(&cinfo, buffer_array, 1);
    }
    jpeg_finish_decompress(&cinfo);
    jpeg_destroy_decompress(&cinfo);
    for (i=0; i<data_size*pixel_size-2; i++) {
    	im[i] = jdata[i];
    }
    fclose(file);
    return 0;
}

int compress(int w, int h, char in_file[], char out_file[]) {
		int quality = 75;
		FILE *buffer = fopen(in_file,"r");
	  FILE *outfile = fopen(out_file, "wb");
	  unsigned char *jdata;
	  struct jpeg_compress_struct cinfo;
		struct jpeg_error_mgr jerr;
		cinfo.err = jpeg_std_error(&jerr);
		jpeg_create_compress(&cinfo);
		jpeg_stdio_dest(&cinfo, outfile);

 	// Parameters of compression
	cinfo.image_width = w;
	cinfo.image_height = h;
	cinfo.input_components = 1;
	cinfo.in_color_space = JCS_GRAYSCALE;

	jpeg_set_defaults(&cinfo);
	jpeg_set_quality (&cinfo, 75, 1);
	jpeg_start_compress(&cinfo, TRUE);

	JSAMPROW row_pointer[1];	 /* row pointer */
 	int i = 0;
 	int row_stride;
 	row_stride = cinfo.image_width * cinfo.input_components;
	unsigned char *b = (unsigned char*) malloc(row_stride);

	while (cinfo.next_scanline < cinfo.image_height) {
		fseek(buffer,cinfo.next_scanline*row_stride,SEEK_SET);
		fread(b,sizeof(unsigned char),row_stride,buffer);
		row_pointer[0] = b;
		jpeg_write_scanlines(&cinfo, row_pointer, 1);
	}

	jpeg_finish_compress(&cinfo);
	jpeg_destroy_compress(&cinfo);

	fclose(buffer);
	fclose(outfile);
  return 0;
}

void ptrToMap(int h, int w, unsigned char *im_ptr, unsigned char **im_map) {
	int i, j;
	for (i=0; i<h; i++) {
		for (j=0; j<w; j++) {
			im_map[i][j] = im_ptr[i*w+j];
		}
	}
}

void mapToPtr(int h, int w, unsigned char **im_map, unsigned char *im_ptr) {
	int i, j;
	for (i=0; i<h; i++) {
		for (j=0; j<w; j++) {
			im_ptr[i*w+j] = im_map[i][j];
		}
	}
}

void sobel(int h, int w, unsigned char** im_in, unsigned char** output_image) {
	int i, j, dx, dy;
	unsigned char max, min;
	max = 0;
	min = 255;

	for (i=0;i<h-2; i++) {
		output_image[i][0] = 0;
		output_image[i][w-1] = 0;
	}

	for (i=1; i<h-1; i++) {
		for (j=1; j<w-1; j++) {
			dy = im_in[i+1][j-1]-im_in[i-1][j-1]+
				2*im_in[i+1][j]-2*im_in[i-1][j]+
				im_in[i+1][j+1]-im_in[i-1][j+1];
			dx = im_in[i-1][j+1]-im_in[i-1][j-1]+
				2*im_in[i][j+1]-2*im_in[i][j-1]+
				im_in[i+1][j+1]-im_in[i+1][j-1];
			output_image[i-1][j] = sqrt((dx*dx)+(dy*dy));
		}
	}

  	for (i = 0; i < h-2; i++){
    	for(j = 0; j < w; j++){
      		if(output_image[i][j] > max){
        		max = output_image[i][j];
      		}
					if(output_image[i][j] < min){
        		min = output_image[i][j];
      		}

    	}
  	}
  	/* Normalizing the array to values between 0-255 */
  	for (i = 0; i < h-2; i++){
    	for(j = 0; j < w; j++){
        	output_image[i][j] = 255*(output_image[i][j]-min)/(max-min);
      	}
    }
}

void prewitt(int h, int w, unsigned char **im_in, unsigned char **out_im) {
	int i, j, dx, dy;
	/* For now, skip the outermost pixels */
	unsigned char min, max;
	min = 255;
	max = 0;

	/* Set edges to 0 */
	for (i=0;i<h-2; i++) {
		out_im[i][0] = 0;
		out_im[i][w-1] = 0;
	}

	for (i=1; i<h-1; i++) {
		for (j=1; j<w-1; j++) {
			dy = im_in[i+1][j-1]-im_in[i-1][j-1]+
				im_in[i+1][j]-im_in[i-1][j]+
				im_in[i+1][j+1]-im_in[i-1][j+1];
			dx = im_in[i-1][j+1]-im_in[i-1][j-1]+
				im_in[i][j+1]-im_in[i][j-1]+
				im_in[i+1][j+1]-im_in[i+1][j-1];
			out_im[i-1][j] = sqrt((dy*dy)+(dx*dx));
			if (out_im[i-1][j] < min) {
				min = out_im[i-1][j-1];
			}
			if (out_im[i-1][j] > max) {
				max = out_im[i-1][j];
			}
		}
	}

	/* Scale the values */
	for (i=0; i<h-2; i++) {
		for (j=0; j<w; j++) {
			out_im[i][j] = 255*(out_im[i][j]-min)/(max-min);
		}
	}
}

void convertToGray(int s,unsigned char * im_in, unsigned char *im_out) {
	int i;
	for (i=0; i<s/3; i++) {
		im_out[i] = 0.2126*im_in[3*i] + 0.7152*im_in[(3*i)+1] + 0.0722*im_in[(3*i)+2];
	}
}

unsigned char** make2DArray(int w, int h) {
		int i;
		unsigned char **theArray;
		theArray = (unsigned char**) malloc(w*sizeof(unsigned char*));
		for (i=0; i<w; i++) {
				theArray[i] = (unsigned char*) malloc(h*sizeof(unsigned char));
		}
		return theArray;
}

/* Main program */
int main(int argc, char **argv) {

	int i, j, width, height, pixel_size, size;
	unsigned char *decompressed_image, *compressed_image, *decompressed_image_local;
	unsigned char *converted_image;

	char filename[] = "world.jpg";
	char decompressed_file[] = "decompressed_file.bin";
	char compressed_file[] = "compressed_file.jpg";

	/* For parallelization */
	int rc, P, rank, tag, subset_height, pixels_left, I,height_offset;

	width = 21600;
	height = 21600;
	pixel_size = 3;
	size = width*height*pixel_size;

	MPI_Status status;
  tag = 100;

  /* Initialize MPI */
  rc = MPI_Init(&argc, &argv);
  rc = MPI_Comm_size(MPI_COMM_WORLD, &P);
  rc = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	double start,end;
	//MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */

    if(rank == 0) {
			start = MPI_Wtime();
		}
  if (width < P) {
	   fprintf(stdout, "The image is to narrow...\n");
	   exit(1);
  }

	subset_height = height/P;		/* load balance */
	pixels_left = height%P; 		/* remaining pixels */
	height_offset = rank*subset_height;

	/* Allocate memory */
	decompressed_image = (unsigned char *) malloc(size*sizeof(unsigned char));
	compressed_image = (unsigned char *) malloc(width*height*sizeof(unsigned char));
	converted_image = (unsigned char *) malloc(width*height*sizeof(unsigned char));

	int done;
	unsigned char** image_map = make2DArray(height,width);

	// Decompress file
	done = decompress(filename,decompressed_image);

	convertToGray(size,decompressed_image,converted_image);
	//printImage(width,height,decompressed_image);
	ptrToMap(height-pixels_left,width,converted_image,image_map);
	height=height-pixels_left;
	/* Three cases*/
	/* If first process */
	if(rank == 0){
		unsigned char **image_map_local;
		unsigned char **detected_map;
		if (P==1) {
			image_map_local = make2DArray(subset_height,width);
			detected_map = make2DArray(subset_height-2,width);
			decompressed_image_local = (unsigned char *) malloc(width*(subset_height)*sizeof(unsigned char));
			for (i=0; i < subset_height; i++){
				for (j=0; j<width; j++){
					image_map_local[i][j] = image_map[i][j];
				}
			}
			//prewitt(subset_height,width,image_map_local,detected_map);
			sobel(subset_height,width,image_map_local,detected_map);
			mapToPtr(subset_height-2,width,detected_map,decompressed_image);
		} else {
			image_map_local = make2DArray(subset_height+1,width);
			detected_map = make2DArray(subset_height-1,width);
			decompressed_image_local = (unsigned char*) malloc(width*(subset_height+1)*sizeof(unsigned char));
			for (i=0; i < subset_height+1; i++){
				for (j=0; j<width; j++){
					image_map_local[i][j] = image_map[i][j];
				}
			}
			//prewitt(subset_height+1,width,image_map_local,detected_map);
			sobel(subset_height+1,width,image_map_local,detected_map);
			mapToPtr(subset_height-1,width,detected_map,decompressed_image);
		}
	} else if (rank == P-1){ 	/* If last process */
		unsigned char **image_map_local = make2DArray(subset_height+1,width);
		unsigned char **detected_map = make2DArray(subset_height-1,width);
	 	decompressed_image_local = (unsigned char*) malloc(width*(subset_height+1)*sizeof(unsigned char));
	 	for (i=height_offset-1; i< height; i++) { // Taking one extra column
		 	for(j=0; j < width ; j++) {
			 	image_map_local[i-height_offset+1][j] = image_map[i][j];
			}
		}
		//prewitt(subset_height+1,width,image_map_local,detected_map);
		sobel(subset_height+1,width,image_map_local,detected_map);
		mapToPtr(subset_height-1,width,detected_map,decompressed_image_local);
 	} else {
		unsigned char **image_map_local = make2DArray(subset_height+2,width);
		unsigned char **detected_map = make2DArray(subset_height,width);

 		decompressed_image_local = (unsigned char*) malloc(width*(subset_height+2)*sizeof(unsigned char));

		for (i = height_offset-1 ; i < height_offset+subset_height+1  ; i++){
			for(j = 0; j < width; j++){
				image_map_local[i-height_offset+1][j] = image_map[i][j];
			}
		}
		//prewitt(subset_height+2,width,image_map_local,detected_map);
		sobel(subset_height+2,width,image_map_local,detected_map);
		mapToPtr(subset_height,width,detected_map,decompressed_image_local);
 	}

	if (P>1) {
		int offset = (subset_height-1)*width;
		if(rank != 0 ){
			if(rank == P-1) {
				MPI_Send(&decompressed_image_local[0],(subset_height-1)*width,MPI_UNSIGNED_CHAR,0,tag,MPI_COMM_WORLD);
			} else {
				MPI_Send(&decompressed_image_local[0],(subset_height)*width,MPI_UNSIGNED_CHAR,0,tag,MPI_COMM_WORLD);
			}
		} else{
			for (i = 1; i<P-1; i++){
					MPI_Recv(&decompressed_image[offset],(subset_height)*width,MPI_UNSIGNED_CHAR,i,tag,MPI_COMM_WORLD,&status);
					offset += subset_height*width;
			}
			MPI_Recv(&decompressed_image[offset],(subset_height-1)*width,MPI_UNSIGNED_CHAR,P-1,tag,MPI_COMM_WORLD,&status);
		}
	}
		if (rank ==0) {
			writeFile(width*(height-2),decompressed_image,decompressed_file);
			done = compress(width, height-2, decompressed_file, compressed_file);

		}
		MPI_Barrier(MPI_COMM_WORLD); /* IMPORTANT */
    end = MPI_Wtime();
 		MPI_Finalize();
		if (rank == 0) { /* use time on master node */
    printf("Runtime = %f\n", end-start);
		}

		free(converted_image);
 		free(compressed_image);
 		free(decompressed_image);
 		free(decompressed_image_local);
		return 0;
}
