#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "seamcarving.h"
#include "c_img.h"

#define MIN(x,y) (((x) < (y)) ? (x) : (y))
/* in other words, can also mean:
int MIN(x, y);

if (x < y){
    MIN(x,y) = x;
} else if (x > y){
    MIN(x,y) = y;
} */

void calc_energy(struct rgb_img *im, struct rgb_img **grad){
    //computes the dual-gradient energy function and place it in the struct rgb_img *grad

    create_img(grad, im->height, im->width); //Allocate memory to store seam energies
    //HEIGHT, WIDTH, AND *RASTER (PART OF struct rgb_img) CAN BE FOUND IN c_img.h

    //Initializing r(x, y), g(x, y), b(x, y).

    //the defined coordinate system is (y, x) for the purpose of making calculations easier
    //an alternate coordinate system would be to use unit vectors, x == i, y == j
    int r_x, r_y, g_x, g_y, b_x, b_y;
    int d_x, d_y;
    int pixel_energy, grad_energy;

    for (int y = 0; y < im->height; y++){
        for (int x = 0; x < im->width; x++){
            //left_img
            //get_pixel used from c_img.h || 4 inputs
            if (x == 0){
                r_x = get_pixel(im, y, (im->width) - 1, 0) - get_pixel(im, y, x + 1, 0);
                g_x = get_pixel(im, y, (im->width) - 1, 1) - get_pixel(im, y, x + 1, 1);
                b_x = get_pixel(im, y, (im->width) - 1, 2) - get_pixel(im, y, x + 1, 2);

            //right_img
            } else if (x == (im->width)-1){
                r_x = get_pixel(im, y, x - 1, 0) - get_pixel(im, y, 0, 0);
                g_x = get_pixel(im, y, x - 1, 1) - get_pixel(im, y, 0, 1);
                b_x = get_pixel(im, y, x - 1, 2) - get_pixel(im, y, 0, 2);
            } else{
                r_x = get_pixel(im, y, x - 1, 0) - get_pixel(im, y, x + 1, 0);
                g_x = get_pixel(im, y, x - 1, 1) - get_pixel(im, y, x + 1, 1);
                b_x = get_pixel(im, y, x - 1, 2) - get_pixel(im, y, x + 1, 2);
            }

            //top_img
            if (y == 0){
                r_y = get_pixel(im, (im->height) - 1, x, 0) - get_pixel(im, y + 1, x, 0);
                g_y = get_pixel(im, (im->height) - 1, x, 1) - get_pixel(im, y + 1, x, 1);
                b_y = get_pixel(im, (im->height) - 1, x, 2) - get_pixel(im, y + 1, x, 2);
            //bot_img
            } else if (y == (im->height) - 1){
                r_y = get_pixel(im, y - 1, x, 0) - get_pixel(im, 0, x, 0);
                g_y = get_pixel(im, y - 1, x, 1) - get_pixel(im, 0, x, 1);
                b_y = get_pixel(im, y - 1, x, 2) - get_pixel(im, 0, x, 2);
            } else{
                r_y = get_pixel(im, y - 1, x, 0) - get_pixel(im, y + 1, x, 0);
                g_y = get_pixel(im, y - 1, x, 1) - get_pixel(im, y + 1, x, 1);
                b_y = get_pixel(im, y - 1, x, 2) - get_pixel(im, y + 1, x, 2);
            }

            d_x = (r_x * r_x) + (g_x * g_x) + (b_x * b_x); //differences in rgb components along the x-axis
            d_y = (r_y * r_y) + (g_y * g_y) + (b_y * b_y); //differences in rgb components along the y-axis

            //energy of pixel (y, x) is
            pixel_energy = sqrt(d_x + d_y);

            //for each pixel, set the r, g, and b channels to the same value (the energy divided by 10 and cast to uint8_t).
            grad_energy = (uint8_t)(pixel_energy / 10);

            //AFTER CALCULATING FOR PIXEL_ENERGY AND GRADIENT_ENERGY
            //set_pixel used from c_img.h || 6 inputs
            set_pixel(*grad, y, x, grad_energy, grad_energy, grad_energy);
        }
    }
}

void dynamic_seam(struct rgb_img *grad, double **best_arr){
    //allocates and computes the dynamic array *best_arr

    *best_arr = malloc((grad->height) * (grad->width) * sizeof(double) + 1); //make room for null pointer: for an array of length n, array has to be length n + 1, 1 space for NULL POINTER
    int width = grad->width, height = grad->height;

    //(*best_arr)[i*width+j] contains minimum cost of a seam from the top of grad to the point P(i, j).
    for (int i = 0; i < width; i++){
        (*best_arr)[i] = grad->raster[3 * i];
    }
    for (int i = 1; i < height; i++){
        for (int j = 0; j < width; j++){
            if (j == 0){
                (*best_arr)[i * width + j] = (double)grad->raster[3 * (i * (grad->width) + j)] + (double)MIN((*best_arr)[(i - 1) * width + j], (*best_arr)[(i - 1) * width + (j + 1)]);
            } else if (j == (width - 1)){
                (*best_arr)[i * width + j] = (double)grad->raster[3 * (i * (grad->width) + j)] + (double)MIN((*best_arr)[(i - 1) * width + j], (*best_arr)[(i - 1) * width + (j - 1)]);
            } else{
                double min = MIN((*best_arr)[(i - 1) * width + j], (*best_arr)[(i - 1) * width + (j - 1)]);
                (*best_arr)[i * width + j] = (double)grad->raster[3 * (i * (grad->width) + j)] + (double)MIN(min, (*best_arr)[(i - 1) * width + (j + 1)]);   
            }
        }
    }
}

void recover_path(double *best, int height, int width, int **path){
    //allocates a path through the minimum seam as defined by the array best

    double min = 10000000000;
    int j = height, min_index = -1;
    *path = malloc(j * sizeof(int) + 1); //make room for NULL POINTER

    for (int i = 0; i < width; i++){
        if (best[(j - 1) * width + i] < min){
            min = best[(j - 1) * width + i];
            min_index = i;
        }
    }

    (*path)[height - 1] = min_index;
    for (int i = j - 1; i > 0; i--){
        //Find the min of the three/two values that are above it
        
        if (min_index == 0){
            min = MIN(best[(i - 1) * width], best[((i - 1) * width) + 1]); 
            if (min == best[(i - 1) * width]){
                min_index = 0;
            } else{
                min_index = 1;
            }
        } else if (min_index == width - 1){
            min = MIN(best[(i - 1) * width + min_index], best[(i - 1) * width + min_index - 1]);
            if (min == best[(i - 1) * width + min_index]){
                min_index = width - 1;
            } else{
                min_index = width - 2;
            }
        } else{
            double min_temporary = MIN(best[(i - 1) * width + min_index], best[(i - 1) * width + (min_index - 1)]);
            min = MIN(best[(i - 1) * width + (min_index + 1)], min_temporary);
            if (min == best[(i - 1) * width + min_index - 1]){
                min_index--; //min_index = min_index - 1
            } if (min == best[(i - 1) * width + (min_index)]){
                min_index = min_index;
            } else{
                min_index++; //min_index = min_index + 1
            }
        }
        (*path)[i - 1] = min_index;
    }
}

void remove_seam(struct rgb_img *src, struct rgb_img **dest, int *path){
    //creates the destination image, and writes to it the source image, with the seam removed
    
    //create_image from c_img.h
    create_img(dest, src->height, src->width - 1);
    //Initialize rgb
    int r, g, b;

    for (int y = 0; y < src->height; y++){
        for (int x = 0; x < src->width; x++){
            if (x == path[y]){
                continue;
            } else{
                r, g, b = get_pixel(src, y, x, 0), get_pixel(src, y, x, 1), get_pixel(src, y, x, 2);
                /* alternatively,
                r = get_pixel(src, y, x, 0);
                g = get_pixel(src, y, x, 1);
                b = get_pixel(src, y, x, 2); 
                */
                if (x < path[y]){
                    set_pixel(*dest, y, x, r, g, b);
                } else{
                    set_pixel(*dest, y, x - 1, r, g, b);
                }
            }
        }
    }
}


