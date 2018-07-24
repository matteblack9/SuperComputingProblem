//
//  main.c
//  KSC2017 Tool
//
//  Created by Chan Park on 2017. 6. 20..
//  Copyright © 2017년 KISTI. All rights reserved.
//

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include <png.h>

enum _RAY_STATUS {HIT, FALL, OUTSIDE, YET, C_NUMBER_OF_STATUS};

png_byte grey[3] = {30, 30, 30};
png_byte yellow[3] = {248, 202, 71};

void error(const char *message)
{
    printf("%s\n", message);
    abort();
}

bool read_data(const char *filename, int *width, int *height, int **status, double **y, double **z)
{
    int i;
    int total_pixel;
    FILE *data_fp = fopen(filename, "rb");
    if (NULL == data_fp) return false;
    if (1 != fread(width  , sizeof(int), 1, data_fp)) return false;
    if (1 != fread(height , sizeof(int), 1, data_fp)) return false;
    if (*width <= 0 || *height <= 0) return false;
    total_pixel = (*width) * (*height);
    *status = malloc(sizeof(int)*total_pixel);
    *y = malloc(sizeof(double)*total_pixel);
    *z = malloc(sizeof(double)*total_pixel);
    if (total_pixel != fread(*status, sizeof(   int), total_pixel, data_fp)) return false;
    if (total_pixel != fread(     *y, sizeof(double), total_pixel, data_fp)) return false;
    if (total_pixel != fread(     *z, sizeof(double), total_pixel, data_fp)) return false;
    fclose(data_fp);
    
    for (i = 0; i < total_pixel; i++) {
        if ((*status)[i] < 0 || (*status)[i] >= C_NUMBER_OF_STATUS) return false;
        if ((*status)[i] == HIT && (isnan((*y)[i]) || isnan((*z)[i]))) return false;
    }
    
    return true;
}

void image_tool(const char *argv[])
{
    int i, j;
    
    const char *data_filename = argv[1];
    int *status;
    double *y;
    double *z;
    
    const char *import_image_filename = argv[2];
    FILE *import_fp;
    png_structp import_png;
    png_infop import_info;
    png_byte import_color_type;
    png_byte import_bit_depth;
    png_bytep *import_row_pointers;
    int import_width, import_height, import_size;
    int import_w, import_h;
    
    const char *export_image_filename = argv[3];
    FILE *export_fp;
    png_structp export_png;
    png_infop export_info;
    png_bytep *export_row_pointers;
    int export_width, export_height;
    int export_w, export_h;
    
    /* Import process */
    import_fp = fopen(import_image_filename, "rb");
    if (!import_fp) error("Import image is not available.");
    
    import_png = png_create_read_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!import_png) error("Import image is not available.");
    
    import_info = png_create_info_struct(import_png);
    if (!import_info) error("Import image is not available.");
    
    if(setjmp(png_jmpbuf(import_png))) error("Import image is not available.");;
    
    png_init_io(import_png, import_fp);
    
    png_read_info(import_png, import_info);
    
    import_width        = png_get_image_width   (import_png, import_info);
    import_height       = png_get_image_height  (import_png, import_info);
    import_color_type   = png_get_color_type    (import_png, import_info);
    import_bit_depth    = png_get_bit_depth     (import_png, import_info);
    
    /* Read any color_type into 8bit depth, RGBA format. */
    /* See http://www.libpng.org/pub/png/libpng-manual.txt */
    
    if (import_bit_depth == 16) {
        png_set_strip_16(import_png);
    }
    
    if (import_color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_palette_to_rgb(import_png);
    }
    
    /* PNG_COLOR_TYPE_GRAY_ALPHA is always 8 or 16bit depth. */
    if (import_color_type == PNG_COLOR_TYPE_GRAY && import_bit_depth < 8) {
        png_set_expand_gray_1_2_4_to_8(import_png);
    }
    
    if(png_get_valid(import_png, import_info, PNG_INFO_tRNS)) {
        png_set_tRNS_to_alpha(import_png);
    }
    
    /* These color_type don't have an alpha channel then fill it with 0xff. */
    if (import_color_type == PNG_COLOR_TYPE_RGB || import_color_type == PNG_COLOR_TYPE_GRAY || import_color_type == PNG_COLOR_TYPE_PALETTE) {
        png_set_filler(import_png, 0xFF, PNG_FILLER_AFTER);
    }
    
    if (import_color_type == PNG_COLOR_TYPE_GRAY || import_color_type == PNG_COLOR_TYPE_GRAY_ALPHA) {
        png_set_gray_to_rgb(import_png);
    }
    
    png_read_update_info(import_png, import_info);
    
    import_row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * import_height);
    for (import_h = 0; import_h < import_height; import_h++) {
        import_row_pointers[import_h] = (png_byte*)malloc(png_get_rowbytes(import_png, import_info));
    }
    
    png_read_image(import_png, import_row_pointers);
    
    fclose(import_fp);
    
    
    
    if (!read_data(data_filename, &export_width, &export_height, &status, &y, &z)) error("Data file is not available.");
    
    
    /* Export preparation */
    export_fp = fopen(export_image_filename, "wb");
    if(!export_fp) error("Export image is not available.");
    
    export_png = png_create_write_struct(PNG_LIBPNG_VER_STRING, NULL, NULL, NULL);
    if (!export_png) error("Export image is not available.");
    
    export_info = png_create_info_struct(export_png);
    if (!export_info) error("Export image is not available.");
    
    if (setjmp(png_jmpbuf(export_png))) error("Export image is not available.");
    
    png_init_io(export_png, export_fp);
    
    /* Output is 8bit depth, RGBA format. */
    png_set_IHDR(
                 export_png,
                 export_info,
                 export_width, export_height,
                 8,
                 PNG_COLOR_TYPE_RGB,
                 PNG_INTERLACE_NONE,
                 PNG_COMPRESSION_TYPE_DEFAULT,
                 PNG_FILTER_TYPE_DEFAULT
                 );
    png_write_info(export_png, export_info);
    
    /* allocation of export array */
    export_row_pointers = (png_bytep*)malloc(sizeof(png_bytep) * export_height);
    for (export_h = 0; export_h < export_height; export_h++) {
        export_row_pointers[export_h] = (png_byte*)malloc(png_get_rowbytes(export_png, export_info));
    }
    
    /* initialize */
    for (export_h = 0; export_h < export_height; export_h++) {
        for (export_w = 0; export_w < export_width; export_w++) {
            for (i = 0; i < 3; i++) {
                export_row_pointers[export_h][3 * export_w + i] = grey[i];
            }
        }
    }
    
    import_size = (import_height < import_width) ? import_height : import_width;
    for (i = 0; i < export_width * export_height; i++) {
        export_h = i / export_width;
        export_w = i % export_width;
        if (HIT == status[i]) {
            import_w = (int)round(+ y[i] * (import_size - 1) + (import_width  - 1) / 2.);
            import_h = (int)round(- z[i] * (import_size - 1) + (import_height - 1) / 2.);
            for (j = 0; j < 3; j++) {
                export_row_pointers[export_h][3 * export_w + j] = import_row_pointers[import_h][4 * import_w + j];
            }
        } else if (FALL == status[i]) {
            for (j = 0; j < 3; j++) {
                export_row_pointers[export_h][3 * export_w + j] = 0;
            }
        }
    }
    
    /* Write process */
    png_write_image(export_png, export_row_pointers);
    png_write_end(export_png, NULL);
    
    fclose(export_fp);
    
    for(import_h = 0; import_h < import_height; import_h++) {
        free(import_row_pointers[import_h]);
    }
    free(import_row_pointers);
    
    free(status);
    free(y);
    free(z);
    
    for(export_h = 0; export_h < export_height; export_h++) {
        free(export_row_pointers[export_h]);
    }
    free(export_row_pointers);
}

void verify_tool(const char *argv[])
{
    int width, height;
    int *status;
    double *y;
    double *z;
    
    if (read_data(argv[1], &width, &height, &status, &y, &z)) {
        printf("No problem.\n");
    } else {
        printf("Data file is not available.\n");
    }
    
    free(status);
    free(y);
    free(z);
}


#define DEALLOC(p) \
    if (NULL != (p)) free(p);

#define FINALIZE \
    DEALLOC(status1) \
    DEALLOC(y1) \
    DEALLOC(z1) \
    DEALLOC(status2) \
    DEALLOC(y2) \
    DEALLOC(z2) \

#define EXIT(msg) \
    printf(msg"\n"); \
    FINALIZE \
    return;


#define SQ(a) ((a)*(a))

void diff_tool(const char *argv[])
{
    int width1, height1, *status1 = NULL;
    double *y1 = NULL, *z1 = NULL;
    int width2, height2, *status2 = NULL;
    double *y2 = NULL, *z2 = NULL;
    
    if (!read_data(argv[1], &width1, &height1, &status1, &y1, &z1)) { EXIT("Data file 1 is not available.") }
    if (!read_data(argv[2], &width2, &height2, &status2, &y2, &z2)) { EXIT("Data file 2 is not available.") }
    if (width1 != width2) { EXIT("Widths are different.") }
    if (height1 != height2) { EXIT("Heights are different.") }
    
    int number_of_pixels = width1 * height1;
    int number_of_different_status = 0;
    for (int i = 0; i < number_of_pixels; i++) {
        if (status1[i] != status2[i]) number_of_different_status++;
    }
    printf("different status pixels / total pixels : %e\n", ((double)number_of_different_status) / number_of_pixels);
    
    double sum_y = 0., sum_z = 0.;
    int number_of_hit_pixels = 0;
    for (int i = 0; i < number_of_pixels; i++) {
        if (status1[i] == HIT && status2[i] == HIT) {
            sum_y += SQ(y1[i] - y2[i]);
            sum_z += SQ(z1[i] - z2[i]);
            number_of_hit_pixels++;
        }
    }
    printf("Root mean square of deviation : y(%e), z(%e)\n", sqrt(sum_y / number_of_hit_pixels), sqrt(sum_z / number_of_hit_pixels));
    
    EXIT("")
}

int main(int argc, const char *argv[]) {
    
    switch(argc) {
        case 2:
            verify_tool(argv);
            break;
        case 3:
            diff_tool(argv);
            break;
        case 4:
            image_tool(argv);
            break;
        default:
            printf("[Data Verification Tool]\nUsage: %s data_file\n\n[Data Diff Tool]\nUsage: %s data_file1 data_file2\n\n[Image Making Tool]\nUsage: %s data_file import_png_image export_png_image\n", argv[0], argv[0], argv[0]);
    }
    
    return 0;
}
