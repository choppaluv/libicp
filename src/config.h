/* This is a file for ICP configuration */

/* # of sample points for ICP */
#define NUM_SAMPLE_POINTS 10000
/* dimension of the points */
#define DIM 3
/* # of max iteration for ICP */
#define MAX_ITER 200
/* threshold for ICP convergence */
#define THRESHOLD 1e-4
/* error metric used in visualization */
#define VISUALIZE_POINT_TO_POINT true
#define VISUALIZE_POINT_TO_PLANE false
#define VISUALIZE_ERROR_METRIC VISUALIZE_POINT_TO_PLANE
/* scale for color map error */
#define ERROR_SCALE 4