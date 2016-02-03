/*	
	N-Body simulation code.
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <mpi.h>

extern double	sqrt(double);
extern double	atan2(double, double);

#define GRAVITY		1.1
#define FRICTION	0.01
#define MAXBODIES	10000
#define DELTA_T		(0.025/5000)
#define	BOUNCE		-0.9
#define	SEED		27102015

typedef struct {
	double x[2];		/* Old and new X-axis coordinates */
	double y[2];		/* Old and new Y-axis coordinates */
	double xf;			/* force along X-axis */
	double yf;			/* force along Y-axis */
	double xv;			/* velocity along X-axis */
	double yv;			/* velocity along Y-axis */
	double mass;		/* Mass of the body */
	double radius;		/* width (derived from mass) */
} bodyType;

bodyType bodies[MAXBODIES];
bodyType bodies_collect[MAXBODIES];
double forces[2 * MAXBODIES];
double forces_collect[2 * MAXBODIES];
double position[2 * MAXBODIES];
double position_collect[2 * MAXBODIES];

int	bodyCt;
int	old = 0;	/* Flips between 0 and 1 */

/*	Macros to hide memory layout
*/
#define	X(B)		bodies[B].x[old]
#define	XN(B)		bodies[B].x[old^1]
#define	Y(B)		bodies[B].y[old]
#define	YN(B)		bodies[B].y[old^1]
#define	XF(B)		bodies[B].xf
#define	YF(B)		bodies[B].yf
#define	XV(B)		bodies[B].xv
#define	YV(B)		bodies[B].yv
#define	R(B)		bodies[B].radius
#define	M(B)		bodies[B].mass

/*	Dimensions of space (very finite, ain't it?)
*/
int		xdim = 0;
int		ydim = 0;

/* variables for MPI 
*/
int 	num_processors;
int 	my_id;
long 	offset;
long 	bodies_per_processor;

void 
check_MPI_error(int returnValue) 
{
	if ( returnValue != MPI_SUCCESS )
	{
		if ( my_id == 0 )
			fprintf(stderr,"MPI error");
		MPI_Finalize();
		exit(1);
	}
}

void
clear_forces(void)
{
	int b;

	/* Clear force accumulation variables */
	for (b=0; b<bodyCt; ++b) {
		forces[2 * b] = (forces[2 * b + 1] = 0);
	}
}

void
compute_forces(void)
{
	int b, c;
	/* Incrementally accumulate forces from each body pair,
	   skipping force of body on itself (c == b)
	*/
	for (b = offset; b < offset + bodies_per_processor; ++b){
		for (c = (b % 2); c < bodyCt; c = c + 2) {
			if (c == b)	{
				c++;
				if (!(c < bodyCt))	{
					break;
				}
			}
			double dx = X(c) - X(b);
			double dy = Y(c) - Y(b);
			double angle = atan2(dy, dx);
			double dsqr = dx*dx + dy*dy;
			double mindist = R(b) + R(c);
			double mindsqr = mindist*mindist;
			double forced = ((dsqr < mindsqr) ? mindsqr : dsqr);
			double force = M(b) * M(c) * GRAVITY / forced;
			double xf = force * cos(angle);
			double yf = force * sin(angle);

			/* Slightly sneaky...
			   force of b on c is negative of c on b;
			*/
			forces[2 * b] += xf;
			forces[2 * b + 1] += yf;
			forces[2 * c] -= xf;
			forces[2 * c + 1] -= yf;
		}
	}
	
	/* sum forces from all nodes and send this info to all nodes */
	check_MPI_error(MPI_Allreduce(&forces, &forces_collect, 2 * bodyCt, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD));

	for(b = 0; b < bodyCt; b++)
	{
		XF(b) = forces_collect[2 * b];
		YF(b) = forces_collect[2 * b + 1];
	}
}

void
compute_velocities(void)
{
	int b;

	/* let each node calculate velocity for their own bodies */
	for (b = offset; b < offset + bodies_per_processor; ++b) {
		double xv = XV(b);
		double yv = YV(b);
		double force = sqrt(xv*xv + yv*yv) * FRICTION;
		double angle = atan2(yv, xv);
		double xf = XF(b) - (force * cos(angle));
		double yf = YF(b) - (force * sin(angle));

		XV(b) += (xf / M(b)) * DELTA_T;
		YV(b) += (yf / M(b)) * DELTA_T;
	}
}

void
compute_positions(void)
{
	int b;

	/* let each node position velocity for their own bodies */
	for (b = offset; b < offset + bodies_per_processor; ++b) {
		double xn = X(b) + (XV(b) * DELTA_T);
		double yn = Y(b) + (YV(b) * DELTA_T);

		/* Bounce of image "walls" */
		if (xn < 0) {
			xn = 0;
			XV(b) = -XV(b);
		} else if (xn >= xdim) {
			xn = xdim - 1;
			XV(b) = -XV(b);
		}
		if (yn < 0) {
			yn = 0;
			YV(b) = -YV(b);
		} else if (yn >= ydim) {
			yn = ydim - 1;
			YV(b) = -YV(b);
		}

		/* Update position */
		position[2 * b] = xn;
		position[2 * b + 1] = yn;
	}

	/* Send all updated positions to each node */
	check_MPI_error(MPI_Allgather(&(position[2*offset]), bodies_per_processor * 2, MPI_DOUBLE, &position_collect, bodies_per_processor * 2, MPI_DOUBLE, MPI_COMM_WORLD));


	for(b = 0; b < bodyCt; b++)
	{
		XN(b) = position_collect[2 * b];
		YN(b) = position_collect[2 * b + 1];
	}
}

/*	Graphic output stuff...
*/

#include <fcntl.h>
#include <sys/mman.h>

int		fsize;
unsigned char	*map;
unsigned char	*image;


unsigned char *
map_P6(char *filename,
int *xdim,
int *ydim)
{
	/* The following is a fast and sloppy way to
	   map a color raw PPM (P6) image file
	*/
	int fd;
	unsigned char *p;
	int maxval;

	/* First, open the file... */
	if ((fd = open(filename, O_RDWR)) < 0) {
		return((unsigned char *) 0);
	}

	/* Read size and map the whole file... */
	fsize = lseek(fd, ((off_t) 0), SEEK_END);
	map = ((unsigned char *)
	       mmap(0,		/* Put it anywhere */
		    fsize,	/* Map the whole file */
		    (PROT_READ|PROT_WRITE),	/* Read/write */
		    MAP_SHARED,	/* Not just for me */
		    fd,		/* The file */
		    0));	/* Right from the start */
	if (map == ((unsigned char *) -1)) {
		close(fd);
		return((unsigned char *) 0);
	}

	/* File should now be mapped; read magic value */
	p = map;
	if (*(p++) != 'P') goto ppm_exit;
	switch (*(p++)) {
	case '6':
		break;
	default:
		goto ppm_exit;
	}

#define	Eat_Space \
	while ((*p == ' ') || \
	       (*p == '\t') || \
	       (*p == '\n') || \
	       (*p == '\r') || \
	       (*p == '#')) { \
		if (*p == '#') while (*(++p) != '\n') ; \
		++p; \
	}

	Eat_Space;		/* Eat white space and comments */

#define	Get_Number(n) \
	{ \
		int charval = *p; \
 \
		if ((charval < '0') || (charval > '9')) goto ppm_exit; \
 \
		n = (charval - '0'); \
		charval = *(++p); \
		while ((charval >= '0') && (charval <= '9')) { \
			n *= 10; \
			n += (charval - '0'); \
			charval = *(++p); \
		} \
	}

	Get_Number(*xdim);	/* Get image width */

	Eat_Space;		/* Eat white space and comments */
	Get_Number(*ydim);	/* Get image width */

	Eat_Space;		/* Eat white space and comments */
	Get_Number(maxval);	/* Get image max value */

	/* Should be 8-bit binary after one whitespace char... */
	if (maxval > 255) {
ppm_exit:
		close(fd);
		munmap(map, fsize);
		return((unsigned char *) 0);
	}
	if ((*p != ' ') &&
	    (*p != '\t') &&
	    (*p != '\n') &&
	    (*p != '\r')) goto ppm_exit;

	/* Here we are... next byte begins the 24-bit data */
	return(p + 1);

	/* Notice that we never clean-up after this:

	   close(fd);
	   munmap(map, fsize);

	   However, this is relatively harmless;
	   they will go away when this process dies.
	*/
}

#undef	Eat_Space
#undef	Get_Number

static inline void
color(int x, int y, int b)
{
	unsigned char *p = image + (3 * (x + (y * xdim)));
	int tint = ((0xfff * (b + 1)) / (bodyCt + 2));

	p[0] = (tint & 0xf) << 4;
	p[1] = (tint & 0xf0);
	p[2] = (tint & 0xf00) >> 4;
}

static inline void
black(int x, int y)
{
	unsigned char *p = image + (3 * (x + (y * xdim)));

	p[2] = (p[1] = (p[0] = 0));
}

void
display(void)
{
	double i, j;
	int b;

	/* For each pixel */
	for (j=0; j<ydim; ++j) {
		for (i=0; i<xdim; ++i) {
			/* Find the first body covering here */
			for (b=0; b<bodyCt; ++b) {
				double dy = Y(b) - j;
				double dx = X(b) - i;
				double d = sqrt(dx*dx + dy*dy);

				if (d <= R(b)+0.5) {
					/* This is it */
					color(i, j, b);
					goto colored;
				}
			}

			/* No object -- empty space */
			black(i, j);

colored:		;
		}
	}
}

void
print(void)
{
	int b;

	for (b=0; b<bodyCt; ++b) {
		printf("%10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n", X(b), Y(b), XF(b), YF(b), XV(b), YV(b));
	}
}

/*	Main program...
*/

int
main(int argc, char **argv)
{
	unsigned int lastup = 0;
	unsigned int secsup;
	int b;
	int steps;
	double rtime;
	double start_time, end_time;

	/* Get Parameters */
	if (argc != 5) {
		fprintf(stderr,
			"Usage: %s num_bodies secs_per_update ppm_output_file steps\n",
			argv[0]);
		exit(1);
	}
	if ((bodyCt = atol(argv[1])) > MAXBODIES ) {
		fprintf(stderr, "Using only %d bodies...\n", MAXBODIES);
		bodyCt = MAXBODIES;
	} else if (bodyCt < 2) {
		fprintf(stderr, "Using two bodies...\n");
		bodyCt = 2;
	}

	/* initialize MPI */
	check_MPI_error(MPI_Init(&argc, &argv));
    check_MPI_error(MPI_Comm_size(MPI_COMM_WORLD, &num_processors));
    check_MPI_error(MPI_Comm_rank(MPI_COMM_WORLD, &my_id));

    /* store arguments */
	secsup = atoi(argv[2]);
	image = map_P6(argv[3], &xdim, &ydim);
	steps = atoi(argv[4]);

	/* Initialize simulation data */
	if(my_id == 0) {
		fprintf(stderr, "Running N-body with %i bodies and %i steps and %i processors\n", bodyCt, steps, num_processors);

		srand(SEED);
		for (b=0; b<bodyCt; ++b) {
			X(b) = (rand() % xdim);
			Y(b) = (rand() % ydim);
			R(b) = 1 + ((b*b + 1.0) * sqrt(1.0 * ((xdim * xdim) + (ydim * ydim)))) /
				(25.0 * (bodyCt*bodyCt + 1.0));
			M(b) = R(b) * R(b) * R(b);
			XV(b) = ((rand() % 20000) - 10000) / 2000.0;
			YV(b) = ((rand() % 20000) - 10000) / 2000.0;
		}
	}

	fprintf(stderr, "my process id: %d \n", my_id);

	/* Get the range of bodies for the current process */
	bodies_per_processor = ceil((double)bodyCt / (double)num_processors);
	offset = bodies_per_processor * my_id;
	/* range of bodies for last processor */
	if ((bodies_per_processor + offset) > bodyCt) {
		bodies_per_processor = bodyCt - offset;
	}

	/* Broadcast all bodies to all processors*/
	check_MPI_error(MPI_Bcast(bodies, 10 * bodyCt, MPI_DOUBLE, 0, MPI_COMM_WORLD ));

	start_time = MPI_Wtime();

	/* Main Loop */
	while (steps--) {
		clear_forces();
		compute_forces();
		compute_velocities();
		compute_positions();

		/* Flip old & new coordinates */
		old ^= 1;

		/*Time for a display update?*/ 
		if (secsup > 0 && (time(0) - lastup) > secsup) {
			display();
			msync(map, fsize, MS_SYNC);	/* Force write */
			lastup = time(0);
		}
	}
	
	end_time = MPI_Wtime();
	
	/* Root machine collects all body information from other machines. */
	check_MPI_error(MPI_Gather(&(bodies[offset]), 10 * bodies_per_processor, MPI_DOUBLE, &bodies_collect, 10 * bodies_per_processor, MPI_DOUBLE, 0, MPI_COMM_WORLD));

	for (b = 0; b < bodyCt; b++) {
		XV(b) = bodies_collect[b].xv;
		YV(b) = bodies_collect[b].yv;
	}

	rtime = end_time - start_time;
	double longest_time;

	/* root node gets time of the machine that took the longest */
	check_MPI_error(MPI_Reduce(&rtime, &longest_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD));

	if ( my_id == 0 ) {
		fprintf(stderr, "N-body took %10.3f seconds\n", longest_time);
		print();
	}

	MPI_Finalize();

	return 0;
}

