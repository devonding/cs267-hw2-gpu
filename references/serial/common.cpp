#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include "common.h"



double size = 0.0;
int num = 0;
double dim = 0.0;
//
//  tuned constants
//
#define density 0.0005
#define mass    0.01
#define cutoff  0.01
#define min_r   (cutoff/100)
#define dt      0.0005

//
//  timer
//
double read_timer()
{
    static bool initialized = false;
    static struct timeval start;
    struct timeval end;
    if( !initialized )
    {
        gettimeofday( &start, NULL );
        initialized = true;
    }
    gettimeofday( &end, NULL );
    return (end.tv_sec - start.tv_sec) + 1.0e-6 * (end.tv_usec - start.tv_usec);
}

//
//  keep density constant
//  
void set_size( int n )
{
    size = sqrt( density * n );
}

//
//  Initialize the particle positions and velocities
//  in this function, it will decide the location of the each partical
//  and their init speed
void init_particles( int n, particle_t *p )
{
    srand48( time( NULL ) );
    
    // init the start location, i.e. center location

    int sx = (int)ceil(sqrt((double)n));
    int sy = (n+sx-1)/sx;
    
    int *shuffle = (int*)malloc( n * sizeof(int) );
    for( int i = 0; i < n; i++ )
        shuffle[i] = i;
    
    for( int i = 0; i < n; i++ ) 
    {
        //
        //  make sure particles are not spatially sorted
        //  this is like to generate a random permutation

        int j = lrand48()%(n-i);
        int k = shuffle[j];  
        shuffle[j] = shuffle[n-i-1];
        //  put the chosen id to the end of the array
        //  thus, k is random sampling from 0 to n - 1

        //
        //  distribute particles evenly to ensure proper spacing
        //
        // here sx and sy are pre-determined the location,
        // then we randomly put all the particals into n space
        p[i].x = size*(1.+(k%sx))/(1+sx);
        p[i].y = size*(1.+(k/sx))/(1+sy);

        //
        //  assign random velocities within a bound
        //

        // the speed is random
        p[i].vx = drand48()*2-1;
        p[i].vy = drand48()*2-1;

        // the region of grid will be inited at grid.init
    }
    free( shuffle );
}

//
//  interact two particles
//  this function would change the acc of the particle
//  here dmin, keep updating is the min dist for particle;
//  davg*1.0 / navg is the avg dist for the particle, also keep updating

void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg)
{
    // get the difference at the location between two particles

    double dx = neighbor.x - particle.x;
    double dy = neighbor.y - particle.y;

    //printf("neighbor is: x: %f ; y: %f \n",neighbor.x, neighbor.y);
    //printf("self is: x: %f ; y: %f \n", particle.x, particle.y);
    
    double r2 = dx * dx + dy * dy;
    //printf("r2 is %f\n", r2);
    //printf("the init a is ax: %f; ay %f\n", particle.ax, particle.ay );

    if( r2 > cutoff*cutoff )
        return; // the force only exists inside the cutoff
	if (r2 != 0)
        {
            //printf("we have force now..\n");
            // here dmin is like a global var, we keep update it..
            // we also keep update davg and navg.. davg is the sum of all the dist
            // the navg is the count of the dist, thus, combining them we can get the avg
	   if (r2/(cutoff*cutoff) < *dmin * (*dmin))
	      *dmin = sqrt(r2)/cutoff;
           (*davg) += sqrt(r2)/cutoff;
           (*navg) ++;
        }
		
    r2 = fmax( r2, min_r*min_r );
    double r = sqrt( r2 );
 
    
	
    //
    //  very simple short-range repulsive force
    //
    double coef = ( 1 - cutoff / r ) / r2 / mass;
    //printf("a:%f; r2: %f \n",coef, r2 );
    particle.ax += coef * dx;
    particle.ay += coef * dy;
}

//
//  integrate the ODE
//
void move( particle_t &p )
{
    //
    //  slightly simplified Velocity Verlet integration
    //  conserves energy better than explicit Euler method
    //
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    p.x  += p.vx * dt;
    p.y  += p.vy * dt;

    //
    //  bounce from walls
    //

    // just do the refelction from the wall
    while( p.x < 0 || p.x > size )
    {
        p.x  = p.x < 0 ? -p.x : 2*size-p.x;
        p.vx = -p.vx;
    }
    while( p.y < 0 || p.y > size )
    {
        p.y  = p.y < 0 ? -p.y : 2*size-p.y;
        p.vy = -p.vy;
    }
}

void moveWithGrid(particle_t &p, subgrid *gridContainer) {
    
    p.vx += p.ax * dt;
    p.vy += p.ay * dt;
    double tmpx = p.x + p.vx * dt;
    double tmpy = p.y + p.vy * dt;
    //printf("the speed is vx: %f, vy %f; the a is ax: %f; ay: %f\n", p.ax, p.ay, p.vx, p.vy);
    
    // here we need to check whether the particle is inside the same subgrid
    
    //
    //  bounce from walls
    //

    // just do the refelction from the wall
    while( tmpx < 0 || tmpx > size )
    {
        tmpx  = tmpx < 0 ? -tmpx : 2*size-tmpx;
        p.vx = -p.vx;
    }
    while( tmpy < 0 || tmpy > size )
    {
        tmpy  = tmpy < 0 ? -tmpy : 2*size-tmpy;
        p.vy = -p.vy;
    }

    int id_old = locationToID(p); // this is the old id which the particle located
    int id_new = locationToID(tmpx, tmpy); // this is the new id ..

    // the id here just means the index in the gridContainer array

    if (id_old != id_new) {
        // then we need to move it..
        gridContainer[id_old].pArray.erase(&p);
        gridContainer[id_old].count --;

        gridContainer[id_new].pArray.insert(&p);
        gridContainer[id_new].count ++;
    }
    p.x = tmpx;
    p.y = tmpy;
}

//
//  I/O routines
//
void save( FILE *f, int n, particle_t *p )
{
    static bool first = true; // whether to print the head of the file
    if( first )
    {
        fprintf( f, "%d %g\n", n, size );
        first = false;
    }
    // print the location of the each partical
    for( int i = 0; i < n; i++ )
        fprintf( f, "%g %g\n", p[i].x, p[i].y );
}

//
//  command line option processing
//
int find_option( int argc, char **argv, const char *option )
{
    for( int i = 1; i < argc; i++ )
        if( strcmp( argv[i], option ) == 0 )
            return i;
    return -1;
}

int read_int( int argc, char **argv, const char *option, int default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return atoi( argv[iplace+1] );
    return default_value;
}

char *read_string( int argc, char **argv, const char *option, char *default_value )
{
    int iplace = find_option( argc, argv, option );
    if( iplace >= 0 && iplace < argc-1 )
        return argv[iplace+1];
    return default_value;
}
