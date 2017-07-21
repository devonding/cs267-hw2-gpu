#ifndef __CS267_COMMON_H__
#define __CS267_COMMON_H__
#include <set>
#include <omp.h>

/** this is the head file for both common.cpp and  grid.cpp **/



inline int min( int a, int b ) { return a < b ? a : b; }
inline int max( int a, int b ) { return a > b ? a : b; }

//
//  saving parameters
//
const int NSTEPS = 10000;
const int SAVEFREQ = 10;

// self added
const double CUTOFF = 0.01;
extern double size;
// I change the location of size from the common.cpp to common.h
// since it's needed in grid

// particle data structure
//
typedef struct
{
  double x;
  double y;
  double vx;
  double vy;
  double ax;
  double ay;
} particle_t;

// I add my struct, the struct for a subgrid, the whole grid is the array of the subgrid

typedef struct
{
  double x;
  double y;
  // the rightup pnt of the grid
  // to hold the particle's pnter, this is pnt to a pnt who pnt to particle_t
  std::set<particle_t*> pArray;
  int count; // record the num of the particles in this struct
  //int space;
  omp_lock_t lock;
} subgrid;
//
//  timing routines
//
double read_timer( );

//
//  simulation routines
//
void set_size( int n );
void init_particles( int n, particle_t *p );
void apply_force( particle_t &particle, particle_t &neighbor , double *dmin, double *davg, int *navg);
void move( particle_t &p );

// I added, this move function will check the location of each partical and reput them into the correct
// place if they change the location
void moveWithGrid(particle_t &p, subgrid *gridContainer);
//
//  I/O routines
//
FILE *open_save( char *filename, int n );
void save( FILE *f, int n, particle_t *p );

//
//  argument processing routines
//
int find_option( int argc, char **argv, const char *option );
int read_int( int argc, char **argv, const char *option, int default_value );
char *read_string( int argc, char **argv, const char *option, char *default_value );



const double DIM = 2.0 * CUTOFF; // DIM is the diermitaro of the subgrid
extern int num;   // the num of the subgrid for one edge, thus, the tot num = num ^2
extern double dim;  // the acutally dim of the size

//---------------------------------------------------------------------//
//function for the grid.cpp

// init the grid
void init_grid();

// judge whether a particle is inside a sugrid
bool isIngrid(particle_t &particle, subgrid &grid);

// put all the particles (n) into the whole grid
void fill_grid(particle_t *particles, int n);

// update the particles, including update their force, and then move them
void update_grid(double *dmin, double *davg, int *navg, particle_t *particles, int n);

// help method to get the id of the subgrid in the subgrid array with given the location of a paritcle
// x, and y are the location of a particle
int locationToID(double x, double y);
// same fun, but with input as a particle
int locationToID(particle_t &particle);

// delete and free the space
void dele_grid();


/*-------------------------------------------------------*/
// debug funcs, just print the set
void debug_printSet(std::set<particle_t*> s);


#endif
