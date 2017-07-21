#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <set>
#include "common.h"


using namespace std;

void init_grid();

int locationToID(particle_t &particle);

int locationToID(double x, double y);


bool isIngrid(particle_t &particle, subgrid &grid);

void fill_grid(particle_t *particles, int n);

void apply_forces(double *dmin, double *davg, int *navg, particle_t *particles, int n, int i);

void move_particles(particle_t *particles, int n, int i);

void dele_grid();


void debug_printSet(std::set<particle_t*> s);
