#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <set>
#include "openmp_common.h"
#include <omp.h>


using namespace std;


// this is the pnter of the subgrid array
// the space will be allocate in the init_grid()
subgrid *gridContainer;



// here we init the grid according to the size
void init_grid() {

  num = ceil(size*1.0 / DIM); // we get the num of the grid for one directions
  dim = size/num; // the acutal size of a subgrid

  gridContainer = new subgrid[num * num];
  // then use a for loop to init the gridContainer

  int i = 0;
  for (i = 0; i < num * num; i++) {
    gridContainer[i].x = i / num * dim;
    gridContainer[i].y = i % num * dim;
    gridContainer[i].count = 0;
    omp_init_lock(&gridContainer[i].lock);
  }

  // just for debug...
  //for (i = 0; i < num * num; i++) {
          //printf("the x of the %d, is %f; y is %f\n", i, gridContainer[i].x, gridContainer[i].y );
  //}
}

int locationToID(particle_t &particle) {
  int xID = particle.x / dim;
  int yID = particle.y / dim;
  return xID * num + yID;
}

int locationToID(double x, double y) {
  int xID = x / dim;
  int yID = y / dim;
  return xID * num + yID;
}


bool isIngrid(particle_t &particle, subgrid &grid) {
  if (particle.x >= grid.x && particle.y >= grid.y && particle.x <= (grid.x + dim) &&
   particle.y <= (grid.y + dim)) return true;
  return false;
}

// add the pnt to the add of the particles to the each subgrid....
// n is the num of the particles in the system
void fill_grid(particle_t *particles, int n) {
  int i = 0;
  for (i = 0; i < n; i ++) {
    // just call the locationToID to get the id of the grid
    int id = locationToID(particles[i]);
    //printf("id is %d\n", id);
    gridContainer[id].pArray.insert(particles + i);
    gridContainer[id].count++;

  }
}

// update the grid

void apply_forces(double *dmin, double *davg, int *navg, particle_t *particles, int n, int i){
  // here we just need to check the 9 subgrids
  std::set<particle_t*>::iterator it;
  it = gridContainer[i].pArray.begin();
  for (; it != gridContainer[i].pArray.end(); ++ it) {
    (**it).ax = 0.0;
    (**it).ay = 0.0;

    int x_s = max(i/num - 1, 0);
    int x_e = min(i/num + 1, num - 1);
    int y_s = max(i%num - 1, 0);
    int y_e = min(i%num + 1, num - 1);

    for (int x_i = x_s; x_i <= x_e; x_i ++) {
      for (int y_i = y_s; y_i <= y_e; y_i ++) {
        int id = x_i*num + y_i;
        std::set<particle_t*>::iterator it2;
        it2 = gridContainer[id].pArray.begin();
        for (; it2 != gridContainer[id].pArray.end(); ++ it2) {
          if(it == it2) continue; // continue if same particle
          apply_force( **it, **it2, dmin, davg, navg);
        }
      }
    }
  }
}

void move_particles(particle_t *particles, int n, int i){
  moveWithGrid(particles[i], gridContainer);
}

void dele_grid() {
  delete[] gridContainer;
}


void debug_printSet(std::set<particle_t*> s) {
  std::set<particle_t*>::iterator it;
  for (it = s.begin(); it != s.end(); ++it) {
    printf("the location of the particle are: x: %f; y: %f\n", (**it).x, (**it).y);
  }
}
