#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>
#include <string.h>
#include <math.h>
#include <set>
#include "common.h"


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

	// here just check and print the location of the subgrid and the particals inside
	/**
	for(i = 0; i < num * num; i ++) {
		printf("the corner is x: %f; y: %f\n", gridContainer[i].x, gridContainer[i].y);
		// then print the particals inside
		debug_printSet(gridContainer[i].pArray);
	}**/
}

// update the grid

void update_grid(double *dmin, double *davg, int *navg, particle_t *particles, int n){
	// here we just need to check the 9 subgrids
	int i = 0;
	for (i = 0; i < num * num; i ++) {
		int k = 0; // k is the iteration of the particals in the subgrid
		std::set<particle_t*>::iterator it;
		it = gridContainer[i].pArray.begin();


		for (; it != gridContainer[i].pArray.end(); ++ it) {

			// we need to assign the a = 0
			(**it).ax = 0.0; 
			(**it).ay = 0.0;

			// then iterate though the 9 surrending subgrids

			int x_s = (i/num > 0) ? (i/num - 1) : 0;
			int x_e = (i/num < num - 1) ? (i/num + 1) : num - 1;
			int y_s = (i%num > 0) ? (i%num - 1) : 0;
			int y_e = (i%num < num - 1) ? (i%num + 1) : num - 1;
			// here x_s and x_e are the start and end index for the horizental axis
			// .....y ..... y....                            ........vertical ....       

			//printf("the surrending is x_s: %d; x_e: %d; y_s: %d; y_e %d; i = %d\n",x_s,x_e,y_s,y_e,i );
			int x_i = x_s;
			for (; x_i <= x_e; x_i ++) {
				int y_i = y_s;
				for (; y_i <= y_e; y_i ++) {
					// the subgrid index = x_i*num + y_i
					int k2 = 0; 
					int id = x_i*num + y_i; // the id is the index for the subgrid
					//printf("id is %d\n", id);
					std::set<particle_t*>::iterator it2;
					it2 = gridContainer[id].pArray.begin();

					for (; it2 != gridContainer[id].pArray.end(); ++ it2) {
						// go though all the particals in the surrending grids
						//particle_t* parti2 = gridContainer[id].pArray[k2];
						if(*it == *it2) continue; // if two particles are the same one, continue
						apply_force( **it, **it2, dmin, davg, navg);
					}
				}
			}
			// then here we move this iterm
		}
	}
	// here we have collected all the force
	for (i = 0; i < n; i ++) {
		moveWithGrid(particles[i], gridContainer);
	}

	// below is just for debug
	/**
	printf("after the movement\n");
	for(i = 0; i < num * num; i ++) {

		printf("the corner is x: %f; y: %f\n", gridContainer[i].x, gridContainer[i].y);
		// then print the particals inside
		debug_printSet(gridContainer[i].pArray);
	}
	**/

}

void dele_grid() {
	
	delete[] gridContainer;
}


void debug_printSet(std::set<particle_t*> s) {
	// here print the location of all the particles in the set
	std::set<particle_t*>::iterator it;
	for (it = s.begin(); it != s.end(); ++it) {
		printf("the location of the particle are: x: %f; y: %f\n", (**it).x, (**it).y);
	} 
}