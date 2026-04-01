#ifndef DLATTICE_H
#define DLATTICE_H

#include <stdlib.h>
#include <emmintrin.h>

#include "optimization.h"
#include "local_types.h"

#define Lx_max 200
#define Ly_max 200
#define Lz_max 200


static uchar fneighbors[2][4][3] = 
{
  { {1,0,0}, {-1,0,0}, {0,1,0},  {0,0,1} },
  { {1,0,0}, {-1,0,0}, {0,-1,0}, {0,0,-1} }
};

// ----------  Elementary cube of diamond lattice ----------------------//

typedef uchar atom_ecube_t;

typedef struct force_t{
  double x, y, z;
} force_ecube_t[8];


// vectors for mapping of 
// packed cube to real coordinates
//  z=0   z=1   z=2   z=3
//  0000  0001  0000  0100
//  0010  0000  1000  0000
//  0000  0100  0000  0001
//  1000  0000  0010  0000
//# (1,2) (3,4) (5,6) (7,8)

static uchar cube_elem_map[4][4][4] = 
{
  // z = 0
  { {1,0,0,0}, {0,0,0,0}, {0,0,2,0}, {0,0,0,0} },
  // z = 1
  { {0,0,0,0}, {0,3,0,0}, {0,0,0,0}, {0,0,0,4} },
  // z = 2
  { {0,0,5,0}, {0,0,0,0}, {6,0,0,0}, {0,7,0,0} },
  // z = 3
  { {0,0,0,0}, {0,0,0,7}, {0,0,0,0}, {0,8,0,0} }
};

static unsigned char elem_cube_map[8][2] = 
{
  {0, 0}, // point #1
  {2, 2}, // point #2
  {1, 1}, // point #3
  {3, 3}, // point #4
  {0, 2}, // point #5
  {2, 0}, // point #6
  {4, 1}, // point #7
  {1, 4} // point #8
};


template<int Lx, int Ly, int Lz, typename Elem> 
class dlattice_t{

public:

  struct ecube_coords
  {
    ulong cx, cy, cz;
    uchar in;
  };

private:

  Elem lattice[Lz/4][Ly/4][Lx/4][8];
  ulong Cx, Cy, Cz;
  __m128i mod_mask;

  /*
   * Convert sequental atom number into diamond
   * lattice local coordinates: 
   * (elementary cube x, y ,z; internal point number [0..7])
   */
  int globseq_2_ecube(ulong n, ecube_coords &c)
  {
    ulong n1, n2;
    ulong z1, y1;
    
    n1 = n % (Cx*Cy*8);
    n2 = n1 % (Cx*Cy*2);
    z1 = n1/(Cx*Cy*2);
    y1 = n2/Cx;
    c.cz = n/(Cx*Cy*8);
    c.cy = y1/2;
    c.cx = n2 % Cx;
    c.in = z1*2 + (y1 % 2);
  }

public:
  
  dlattice_t()
  {
    Cx = Lx/4;
    Cy = Lx/4;
    Cz = Lz/4;
    mod_mask = _mm_set1_epi32(0x3);
  }
    
  inline Elem &operator () (uint x, uint y, uint z) {

    // indexes of corresponding elementary 
    // cube in the grid
    uint lx = x/4, ly = y/4, lz = z/4; 
    // indexes of the point inside 
    // elementary cube
    uint idx1 = cube_elem_map[z%4][y%4][x%4] - 1;
    
    // In debug case
    //if( lx > Lx || ly > Ly || lz > Lz || idx < 0 ){
    //  printf("%s: bad indexes!\n", __FUNCTION__);
    //  abort(); 
    //}
    //printf("[(%d,%d,%d),%d] ",lx,ly,lz,idx);

    return lattice[lz][ly][lx][idx];
/*

    lcoord c = { {x, y, z} }, cube, inc;
    cube.vect = _mm_srli_epi32(c.vect, 2);
    inc.vect = _mm_and_si128(c.vect, mod_mask);
    uint idx = cube_elem_map[inc.coord.z][inc.coord.y][inc.coord.x] - 1;


    printf("(%d, %d, %d): (%d,%d,%d), (%d,%d,%d), %d\n",x,y,z, 
      cube.coord.x, cube.coord.y, cube.coord.z,
      inc.coord.x, inc.coord.y, inc.coord.z,
      idx);
    fflush(stdout);
    
    if( lx != cube.coord.x || ly != cube.coord.y || lz !=cube.coord.z || idx1 != idx ){
      abort();
    }
    
    return lattice[cube.coord.z][cube.coord.y][cube.coord.x][idx];
*/
  }
  
  // Fill cubes one by one
  void test_fill()
  {
    int counter = 1;
    int x, y, z, i;
    for(z=0;z<Cz;z++)
      for(y=0; y<Cy; y++)
        for(x=0; x<Cx; x++){
          for(i=0;i<8;i++)
            lattice[z][y][x][i] = counter++;
        }
  }

};

#endif