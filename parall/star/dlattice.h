#ifndef DLATTICE_H
#define DLATTICE_H

#include <stdlib.h>
#include <emmintrin.h>
#include <xmmintrin.h>
#include <unistd.h>

#include "optimization.h"
#include "local_types.h"

#define dir_number 16
static int dx_neighbor[dir_number]  __attribute__ ( (aligned(16)) ) = { 1, 1,-1,-1, 0, 2, 2, 0, 2, 2, 0,-2,-2, 0,-2,-2}; // x-координаты соседей
static int dy_neighbor[dir_number]  __attribute__ ( (aligned(16)) )= { 1,-1, 1,-1, 2, 0, 2,-2, 0,-2, 2, 0, 2,-2, 0,-2}; // y-координаты соседей
static int dz_neighbor[dir_number]  __attribute__ ( (aligned(16)) )= { 1,-1,-1, 1, 2, 2, 0,-2,-2, 0,-2,-2, 0, 2, 2, 0}; // z-координаты соседей
static int pr_neighbor[dir_number]  __attribute__ ( (aligned(16)) )= {-1,-1,-1,-1, 0, 0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3}; // номер предыдущего

typedef struct {
  int x[dir_number] __attribute__ ( (aligned(16)) );
  int y[dir_number] __attribute__ ( (aligned(16)) );
  int z[dir_number] __attribute__ ( (aligned(16)) );
} neighbors_t;

#define LPTR(type,LX,LY,newname,oldname) type (*newname)[LY][LX] \
      __attribute__ ( (aligned(16)) ) = (type (*)[LY][LX])oldname

template<typename Elem> 
class dlattice_t{

public:

  Elem *lat;
  ulong Lx, Ly, Lz;

  __m128i v_0, v_1, v_2, v_3;

public:

  inline ulong lsize()
  {
    return Lx * Ly * Lz * sizeof(Elem);
  }

  dlattice_t()
  {
    Lx = 0;
    Ly = 0;
    Lz = 0;
    lat = NULL;

    v_0 = _mm_set1_epi32(0);
    v_1 = _mm_set1_epi32(1);
    v_2 = _mm_set1_epi32(2);
    v_3 = _mm_set1_epi32(3);
  }

  void reset()
  {
    if( lat )
        memset(lat, 0, lsize());
  }

  void reset(Elem val)
  {
    int x, y, z;
    if( lat == NULL )
      return;
    LPTR(Elem, Lx, Ly, lattice, lat);
    for(z = 0; z < Lz/4; z++){
      for(y = 0; y < Ly/2; y++){
        for(x = 0; x < Lx; x++){
          lattice[z][y][x] = val;
        }
      }
    }
  }

  
  int tune(int _Lx, int _Ly, int _Lz)
  {
     if((_Lx % 4) || (_Ly % 4) || (_Lz % 4)) {
        printf("(Lx %% 4 != 0) OR (Ly %% 4 != 0) OR (Lz %% 4 != 0): (%ld,%ld,%ld)\n", Lx, Ly, Lz);
        abort();
     }
                                                    
     Lx = _Lx;
     Ly = _Ly;
     Lz = _Lz;

     if(lat)
       free(lat);

     if(posix_memalign((void**) &lat, 16, lsize())) {
       printf("Cannot allocate %lu bytes\n", lsize());
       abort();
     }
     reset();
     return 0;
  }                                                                                                                                                          

  static inline void convert(int x, int y, int z, int &x1, int &y1, int &z1)
  {
    x1 = x;
    /*
    if( unlikely(y < 0) )
      y--;
    y1 = (y/4)*2 + (y%4)/2;
    */
    y1 = ( (y >> 2) << 1) + ((y & 0x3) >> 1);
    
    /*
    if( unlikely(z < 0) )
      z -= 3;
    z1 = z/4;
    */
    
    z1 = z >> 2;
  }

  static inline void r_convert(int x1, int y1, int z1, int &x, int &y, int &z)
  {
    // coordinates of x, y, z inside elementary cube
    int ix, iy, iz;
    // shifts of coordinates to put elementary cube to (0,0,0)
    int sx=x1, sy = y1;
   
    // if cube has negative coordinates it has to be shifted one 
    // cube more sinse -3/4 = 0 and 3/4 = 0
    // cube with negative X=-3 has number -1
    // cube with positive X=3 has number 0
    if( sx < 0 )
      sx -= 3;
    sx = (sx/4)*4;
    
    if( sy < 0 )
      sy -= 1;
    sy = (sy/2)*2;
    
    // Count incube x coordinate
    ix = x1 - sx;
    x = ix;
    // Count incube y coordinate
    // y is reduced by 2 (y/2) for packing purposes so ve loose 
    // y%2 component. But knowing that x and y have equal parity
    // we can use x%2 to find out y%2
    //  ecube     y   y/2  x%2  [(y/2)*2 + x%2]
    // (0,0,0)    0    0    0          0
    // (0,2,2)    2    1    0          2
    // (1,1,1)    1    0    1          1
    // (1,3,3)    3    1    1          3
    // (2,0,2)    0    0    0          0
    // (2,2,0)    2    1    0          2
    // (3,1,3)    1    0    1          1
    // (3,3,1)    3    1    1          3
    iy = y1 - sy;
    y = 2*(iy%2) + ix%2;
    
    // Diamond lattice elementary cube has strong
    // dep's between (x,y) coordinates and z coord
    // on odd-z level coordinates are also odd
    // on even level - they are even
    // Accordint to the table below for z=0,1 (x-y) is zero
    // for z=2,3 (x-y) = +-2 (for module do (x-y)(x-y)/2
    // next to consider [+1] for z=1,3 we can use parity of one
    // if coordinates, for example x%2
    //  ecube     z   x-y  [ (x-y)(x-y)/2 + x%2 ]
    // (0,0,0)    0    0            0
    // (0,2,2)    2   -2            2
    // (1,1,1)    1    0            1
    // (1,3,3)    3   -2            3
    // (2,0,2)    2    2            2
    // (2,2,0)    0    0            0
    // (3,1,3)    3    2            3
    // (3,3,1)    1    0            1
    
    z = (y-x)*(y-x)/2 + x%2;
    
    // shift all coordinates to appropriate positions
    x += sx;
    y += sy*2;
    z += z1*4;
  }

  inline Elem &operator () (uint x, uint y, uint z) {
//  inline Elem &o(uint x, uint y, uint z) {

    // indexes of corresponding elementary 
    // cube in the grid
    int lx, ly, lz;
    LPTR(Elem, Lx, Ly, lattice, lat);
    
    convert(x, y, z, lx, ly, lz);

    return lattice[lz][ly][lx];
  }

  inline Elem &l(int x, int y, int z) {
    // access to elements using packed coordinates
    LPTR(Elem, Lx, Ly, lattice, lat);
    if( unlikely(y>Lx/2 || z>Lz/4) ){
      printf("Use original coordinates instead of packed!\n");
      abort();
    }
    return lattice[z][y][x];
  }

  inline  Elem &l(int c[3]) {
    // access to elements using packed coordinates
    LPTR(Elem, Lx, Ly, lattice, lat);
    return lattice[c[2]][c[1]][c[0]];
  }




  inline void boundary(int &x, int &y, int &z) { // изменить (x,y,z), если нужно, так чтобы попасть внутрь ящика
    if ( unlikely(x<0) )   x+=Lx;
    if ( unlikely(x>=Lx) ) x-=Lx;
    if ( unlikely(y<0) )   y+=Ly;
    if ( unlikely(y>=Ly) ) y-=Ly;
    if ( unlikely((z<0) || (z>=Lz)) ) {
      fprintf(stderr,"Error!!! z out of range (%d,%d,%d)\n",x,y,z);
      abort();
    }
  }
  
  inline void boundary_v(neighbors_t &nbs) // изменить (x,y,z), если нужно, так чтобы попасть внутрь ящика
  {
    static __m128i low_v = _mm_set1_epi32(0);
    static __m128i upx_v = _mm_set_epi32(Lx,Lx,Lx,Lx);
    static __m128i upy_v = _mm_set_epi32(Ly,Ly,Ly,Ly);
    static __m128i upz_v = _mm_set_epi32(Lz,Lz,Lz,Lz);
    __m128i cklow, ckup;
    __m128i zck[4];
    static int count = 0;
    
    __m128i *xptr = (__m128i *)nbs.x;
    __m128i *yptr = (__m128i *)nbs.y;
    __m128i *zptr = (__m128i *)nbs.z;
    
    for(int i=0;i<4;i++){
      // check bounds of X
      cklow = _mm_cmplt_epi32(  *( xptr + i), low_v);
      *(xptr + i) = _mm_add_epi32( *( xptr + i), _mm_and_si128( cklow, upx_v) ) ;
      
      ckup = _mm_cmplt_epi32( upx_v, *( xptr + i));
      ckup = _mm_add_epi32( ckup, _mm_cmpeq_epi32( upx_v, *( xptr + i)) );
      *( xptr + i) = _mm_sub_epi32( *(xptr + i), _mm_and_si128( ckup, upx_v) );

      // check bounds of Y
      cklow = _mm_cmplt_epi32(  *(yptr + i), low_v);
      *( yptr + i) = _mm_add_epi32( *( yptr + i), _mm_and_si128( cklow, upy_v) ) ;
      
      ckup = _mm_cmplt_epi32( upy_v, *(yptr + i) );
      ckup = _mm_add_epi32( ckup, _mm_cmpeq_epi32( upy_v, *( yptr + i)) );
      *(yptr + i) = _mm_sub_epi32( *(yptr + i), _mm_and_si128( ckup, upy_v) );
    }
  }  

  inline int neighbors(int x, int y, int z, neighbors_t &nbs)
  {
  /*
    // DEBUG
    neighbors_t nbs1;
    memset(&nbs,0,sizeof(nbs) );
    memset(&nbs1,0,sizeof(nbs1) );
    //DEBUG 
    int factor = neighbors_s(x, y, z, nbs);

    //DEBUG
    neighbors_s1(x, y, z, nbs1);
    if( memcmp(&nbs, &nbs1, sizeof(nbs1) ) )
      abort(); 
    //DEBUG
   */
    return neighbors_s1(x, y, z, nbs);
    
  }
  
  int neighbors_s(int x, int y, int z, neighbors_t &nbs)
  {
    int factor, dir;
    if ( z%2==0 )
      factor=1;
    else 
      factor=-1;
  
    for (dir=0; dir<dir_number; dir++) {
      nbs.x[dir] = x+factor*dx_neighbor[dir];
      nbs.y[dir] = y+factor*dy_neighbor[dir];
      nbs.z[dir] = z+factor*dz_neighbor[dir];
      boundary(nbs.x[dir], nbs.y[dir], nbs.z[dir]);
    }
    return factor;
  }

  int neighbors_s1(int x, int y, int z, neighbors_t &nbs)
  {
    int factor, dir;
    if ( z%2==0 )
      factor=1;
    else 
      factor=-1;
  
    for (dir=0; dir<dir_number; dir++) {
      nbs.x[dir] = x+factor*dx_neighbor[dir];
      nbs.y[dir] = y+factor*dy_neighbor[dir];
      nbs.z[dir] = z+factor*dz_neighbor[dir];
    }
    
    boundary_v(nbs);
    return factor;
  }

  void one_neighbor(int x, int y, int z, int dir, int &x2, int &y2, int &z2)
  {
    int factor; 
    if (z%2==0)
      factor=1;
    else
      factor=-1;

    x2 = x + factor*dx_neighbor[dir];
    y2 = y + factor*dy_neighbor[dir];
    z2 = z + factor*dz_neighbor[dir];
    boundary(x2,y2,z2);
  }


};

#endif
