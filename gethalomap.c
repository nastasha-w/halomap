#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifdef _OPENMP
#include <omp.h>
#endif

// dimension of position arrays (fast index)
#define DIMENSION 3 

/*
Functions to return a 2d boolean array of which parts of a 2d map
belong to a set of cilindrical/spherical halos with positions
specified
assumes square pixels
includes periodic boundaries
assumes pixel 0,0 is at position 0, 0
*/


/* 
--------------
| prototypes |
--------------
*/
void gethalomap_2d_incl(int xpix, int ypix, float pixsize,
                        float *positions, float* radii, int lenpos,
                        int Axis1, int Axis2, int Axis3,
                        int periodic, float size,
                        int* outmap);
void gethalomap_2d_inexcl(int xpix, int ypix, float pixsize,
                         float *inpositions, float *inradii, int inlenpos,
                         float *expositions, float *exradii, int exlenpos, int normdiff,
                         int Axis1, int Axis2, int Axis3,
                         int periodic, float size,
                         int *outmap);
float diff_periodic(float x1, float x2, float size);

/*
--------
| main |
--------
Contains only test cases
*/



int main(void) 
{
    printf("Hello world");
    return 0;
}





/*
--------------------------------
| the functions doing the work |
--------------------------------
*/

float diff_periodic(float x1, float x2, float size)
{
    float halfsize = 0.5*size;
    if(x1-x2 <= halfsize)
    {   
        if(x2-x1 <= halfsize) 
        {
            return x1-x2;	
        }
        else
        {
            return x1-x2 + size;       
        }
    } 
    else
    {
        return x1-x2 - size;    
    }
}


/* 2d information only, no excluded overlap regions 
   assumes all radii are < box size
   xpix: 
       number of pixels in the final map, along Axis1
   ypix: 
       number of pixels in the final map, along Axis2
   pixsize:
       size of each pixel [position units], same in both dimensions
   positions: 
       positions of the haloes [centers, coordinates x, y, z are the fast 
       index]
   radii: 
       radii of the haloes to include [position units]
   lenpos:    
       number of haloes in the positions/radii arrays
   Axis1: 
       X axis in the final map; corresponds the positions x, y, or z 
       (0, 1, 2 respectively) axis 
   Axis2: 
       Y axis in the final map; corresponds the positions x, y, or z 
       (0, 1, 2 respectively) axis 
   Axis3: 
       Z axis (sightline direction) in the final map; corresponds the 
       positions x, y, or z (0, 1, 2 respectively) axis
   periodic:
       are the edges of the output map perdiodic (1) or not (0)
   size:
       size of the periodic volume [position units]
   
   outmap:
       array that will hold the map of 0/1 (False/True) values, indicating
       whether each pixel is within an impact parameter <radius> of one of
       the haloes
*/
void gethalomap_2d_incl(int xpix, int ypix, float pixsize,
                        float *positions, float* radii, int lenpos,
                        int Axis1, int Axis2, int Axis3,
                        int periodic, float size,
                        int* outmap)
         
{   // outmap -> bool-like array: 1 = in region, 0 = outside region; y is the fast index   
    printf("---- gethalomap_2d_incl ----\n");
    #ifdef _OPENMP
    #pragma omp parallel 
    {
        int size = omp_get_num_threads();
        int rank = omp_get_thread_num();
        if(rank==0){printf("called with OpenMP on %i cores\n",size);}
    }
    #else
    printf("serial version\n");
    #endif 
    printf("input parameters:\n");
    printf("xpix:\t%i\typix:\t%i\n", xpix, ypix);
    printf("pixsize:\t%f\n", pixsize);
    printf("size:\t%f\n", size);
    printf("periodic:\t%i\n", periodic);
    printf("Axis1:\t%i\tAxis2:\t%i\tAxis3:\t%i\n", Axis1, Axis2, Axis3);
    printf("lenpos:\t%i\n", lenpos);
    //printf("test: %i mod %i = %i\n", -2, 3, -2 % 3);
    //printf("test: %i mod %i = %i\n", 1, 3, 1 % 3);
    //printf("test: %i mod %i = %i\n", 5, 3, 5 % 3);
    
    // add omp parallel loop option; different threads would set to same value -> no need to
    // check race conditions
    long int i;
    for(i=0; i<xpix*ypix; i++)
    {
        outmap[i] = 0;
    }
    #ifdef _OPENMP  
    #pragma omp parallel for
    #endif
    for(i=0; i<lenpos; i++)
    {
        // pixel position = (pixel index + 0.5) * pixsize
        float xpos = positions[i*DIMENSION + Axis1];
        float ypos = positions[i*DIMENSION + Axis2];
        float radius  = radii[i];
        //printf("selecting x element: %li\n", i*DIMENSION + Axis1);
        //printf("selecting y element: %li\n", i*DIMENSION + Axis2);
        //printf("ignoring  z element: %li\n", i*DIMENSION + Axis3);
        //printf("Using x, y, r = %f, %f, %f\n", xpos, ypos, radius);
        float r2 = radius*radius;
        long int xpmin = floor((xpos - radius) / pixsize - 0.5);
        long int xpmax =  ceil((xpos + radius) / pixsize - 0.5);
        long int ypmin = floor((ypos - radius) / pixsize - 0.5);
        long int ypmax =  ceil((ypos + radius) / pixsize - 0.5);
        /* printf("xpmin, xpmax, ypmin, ypmax: %li, %li, %li, %li\n", xpmin, xpmax, ypmin, ypmax); */
       
        if(periodic == 0)  // deal with edge overlaps of halo i; periodic indices are taken modulo size later
        {
            if(xpmax >= xpix){xpmax = xpix-1;}
            if(ypmax >= ypix){ypmax = ypix-1;}
            if(xpmin < 0){xpmin = 0;}
            if(ypmin < 0){ypmin = 0;}
        }
        long int _xi, _yi;
        if(periodic)
        {
            for(_xi=xpmin; _xi<xpmax+1; _xi++)
                for(_yi=ypmin; _yi<ypmax+1; _yi++)
                {   
                    long int xi, yi;
                    xi = (_xi + xpix) % xpix; // modulo gives remainders in C; this only works if halos are always < box, but that should not be an issue
                    yi = (_yi + ypix) % ypix;             
                
                    float xoff = diff_periodic((xi + 0.5) * pixsize, xpos, size);
                    float yoff = diff_periodic((yi + 0.5) * pixsize, ypos, size); 
                    if(xoff*xoff + yoff*yoff <= r2)
                    {
                        outmap[xi*ypix +yi] = 1;
                        //printf("Setting %li, %li to True\n", xi, yi);
                    }
                }
         } // end of if periodic
         else
         {
         for(_xi=xpmin; _xi<xpmax+1; _xi++)
                for(_yi=ypmin; _yi<ypmax+1; _yi++)
                {                
                    float xoff = diff_periodic((_xi + 0.5) * pixsize, xpos, size);
                    float yoff = diff_periodic((_yi + 0.5) * pixsize, ypos, size); 
                    if(xoff*xoff + yoff*yoff <= r2)
                    {
                        outmap[_xi * ypix + _yi] = 1;
                        //printf("Setting %li, %li to True\n", xi, yi);
                    }
                }
         
         }
        
    }
    printf("---- gethalomap_2d_incl finished----\n");
}

/* 2d information only, excluded overlap regions. Note that excluded and included should not overlap;
   if they do, the program will try to include and exlcude the same halo, and behaviour will be 
   undefined 
   in case of a `tie break', the normed/unnormed comparison is used in the sunnormed/normed case 
   
   inputs match the _incl version, except that:
   - the halo data include in- and ex- versions, for the haloes to include and
     exclude
   normdiff:
       when a halo to include and a halo to exlcude overlap, the overlapping
       pixels can be attributed to a halo based on one of two criteria.
       for normdiff 0 (False), the pixel is attributed to whichever halo has
       its center closest to the pixel
       for normdiff 1 (True), the pixel is attributed to whichever halo has
       a smaller <impact parameter to halo centre> / halo radius
       If the distances are equal by the chosen metric, the other metric is
       used as a tiebreaker.   
 */
void gethalomap_2d_inexcl(int xpix, int ypix, float pixsize,
                         float *inpositions, float *inradii, int inlenpos,
                         float *expositions, float *exradii, int exlenpos, int normdiff,
                         int Axis1, int Axis2, int Axis3,
                         int periodic, float size,
                         int *outmap)
                        
{
    printf("---- gethalomap_2d_inexcl ----\n");
    // bool-like array: 1 = in region, 0 = outside region; y is the fast index
    long int i;
    for(i=0; i<xpix*ypix; i++)
    {
    	  outmap[i] = 0; // default: not in a halo
    }
    // assumes excluded position list is small enough for int indices;
    // buffer contains only the halos that overlap with halo i; 
    // index zero is for the halo that overlapped with the last pixel 
    #ifdef _OPENMP
    #pragma omp parallel 
    {
    int numt = omp_get_num_threads();
    int rank = omp_get_thread_num();
    if(rank==0){printf("called with OpenMP on %i cores\n", numt);}
    #else
    printf("serial version\n");
    int rank = 0;
    #endif 
    if(rank==0)
    {    
        printf("input parameters:\n");
        printf("xpix:\t%i\typix:\t%i\n", xpix, ypix);
        printf("pixsize:\t%f\n", pixsize);
        printf("size:\t%f\n", size);
        printf("periodic:\t%i\n", periodic);
        printf("Axis1:\t%i\tAxis2:\t%i\tAxis3:\t%i\n", Axis1, Axis2, Axis3);
        printf("inlenpos:\t%i\n", inlenpos);
        printf("exlenpos:\t%i\n", exlenpos);
        printf("in/exclude pixels based on distance/radius:\t%i\n", normdiff);
    }
    int *exbuff = malloc((exlenpos + 1)*sizeof(int)); 
    // add omp parallel loop option; different threads would set to same value -> no need to
    // check race conditions. Do need to use separate overlapping halo buffers in that case.
    #ifdef _OPENMP
    #pragma omp for
    #endif
    for(i=0; i<inlenpos; i++)
    {
    	  //printf("thread %i handling halo %li\n", rank, i);
        // pixel position = (pixel index + 0.5) * pixsize
        float xpos = inpositions[i*DIMENSION + Axis1];
        float ypos = inpositions[i*DIMENSION + Axis2];
        float radius  = inradii[i];
        // make list of exclude halos overlapping with halo i
        int j;
        for(j=0; j < exlenpos + 1; j++)
        {
            exbuff[j] = -1; 	  
        }
        int lenbuff = 0;
        for(j=0; j<exlenpos; j++)
        {
            float xdiff = diff_periodic(xpos, expositions[j*DIMENSION + Axis1], size);
            float ydiff = diff_periodic(ypos, expositions[j*DIMENSION + Axis2], size);
            float exradius = exradii[j];
            if(xdiff * xdiff + ydiff * ydiff <= (radius + exradius) * (radius + exradius))
            {
                lenbuff += 1;
                exbuff[lenbuff] = j;
            } 	  
        }
        //printf("Found %i halos overlapping halo %li\n", lenbuff, i);
        long int xpmin = floor((xpos - radius) / pixsize - 0.5);
        long int xpmax =  ceil((xpos + radius) / pixsize - 0.5);
        long int ypmin = floor((ypos - radius) / pixsize - 0.5);
        long int ypmax =  ceil((ypos + radius) / pixsize - 0.5);

        if(periodic == 0) // cap edge overlaps for halos
        {
            if(xpmax >= xpix){xpmax = xpix-1;}
            if(ypmax >= ypix){ypmax = ypix-1;}
            if(xpmin < 0){xpmin = 0;}
            if(ypmin < 0){ypmin = 0;}
        }
        long int _xi, _yi;
        int excl; 
        float r2 = radius*radius;
        for(_xi=xpmin; _xi<xpmax+1; _xi++)
            for(_yi=ypmin; _yi<ypmax+1; _yi++)
            {   
                long int xi, yi;
                if(periodic)
                {
                    xi = (_xi + xpix) % xpix;
                    yi = (_yi + ypix) % ypix;              
                }
                else
                {
                    xi = _xi;
                    yi = _yi;
                }
                if(outmap[xi*ypix + yi] == 1) // already set for a different include halo
                {
                	  continue;
                }
                float xoff = diff_periodic((xi + 0.5)*pixsize, xpos, size);
                float yoff = diff_periodic((yi + 0.5)*pixsize, ypos, size);  
                float rpos2 = xoff*xoff + yoff*yoff;
                if(rpos2 <= r2)
                {           	
                    excl = 0;  
                    float rexpos2, nrpos2, nrexpos2, exxoff, exyoff;
                    for(j=0; j<lenbuff+1; j++)
                    {   
                        if(exbuff[j] == -1)
                        {
                        	  continue; // first index of buffer is undefined 
                        }
                        exxoff = diff_periodic((xi + 0.5) * pixsize, expositions[exbuff[j] * DIMENSION + Axis1], size);
                        exyoff = diff_periodic((yi + 0.5) * pixsize, expositions[exbuff[j] * DIMENSION + Axis2], size);
                        rexpos2 = exxoff * exxoff + exyoff * exyoff;
                        nrexpos2 = rexpos2 / (exradii[exbuff[j]] * exradii[exbuff[j]]);
                        nrpos2 = rpos2/r2;                        	 
                        if(normdiff && (nrpos2 > nrexpos2 || (nrexpos2 == nrpos2 && rpos2 > rexpos2) )) // we've found a different halo the pixel should belong to
                        {
                        	 excl = 1;
                        	 exbuff[0] = exbuff[j];
                        	 break;
                        }
                        else if(!normdiff && (rpos2 > rexpos2 || (rexpos2 == rpos2 && nrpos2 > nrexpos2)))
                        {
                        	 excl = 1;
                        	 exbuff[0] = exbuff[j];
                        	 break;
                        }
                    }
                    if(excl == 0)
                    {                 	
                        //printf("Including pixel %li, %li\n", xi, yi);
                        outmap[xi*ypix +yi] = 1;
                    }                
                } // end of if rpos < r2
            } // end of halo x/y pixel loops
    } // end of include halo loop
    free(exbuff);
    // check race conditions. Do need to use separate overlapping halo buffers in that case.
    #ifdef _OPENMP
    } // end of OMP parallel
    #endif
     printf("---- gethalomap_2d_inexcl finished----\n");
}


/*
N nearest neighbour distance to a pixel (given a set of galaxies)
NVM, it looks like this is a whole thing, with dedicated python libraries
determine largest nn distance between galaxies (sets search radius for pixels):
double maxnndg = 0.;
loop 1 galaxies:
   sort galaxies - galaxy 1 distance (qsort?)
   nnd_this = get element no. N (element 0 will be zero: galaxy distance to itself)
   if nnd_this > maxnng:
      maxnndg = nnd_this

initialize output map
initialize size-N array
loop over pixel blocks:
     get all galaxies within margin * maxnndg of the edges
     check if the range is large enough: for all the corners, there
       should be at least N galaxies within the margin region (otherwise,
       the nearest N galaxies in the block sample might not be the nearest
       N overall). Actually, would need to check the edges...
     for each pixel in the block:
         sort selected galaxies by distance to pixel
         get distance N (element N-1)
         put into output array            
         (use in-place sort and re-use array; should be pretty well sorted then, mostly)
          
*/
 
 

