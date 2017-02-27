/************************************************************
 * partition.c  main program for testing out partition fn   *
 *                                                          *
 * Created by William DeMeo                                 *
 * on 2011.09.24                                            *
 *                                                          *
 ************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <prototypes.h>

int partition(double *v, double start, double end, double step);

void main()
{
  double *sigLevels, *ptr;
  int i, j, len=0, place=1;
  //sigLevels = [0.001 (0.005:0.005:0.10) (0.125:0.025:0.20) (0.80:0.025:0.875) (0.90:0.005:0.995) 0.999];
  double sigs[4][3] = {{0.005,0.005,0.10}, {0.125,0.025,0.20}, {0.80,0.025,0.875}, {0.90,0.005,0.995}};

  sigLevels = dmalloc(50);
  sigLevels[0] = 0.001;
  ptr = sigLevels+1;
  for (i=0; i<4; ++i){
    printf("\nsigs[%d] = %f, %f, %f", i, sigs[i][0],sigs[i][2],sigs[i][1]);
    len = partition(ptr,sigs[i][0],sigs[i][2],sigs[i][1]);
    assert(len > 1);
    for (j=0;j<len;++j){
      printf("\nsigLevels[%d] : %f", place+j, sigLevels[place+j]);
    }
    ptr += len;
    place = place+len;
  }
  sigLevels[49]=0.999;
}

int partition(double *v, double start, double end, double step){
  int len, i;

  printf("\n(end - start)/step + 1 = %f",((end - start)/step)+1.0);
  len= (int)(((end - start)/step)+1.0);
  printf("\nlen = %d",len);

  for (i=0;i<len;++i) {
    v[i] = start+(step*i);
    printf("\nv[%d] = %f",i,v[i]);
  }
  assert(v[len-1]==end);
  return len;
}
    
  
  
