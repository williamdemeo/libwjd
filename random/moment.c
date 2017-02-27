
/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   moment.c
  
   Three routines for for computing the mean and variance
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
/* dmoment: The desk calculator algorithm
   arguments:
              n = length of data[]
              data = a nx1 array of doubles
              ave =(on entry)= the address of *ave (call by reference)
              ave =(on exit)= average of data[]
              var =(on entry)= the address of *var (call by reference)
              var =(on exit)= variance of data[]
              */
void dmoment(long n, double * data, double *ave, double *var)
{
     long i;
  
     *ave=0; *var=0;
  
     for(i=0;i<n;i++){
          *ave += data[i];
          *var += data[i]*data[i];
     }
     *ave /= (double) n;
     *var = (*var - (double)n*((*ave)*(*ave)))/(double)(n-1);
}

/* pmoment: The provisional means algorithm
   arguments:
              n = length of data[]
              data = a nx1 array of doubles
              ave =(on exit)= the average of data[]
              var =(on exit)= the variance of data[]
              */
void pmoment(long n, double * data, double *ave, double *var)
{
     long i;

     *ave = data[0];
     *var=0;
  
     for(i=1;i<n;i++){
          *var +=((double)i/(double)(i+1))*(data[i]-*ave)*(data[i]-*ave);
          *ave *= ((double)i/(double)(i+1));
          *ave += data[i]/(double)(i+1);
     }
     *var /= (double)(n-1);
}

/* cmoment: centering around the first observation 
   arguments:
              n  length of data vector
              data = a nx1 array of doubles
              ave =(on exit)= the average of data[]
              var =(on exit)= the variance of data[]
              */
void cmoment(long n, const double *data, double *ave, double *var)
{
     long i;
  
     *ave=0;
     *var=0;
     for(i=1;i<n;i++){
          *ave += data[i] - data[0];
          *var += (data[i] - data[0])*(data[i] - data[0]);
     }
     *ave /= (double)n;
     *var = (*var - (double)n * ((*ave)*(*ave)))/(double)(n-1);
     *ave += data[0];
}

/* moment: centering around the first observation 
           (same as cmomment, except n is passed by reference)
   arguments:
              n  length of data vector
              data = a nx1 array of doubles
              ave =(on exit)= the average of data[]
              var =(on exit)= the variance of data[]
              */
void moment(const long *n, const double *data, double *ave, double *var)
{
     long i;
  
     *ave=0;
     *var=0;
     for(i=1;i<*n;i++){
          *ave += data[i] - data[0];
          *var += (data[i] - data[0])*(data[i] - data[0]);
     }
     *ave /= (double)(*n);
     *var = (*var - (double)(*n) * ((*ave)*(*ave)))/(double)(*n-1);
     *ave += data[0];
}


