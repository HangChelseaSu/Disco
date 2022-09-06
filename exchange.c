
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "paul.h"
#include "profiler.h"

struct cell_lite{
   double prim[NUM_Q];
   double cons[NUM_Q];
   double RKcons[NUM_Q];
   double piph;
   double dphi;
   double wiph;
   double Phi[NUM_FACES];
   double RK_Phi[NUM_FACES];
};

#if USE_MPI
void generate_mpi_cell( MPI_Datatype * cell_mpi ){

   struct cell_lite test;
   int count = 8;
   int blocksize[]      = {NUM_Q,NUM_Q,NUM_Q,1,1,1,NUM_FACES,NUM_FACES};
   MPI_Datatype types[] = {MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE,MPI_DOUBLE};
   MPI_Aint offsets[8];

   offsets[0] = (char *)&(test.prim)   - (char *)(&test);
   offsets[1] = (char *)&(test.cons)   - (char *)(&test);
   offsets[2] = (char *)&(test.RKcons) - (char *)(&test);
   offsets[3] = (char *)&(test.piph)   - (char *)(&test);
   offsets[4] = (char *)&(test.dphi)   - (char *)(&test);
   offsets[5] = (char *)&(test.wiph)   - (char *)(&test);
   offsets[6] = (char *)&(test.Phi)    - (char *)(&test);
   offsets[7] = (char *)&(test.RK_Phi) - (char *)(&test);

   MPI_Type_create_struct( count , blocksize , offsets , types , cell_mpi );
   MPI_Type_commit( cell_mpi );

}

void copy_cell_to_lite( struct domain *theDomain, int jk, int i,
                        struct cell_lite * cl ){

   double *prim = &(theDomain->prim[jk][NUM_Q*i]);
   double *cons = &(theDomain->cons[jk][NUM_Q*i]);
   double *RKcons = &(theDomain->RKcons[jk][NUM_Q*i]);
   double *Phi = &(theDomain->Phi[jk][NUM_FACES*i]);
   double *RK_Phi = &(theDomain->RK_Phi[jk][NUM_FACES*i]);
  
   memcpy( cl->prim   , prim   , NUM_Q*sizeof(double) ); 
   memcpy( cl->cons   , cons   , NUM_Q*sizeof(double) ); 
   memcpy( cl->RKcons , RKcons , NUM_Q*sizeof(double) );
   memcpy( cl->Phi    , Phi    , NUM_FACES*sizeof(double) );
   memcpy( cl->RK_Phi , RK_Phi , NUM_FACES*sizeof(double) );
   cl->piph   = theDomain->piph[jk][i];
   cl->dphi   = theDomain->dphi[jk][i];
   cl->wiph   = theDomain->wiph[jk][i]; 

}

void copy_lite_to_cell( struct cell_lite * cl , struct domain * theDomain,
                        int jk, int i){ 

   double *prim = &(theDomain->prim[jk][NUM_Q*i]);
   double *cons = &(theDomain->cons[jk][NUM_Q*i]);
   double *RKcons = &(theDomain->RKcons[jk][NUM_Q*i]);
   double *Phi = &(theDomain->Phi[jk][NUM_FACES*i]);
   double *RK_Phi = &(theDomain->RK_Phi[jk][NUM_FACES*i]);

   memcpy( prim   , cl->prim   , NUM_Q*sizeof(double) ); 
   memcpy( cons   , cl->cons   , NUM_Q*sizeof(double) ); 
   memcpy( RKcons , cl->RKcons , NUM_Q*sizeof(double) );
   memcpy( Phi    , cl->Phi    , NUM_FACES*sizeof(double) );
   memcpy( RK_Phi , cl->RK_Phi , NUM_FACES*sizeof(double) );
   theDomain->piph[jk][i]   = cl->piph;
   theDomain->dphi[jk][i]   = cl->dphi;
   theDomain->wiph[jk][i]   = cl->wiph;

}

void generate_sendbuffer( struct domain * theDomain , int rnum , int znum , int dim , int * nijk , int * indexL , int * indexR , struct cell_lite * pl , struct cell_lite * pr , int dn1 , int dn2 , int mode ){

   int Periodic;
   if( dim == 0 )
       Periodic = theDomain->theParList.R_Periodic;
   else
       Periodic = theDomain->theParList.Z_Periodic;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int i,j,k;

   int iL = 0;
   int iR = 0;
   for( k=0 ; k<znum ; ++k ){
      nijk[1]=k;
      for( j=0 ; j<rnum ; ++j ){
         nijk[0]=j;
         nijk[dim] += dn1;

         int jk = nijk[0]+Nr*nijk[1];
         for( i=0 ; i<Np[jk] ; ++i ){
            if( mode==1 ){
               copy_cell_to_lite( theDomain, jk, i , pl+iL );
            }else if( mode==2 && ( dim_rank[dim] != 0 || Periodic ) ){
               copy_lite_to_cell( pl+iL , theDomain, jk, i );
            }
            ++iL;
         }

         nijk[dim] += dn2-dn1;

         jk = nijk[0]+Nr*nijk[1];
         for( i=0 ; i<Np[jk] ; ++i ){
            if( mode==1 ){
               copy_cell_to_lite( theDomain, jk, i , pr+iR );
            }else if( mode==2 && ( dim_rank[dim] != dim_size[dim]-1 || Periodic ) ){
               copy_lite_to_cell( pr+iR , theDomain, jk, i );
            }
            ++iR;
         }

         nijk[dim] -= dn2;
      }
   }

   *indexL = iL;
   *indexR = iR;
}

void resize_cell_data(struct domain *theDomain, int jk, int Np)
{
    if(Np == theDomain->Np[jk])
        return;

    int p;
    for(p=0; p<theDomain->N_data; p++)
    {
        double **field = theDomain->data[p];
        int size = theDomain->data_len[p];
        field[jk] = (double *)realloc(field[jk], Np * size * sizeof(double));
    }

    theDomain->Np[jk] = Np;
}

void generate_intbuffer( struct domain * theDomain , int rnum , int znum , int dim , int * nijk , int * indexL , int * indexR , int * Npl , int * Npr , int dn1 , int dn2 , int mode ){

   int Periodic;
   if( dim == 0 )
       Periodic = theDomain->theParList.R_Periodic;
   else
       Periodic = theDomain->theParList.Z_Periodic;
   int * dim_rank = theDomain->dim_rank;
   int * dim_size = theDomain->dim_size;
   int Nr = theDomain->Nr;
   int * Np = theDomain->Np;
   int j,k; 

   int iL = 0;
   int iR = 0; 
   for( k=0 ; k<znum ; ++k ){
      nijk[1]=k;
      for( j=0 ; j<rnum ; ++j ){
         nijk[0]=j;
         nijk[dim] += dn1; 

         int jk = nijk[0]+Nr*nijk[1];
         if( mode==1 ){
            Npl[ iL ] = Np[jk];
         }else if( mode==2 
                 && ( dim_rank[dim] != 0 || Periodic ) ){ 
            resize_cell_data(theDomain, jk, Npl[iL]);
         }    
         ++iL;

         nijk[dim] += dn2-dn1;

         jk = nijk[0]+Nr*nijk[1];
         if( mode==1 ){
            Npr[ iR ] = Np[jk];
         }else if( mode==2 
                 && ( dim_rank[dim] != dim_size[dim]-1 || Periodic ) ){
            resize_cell_data(theDomain, jk, Npr[iR]);
         }    
         ++iR;

         nijk[dim] -= dn2; 

      }    
   }
   *indexL = iL; 
   *indexR = iR;

}
#endif

void exchangeData( struct domain * theDomain , int dim ){
#if USE_MPI
   MPI_Datatype cell_mpi = {0}; 
   generate_mpi_cell( &cell_mpi );

   MPI_Comm grid_comm = theDomain->theComm;
   int * left_rank = theDomain->left_rank;
   int * right_rank = theDomain->right_rank;
   int Ng = theDomain->Ng;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;

   int tag = 0;
   MPI_Status status;
   int nijk[2];
   int rnum = Nr;
   int znum = Nz;
   int NN;

   if( dim == 0 ){
      NN = Nr;
      rnum = Ng;
   }else{
      NN = Nz;
      znum = Ng;
   }

   int indexL,indexR;
   int send_sizeL = 0;
   int send_sizeR = 0;
   int recv_sizeL = 0;
   int recv_sizeR = 0;

////////////////
//Send Np[jk]...
////////////////

   prof_tick(theDomain->prof, PROF_EXCH_NP_COUNT1);

//Count the number of Np's to send...
   generate_intbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , NULL    , NULL    , Ng , NN-2*Ng , 0 );
   prof_tock(theDomain->prof, PROF_EXCH_NP_COUNT1);
   prof_tick(theDomain->prof, PROF_EXCH_NP_COMM1);

   send_sizeL = indexL;
   send_sizeR = indexR;
//Tell your neighbor how many to expect...
   MPI_Sendrecv( &send_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag   , 
                 &recv_sizeR , 1 , MPI_INT , right_rank[dim] , tag  , grid_comm , &status);
   MPI_Sendrecv( &send_sizeR , 1 , MPI_INT , right_rank[dim] , tag+1 , 
                 &recv_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag+1, grid_comm , &status);
   
   prof_tock(theDomain->prof, PROF_EXCH_NP_COMM1);
   prof_tick(theDomain->prof, PROF_EXCH_NP_COUNT2);

   int Nl_send[send_sizeL];
   int Nr_send[send_sizeR];
   int Nl_recv[recv_sizeL];
   int Nr_recv[recv_sizeR];
//Build up list of ints to send...
   generate_intbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , Nl_send , Nr_send , Ng , NN-2*Ng , 1 );
   prof_tock(theDomain->prof, PROF_EXCH_NP_COUNT2);
   prof_tick(theDomain->prof, PROF_EXCH_NP_COMM2);
//Send!
   MPI_Sendrecv( Nl_send , send_sizeL , MPI_INT ,  left_rank[dim] , tag+2 ,
                 Nr_recv , recv_sizeR , MPI_INT , right_rank[dim] , tag+2, grid_comm , &status);
   MPI_Sendrecv( Nr_send , send_sizeR , MPI_INT , right_rank[dim] , tag+3 ,
                 Nl_recv , recv_sizeL , MPI_INT ,  left_rank[dim] , tag+3, grid_comm , &status);
   prof_tock(theDomain->prof, PROF_EXCH_NP_COMM2);
   prof_tick(theDomain->prof, PROF_EXCH_NP_FIN);
//Now take the list of ints and put them where they belong...
   generate_intbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , Nl_recv , Nr_recv , 0  , NN-Ng   , 2 );
   
   prof_tock(theDomain->prof, PROF_EXCH_NP_FIN);

////////////
//Send Cells
////////////

   prof_tick(theDomain->prof, PROF_EXCH_PREP);

//Count the number of cell_lites to send...
   generate_sendbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , NULL    , NULL    , Ng , NN-2*Ng , 0 );
   send_sizeL = indexL;
   send_sizeR = indexR;
//Tell you neighbor how many cells to expect...
   MPI_Sendrecv( &send_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag   ,
                 &recv_sizeR , 1 , MPI_INT , right_rank[dim] , tag  , grid_comm , &status);
   MPI_Sendrecv( &send_sizeR , 1 , MPI_INT , right_rank[dim] , tag+1 ,
                 &recv_sizeL , 1 , MPI_INT ,  left_rank[dim] , tag+1, grid_comm , &status);

   // Old version failed in some 3D MHD cases when allocation became too
   // large for the stack.  Changed to heap (ie. malloc() ) to fix this.
   //struct cell_lite pl_send[send_sizeL];
   //struct cell_lite pr_send[send_sizeR];
   //struct cell_lite pl_recv[recv_sizeL];
   //struct cell_lite pr_recv[recv_sizeR];
   struct cell_lite *pl_send = (struct cell_lite *)malloc(
                                        sizeof(struct cell_lite) * send_sizeL);
   struct cell_lite *pr_send = (struct cell_lite *)malloc(
                                        sizeof(struct cell_lite) * send_sizeR);
   struct cell_lite *pl_recv = (struct cell_lite *)malloc(
                                        sizeof(struct cell_lite) * recv_sizeL);
   struct cell_lite *pr_recv = (struct cell_lite *)malloc(
                                        sizeof(struct cell_lite) * recv_sizeR);
//Build up list of cells to send...
   generate_sendbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , pl_send , pr_send , Ng , NN-2*Ng , 1 );

   prof_tock(theDomain->prof, PROF_EXCH_PREP);
   prof_tick(theDomain->prof, PROF_EXCH_COMM);

//Send!
   MPI_Sendrecv( pl_send , send_sizeL , cell_mpi ,  left_rank[dim] , tag+2 ,
                 pr_recv , recv_sizeR , cell_mpi , right_rank[dim] , tag+2, grid_comm , &status);
   MPI_Sendrecv( pr_send , send_sizeR , cell_mpi , right_rank[dim] , tag+3 ,
                 pl_recv , recv_sizeL , cell_mpi ,  left_rank[dim] , tag+3, grid_comm , &status);

   prof_tock(theDomain->prof, PROF_EXCH_COMM);
   prof_tick(theDomain->prof, PROF_EXCH_FIN);

//Now take the list of cells and put them into the appropriate locations...
   generate_sendbuffer( theDomain , rnum , znum , dim , nijk , &indexL , &indexR , pl_recv , pr_recv , 0  , NN-Ng   , 2 );
   
   prof_tock(theDomain->prof, PROF_EXCH_FIN);

   free(pl_send);
   free(pr_send);
   free(pl_recv);
   free(pr_recv);

   MPI_Type_free( &cell_mpi );
#endif
}

void exchangePlanets(struct domain *theDomain)
{
    /*
     * This reduces (sums) the gas_track integrals to allow live planet motion
     * to work.
     */
#if USE_MPI

    int buf_size = theDomain->Npl * NUM_PL_INTEGRALS;

    MPI_Allreduce(MPI_IN_PLACE, theDomain->pl_gas_track, buf_size,
                  MPI_DOUBLE, MPI_SUM, theDomain->theComm);

#endif
    theDomain->planet_gas_track_synced = 1;
}

