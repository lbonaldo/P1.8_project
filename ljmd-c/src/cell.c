#include <cell.h>
#include <stddef.h>
#include <stdlib.h>

/* to allocate the list of pairs */
void allocate_pairs(mdsys_t *sys){

  int i, j, k , m;
  int iprec, inext, jprec, jnext, kprec, knext;
  int cidx, cstart;
  int N = sys->ncells;
  int *plist = sys->plist;
  int tmp[27]; //temporary plist for a single cell

  for(k=0; k<N; k++){
    kprec = (N+k-1)%N; knext = (k+1)%N;
    for(i=0; i<N; i++){
      iprec = (N+i-1)%N; inext = (i+1)%N;
      for(j=0; j<N; j++){
	jprec = (N+j-1)%N; jnext = (j+1)%N;
	cidx = j + i*N + k*N*N;
	//allocate the 27 cells
	tmp[0] = cidx;
	tmp[1] = cidx - j + jprec;
	tmp[2] = cidx - j + jnext;
	for(m=0; m<3; m++){
	  tmp[m+3] = tmp[m] - i*N + iprec*N;
	  tmp[m+6] = tmp[m] - i*N + inext*N;
	}
	for(m=0; m<9; m++){
	  tmp[m+9] = tmp[m] - k*N*N + kprec*N*N;
	  tmp[m+18] = tmp[m] - k*N*N + knext*N*N;
	}
	//allocate the 27 cells
	cstart = cidx*27;
	for(m=0; m<27; m++){
	  plist[2*(cstart+m)] = cidx;
	  plist[2*(cstart+m)+1] = tmp[m];
	}
      }
    }
  }
}

/* to allocate atoms positions inside the cells*/
void allocate_cells(mdsys_t *sys){
  
  int i, j, A;
  int c[3], cidx; //cell indices
  double linv, N;
  double r[3]; //array of positions;

  empty_cells(sys);
  
  for(i=0; i<sys->natoms; i++){
    r[0] = sys->rx[i];
    r[1] = sys->ry[i];
    r[2] = sys->rz[i];
    //renormalize positions
    for(j=0; j<3; j++){
      while(r[j] < 0)
	r[j] += sys->box;
      while(r[j] >= sys->box)
	r[j] -= sys->box;
    }

    //allocate the corresponding cell
    N = sys->ncells; linv = 1/(sys->clen);
    for(j=0; j<3; j++){
      c[j] = r[j]*linv;
      if(r[j]*linv > N)
	c[j] -= 1;
    }
    cidx = c[2]*N*N + c[0]*N + c[1];
    A = sys->catoms[cidx] + 1;

    if( sizeof(sys->clist[cidx])/sizeof(int) < A )
      sys->clist[cidx] = (int *)realloc(sys->clist[cidx], A*sizeof(int));

    sys->clist[cidx][A-1] = i;
    sys->catoms[cidx] += 1;
  }
}

void empty_cells(mdsys_t *sys){
  int i, j;
  for(i=0; i<sys->Ncells; i++){
    sys->catoms[i] = 0;
    for(j=0; j<sizeof(sys->clist[i])/sizeof(int); j++)
      sys->clist[i][j] = 0;
  }
}
