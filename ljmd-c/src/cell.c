#include <cell.h>
#include <stddef.h>

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
	//allocate the 26 cells
	cstart = cidx*52;
	for(m=0; m<26; m++){
	  plist[cstart+m] = cidx;
	  plist[cstart+m+1] = tmp[m+1];
	}
      }
    }
  }
}

/* to allocate atoms positions inside the cells*/
void allocate_cells(mdsys_t *sys){
  
  int i, j;
  int c[3], cidx; //cell indices
  double linv, N;
  double r[3]; //array of positions;

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
    if(sys->clist[cidx]==NULL){
      cell_t tmp;
      tmp.natoms=1;
      tmp.idxlist=(int*)malloc(sizeof(int));
      tmp.idxlist[0]=i;
      sys->clist[cidx]=&tmp;
      //printf("%d\n",sys->clist[cidx]->natoms);
      //sys->clist[cidx]->idxlist[sys->clist[cidx]->natoms] = i;
      //sys->clist[cidx]->natoms += 1;
      //printf("%d\t%d\n",i,cidx);
    }
    else{
       printf("%d\n",sys->clist[cidx]->natoms);
      cell_t tmp;
      tmp.natoms=(sys->clist[cidx]->natoms) + 1;
      tmp.idxlist=(int*)malloc(tmp.natoms*sizeof(int));
      //for(j=0; j<tmp.natoms-3; j++)
      //printf("%d\n",sys->clist[cidx]->idxlist[j]);
       	//tmp.idxlist[j]=sys->clist[cidx]->idxlist[j];
      //tmp.idxlist[tmp.natoms-2]=i;
      //free(sys->clist[cidx]);
      //sys->clist[cidx]=&tmp;
      //printf("%d\t%d\n",i,cidx);
    }
  }
}
