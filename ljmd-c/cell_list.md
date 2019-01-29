# Cell List

## Scratch notes

1. Create a `struct _cell` which contains:

   - `int natoms` to save the total number of atoms; 
   - `int *idxlist` to store the atom indices belonging to the current cell.

2. Update `struct _syst` with:

   - `_cell *clist` to store the cell ID's;
   - `int *plist` to store the list of pairs;
   - `double *clen` to save the length of the cells;
   - `int ncells` to save the number of cells IN A SINGLE DIRECTION;
   - NOTICE: put an assert n the main to check `sys->clen >= sys->rcut `;

3. Compute the total number of cells;

   - divide `sys->box / sys->clen`;
   - if `remainder != 0` enlarge the last cell;
   - therefore `sys->ncells = int{sys->bos/sys->clen}`.

4. Prepare `sys->plist` with the list of close pairs:

   - each cell has `26` close cells;
   - if the cell is on the border, the close cell is the one on the other side; account for this case.

5. Create a function `allocate_cells` which allocates atoms into cells at every step (maybe I can split in two functions?):

   ```c
   void allocate_cells(_mdsys sys){
   	int i, j;
       int ix, iy, iz, cidx; //cell indices
   	double rx, ry, rz, linv, N;
   	double *r; //array of positions;
       _cell *cell; //cell pointer;
   	for(i=0; i < sys->natoms; i++){
       	rx = sys->rx[i]; r[0] = rx;
       	ry = sys->ry[i]; r[1] = ry;
       	rz = sys->rz[i]; r[2] = rz;
           //renormalize position
       	for(j=0; j<3; j++){
           	while(r[j] < 0)
               	r[j] += sys->box;
           	while(r[j] >= sys->box)
               	r[j] -= sys->box;
       	}
           //allocate into the corresponding cell
           N = sys->ncells; linv = 1/sys->clen;
           cx = rx*linv - _Bool{rx*linv > N};
           cy = ry*linv - _Bool{ry*linv > N}; 
           cz = rz*linv - _Bool{rz*linv > N};
           cidx = cz*N*N + cx*N + cy; //cell index
           cell = sys->clist[cidx]; //cell identifier in the struct
           cell->idxlist[cell->natoms + 1] = i;
           cell->natoms += 1;
   	}
   }
   ```
