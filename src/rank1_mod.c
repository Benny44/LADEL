#include "types.h"
#include "global.h"
#include "constants.h"
#include "rank1_mod.h"

ladel_int ladel_add_nonzero_pattern_to_L(ladel_col *col, ladel_set *set)
{

}

ladel_set *ladel_init_set(ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set)
{
    ladel_set *set = ladel_malloc(1, sizeof(ladel_set));
    set->set = set_vals;
    set->size_set = size_set;
    set->max_size_set = max_size_set;
}

ladel_col *ladel_init_col(ladel_int *i, ladel_double *x, ladel_int nz , ladel_int nzmax)
{
    ladel_col *col = ladel_malloc(1, sizeof(ladel_col));
    col->i = i;
    col->x = x;
    col->nz = nz;
    col->nzmax = nzmax;
}

void ladel_set_col(ladel_col* col,ladel_int *i, ladel_double *x, ladel_int nz , ladel_int nzmax)
{
    col->i = i;
    col->x = x;
    col->nz = nz;
    col->nzmax = nzmax;
}

ladel_int ladel_set_union(ladel_set *first_set, ladel_set *second_set, ladel_set *difference, ladel_int* offset, ladel_int* insertions)
{
    ladel_int *set1 = first_set->set; 
    ladel_int size_set1 = first_set->size_set; 
    ladel_int max_size_set1 = first_set->max_size_set;
    ladel_int *set2 = second_set->set;
    ladel_int size_set2 = second_set->size_set;

    ladel_int *dif = difference->set;
    // difference->size_set = 0;

    ladel_int index1 = 0, index2, row1 = set1[0], row2, index, index_dif = 0;

    /* Construct difference set and offsets -----------------------------------*/ 
    for (index2 = 0; index2 < size_set2; index2++)
    {
        row2 = set2[index2];
        for (; index1 < first_set->size_set && row1 < row2; index1++) 
        {
            row1 = set1[index1];
            offset[index1] = index_dif;
            if (row1 >= row2) break; 
        }
        if (row1 > row2)/*add elem of set2 in place index1 in set1*/
        {
            dif[index_dif] = row2;
            index_dif++;

            size_set1++;
            if (size_set1 == max_size_set1) return MAX_SET_SIZE_EXCEEDED;
        }
        else if (row1 < row2) /*append the rest of set2 to the end of set1*/
        {
            for (; index2 < size_set2; index2++, size_set1++, index_dif++)
            {
                if (size_set1 == max_size_set1) return MAX_SET_SIZE_EXCEEDED;
                dif[index_dif] = set2[index2];
                insertions[index_dif] = index1+index_dif;
            }
        }
    }
    
    if (index_dif == 0) return SET_HAS_NOT_CHANGED;
    
    for (; index1 < first_set->size_set; index1++) offset[index1] = index_dif;
    difference->size_set = index_dif;

    /* Merge difference into the first set ------------------------------------*/

    /* Move original set1 values to the correct positions */
    for (index1 = first_set->size_set-1; index1 >= 0; index1--) set1[index1+offset[index1]] = set1[index1];

    /* Compute the positions for the new elements */
    index_dif = 0;
    for (index1 = 0; index1 < first_set->size_set; index1++)
        for (; index_dif < offset[index1]; index_dif++) 
            insertions[index_dif] = index1 + index_dif;

    /* Insert the new elements */
    for (index_dif = 0; index_dif < difference->size_set; index_dif++) set1[insertions[index_dif]] = dif[index_dif];
    
    first_set->size_set = size_set1;
    return SET_HAS_CHANGED;
}


ladel_int ladel_rank1_update(ladel_factor **LD, ladel_symbolics **sym, ladel_sparse_matrix *W, ladel_int col_in_W)
{

}