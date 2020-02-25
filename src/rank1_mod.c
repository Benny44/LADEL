#include "types.h"
#include "global.h"
#include "constants.h"
#include "rank1_mod.h"

ladel_int ladel_add_nonzero_pattern_to_L(ladel_int *Li, ladel_double *Lx, ladel_int *nz_L, ladel_int nzmax_L, ladel_int *set, ladel_int size_set)
{

}

ladel_set *ladel_init_set(ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set)
{
    ladel_set *set = ladel_malloc(1, sizeof(ladel_set));
    set->set = set_vals;
    set->size_set = size_set;
    set->max_size_set = max_size_set;
}

ladel_int ladel_set_union(ladel_set *first_set, ladel_set *second_set)
{
    ladel_int *set1 = first_set->set; 
    ladel_int size_set1 = first_set->size_set; 
    ladel_int max_size_set1 = first_set->max_size_set;
    ladel_int *set2 = second_set->set;
    ladel_int size_set2 = second_set->size_set;

    ladel_int changed = SET_HAS_NOT_CHANGED;
    ladel_int index1 = 0, index2, row1 = set1[0], row2, index, prev_row, temp;
    for (index2 = 0; index2 < size_set2; index2++)
    {
        row2 = set2[index2];
        for (; index1 < size_set1 && row1 < row2; index1++) 
        {
            row1 = set1[index1];
            if (row1 >= row2) break;
        }
        if (row1 < row2) /*append the rest of set2 to the end of set1*/
        {
            changed = SET_HAS_CHANGED; 
            for (; index2 < size_set2; index2++, index1++)
            {
                if (index1 == max_size_set1) return MAX_SET_SIZE_EXCEEDED;
                set1[index1] = set2[index2];
            }
            size_set1 = index1;
        } else if (row1 == row2) ;/*no addition needed*/
        else /*add elem of set2 in place index1 in set1*/
        {
            changed = SET_HAS_CHANGED;
            (size_set1)++;
            if (size_set1 == max_size_set1) return MAX_SET_SIZE_EXCEEDED;
            prev_row = row1;
            set1[index1] = row2;
            for (index = index1+1; index < size_set1; index++)
            {
                temp = set1[index];
                set1[index] = prev_row;
                prev_row = temp;
            }
        }
    }
    first_set->size_set = size_set1;
    return changed;
}


ladel_int ladel_rank1_update(ladel_factor **LD, ladel_symbolics **sym, ladel_sparse_matrix *W, ladel_int col_in_W)
{

}