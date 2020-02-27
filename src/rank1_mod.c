#include "types.h"
#include "global.h"
#include "constants.h"
#include "rank1_mod.h"

ladel_int ladel_add_nonzero_pattern_to_col_of_L(ladel_sparse_matrix* L, ladel_int col, ladel_set *col_set, ladel_set *set, ladel_set *difference, ladel_int* offset, ladel_int* insertions)
{
    ladel_int start = L->p[col], status;
    ladel_set_set(col_set, L->i + start, L->nz[col], L->p[col+1] - L->p[col]);
    status = ladel_set_union(col_set, set, difference, offset, insertions);
    
    /* For now it is assumed the user has allocated enough space. If not, the error is passed on.
    Technically, some reallocation code could go here to include cases where one doesn't know the
    maximum number of nonzeros in each column of L. */
    if (status == MAX_SET_SIZE_EXCEEDED) return MAX_SET_SIZE_EXCEEDED;
    
    /* Initialize new nonzero values in L to zero */
    ladel_int index;
    for (index = 0; index < difference->size_set; index++) L->x[start+insertions[index]] = 0;
    L->nz[col] = col_set->size_set;

    return status;
}

void ladel_set_set(ladel_set *set, ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set)
{
    set->set = set_vals;
    set->size_set = size_set;
    set->max_size_set = max_size_set;
}

ladel_set *ladel_init_set(ladel_int *set_vals, ladel_int size_set, ladel_int max_size_set)
{
    ladel_set *set = ladel_malloc(1, sizeof(ladel_set));
    ladel_set_set(set, set_vals, size_set, max_size_set);
}

void ladel_shallow_copy_set(ladel_set *set_out, ladel_set *set_in)
{
    set_out->set = set_in->set;
    set_out->size_set = set_in->size_set;
    set_out->max_size_set = set_in->max_size_set;
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


ladel_int ladel_rank1_update(ladel_factor *LD, ladel_symbolics *sym, ladel_sparse_matrix *W, ladel_int col_in_W)
{
    ladel_int *etree = sym->etree;
    ladel_sparse_matrix *L = LD->L;
    ladel_double *Dinv = LD->Dinv;

    /* TODO: account for the permutation! */
    
    ladel_int ncol = L->ncol, col, row, index, index_L, size_W = W->p[col_in_W+1] - W->p[col_in_W];
    ladel_int changed = SET_HAS_NOT_CHANGED, changed_W, changed_child;
    /* TODO: capture all the memory allocations in a ladel_work struct */
    ladel_set *set_W = ladel_init_set(W->i + W->p[col_in_W], size_W, size_W);
    ladel_set *set_L = ladel_malloc(1, sizeof(ladel_set));

    ladel_set *difference = ladel_malloc(1, sizeof(ladel_set));
    ladel_int *temp_ncol = ladel_malloc(ncol, sizeof(ladel_int));
    ladel_set *difference2 = ladel_init_set(temp_ncol, 0, ncol);
    ladel_set *difference_child = ladel_malloc(1, sizeof(ladel_set));
    ladel_int *offset = ladel_malloc(ncol, sizeof(ladel_int));
    ladel_int *insertions = ladel_malloc(ncol, sizeof(ladel_int));
    ladel_double *W_col = ladel_calloc(ncol, sizeof(ladel_double));
    for (index = W->p[col_in_W]; index < W->p[col_in_W+1]; index++) 
        W_col[W->i[index]] = W->x[index];

    ladel_double alpha, alpha_new, gamma, w, dinv;
    ladel_int child, old_parent;
    for (index = W->p[col_in_W]; index < W->p[col_in_W+1]; index++)
    {
        col = W->i[index];
        changed = ladel_add_nonzero_pattern_to_col_of_L(L, col, set_L, set_W, difference, offset, insertions);
        if (changed == MAX_SET_SIZE_EXCEEDED) return FAIL;

        /* numerical update */
        w = W_col[col];
        dinv = Dinv[col];
        alpha_new = alpha + w*w*dinv;
        gamma = w*dinv/alpha_new; /*if alpha_new == o then matrix not full rank */
        Dinv[col] = alpha_new/alpha;
        for (index_L = L->p[col]; index_L < L->p[col]+L->nz[col]; index_L++)
        {
            row = L->i[index_L];
            W_col[row] -= w*L->x[index_L];
            L->x[index_L] += gamma*W_col[row];
        }

        if (changed == SET_HAS_CHANGED)
        {
            child = col;
            old_parent = etree[col];
            col = etree[col] = L->i[L->p[col]]; /*new_parent, there must be one since the set has changed*/

            /* prepare the difference in the next col due to the child */
            if (col == old_parent)
                ladel_shallow_copy_set(difference_child, difference);
            else
                ladel_set_set(difference_child, L->i+L->p[child], L->nz[child], L->p[child+1] - L->p[child]);
            
            break;
        }
    }

    if (changed == SET_HAS_CHANGED)
    {
        while (TRUE)
        {
            changed_child = ladel_add_nonzero_pattern_to_col_of_L(L, col, set_L, difference_child, difference, offset, insertions);
            changed_W = ladel_add_nonzero_pattern_to_col_of_L(L, col, set_L, set_W, difference_child, offset, insertions);
            
            
            /* numerical update */
            w = W_col[col];
            dinv = Dinv[col];
            alpha_new = alpha + w*w*dinv;
            gamma = w*dinv/alpha_new; /*if alpha_new == o then matrix not full rank */
            Dinv[col] = alpha_new/alpha;
            for (index_L = L->p[col]; index_L < L->p[col]+L->nz[col]; index_L++)
            {
                row = L->i[index_L];
                W_col[row] -= w*L->x[index_L];
                L->x[index_L] += gamma*W_col[row];
            }

            child = col;
            old_parent = etree[col];
            if (L->nz[col] == 0) break; /* There is no subdiag entry, so no new parent */
            col = etree[col] = L->i[L->p[col]]; /*new_parent*/

            /* prepare the difference in the next col due to the child */
            if (col == old_parent)
            {
                /* NB difference_child is already the result of the change in W, therefore it does
                not need to be adjusted when only changed_W == SET_HAS_CHANGED */
                if (changed_W == SET_HAS_CHANGED && changed_child == SET_HAS_CHANGED)
                    ladel_set_union(difference_child, difference, difference2, offset, insertions);
                else if (changed_child == SET_HAS_CHANGED)
                    ladel_shallow_copy_set(difference_child, difference);
                else if (changed_W == SET_HAS_NOT_CHANGED)
                    difference_child->size_set = 0; /* difference_child is empty */
            }     
            else
                ladel_set_set(difference_child, L->i+L->p[child], L->nz[child], L->p[child+1] - L->p[child]);
            
             
        }
    }

}
