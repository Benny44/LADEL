#include "types.h"
#include "global.h"

#define STACK_NOT_EMPTY (top >= 0)

static ladel_int ladel_depth_first_search(  ladel_int root, 
                                            ladel_int post_index, 
                                            ladel_int *first_child, 
                                            ladel_int *next_child, 
                                            ladel_int *stack, 
                                            ladel_int *postorder)
{
    ladel_int current_node, current_child, top = 0;
    
    stack[top] = root;
    while (STACK_NOT_EMPTY)
    {
        current_node = stack[top];
        current_child = first_child[current_node];
        if (current_child == NONE)
        {
            top--;
            postorder[post_index] = current_node;
            post_index++;
        } else
        {
            first_child[current_node] = next_child[current_child];
            top++;
            stack[top] = current_child;
        }    
    }
}

ladel_int ladel_postorder(ladel_sparse_matrix *M, ladel_symbolics *sym)
{
    ladel_int *etree = sym->etree, *postorder = sym->postorder;
    if (!M || !etree) return FAIL;
    ladel_int col, ncol = M->ncol, top, prev_root = NONE, post_index = 0;
    ladel_int *work, *first_child, *next_child, *stack, *roots;
    work = ladel_malloc(ncol*4, sizeof(ladel_int));
    if (!work) return FAIL;
    first_child = work;
    next_child = work + ncol;
    stack = work + 2*ncol;
    roots = work + 3*ncol;

    for (col = 0; col < ncol; col++) first_child[col] = NONE;
    for (col = ncol-1; col >= 0; col--)
    {
        if (IS_ROOT(col, etree))
        {
            roots[col] = prev_root;
            prev_root = col;
        } else
        {
            next_child[col] = first_child[etree[col]];
            first_child[etree[col]] = col;
        }
    }
    for (col = prev_root; col != NONE; col = roots[col])
    {
        post_index = ladel_depth_first_search(col, post_index, first_child, next_child, stack, postorder);
    }

    ladel_free(work);
    return SUCCESS;
}