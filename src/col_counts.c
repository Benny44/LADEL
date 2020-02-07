#include "types.h"
#include "constants.h"
#include "transpose.h"
#include "global.h"

#define FIRST_LEAF 1
#define SUBSEQUENT_LEAF 2
#define NOT_A_LEAF -1

ladel_int ladel_least_common_ancestor(ladel_int subtree_root, ladel_int node, ladel_int *first_descendant,
                                        ladel_int *max_first_descendant, ladel_int *prev_leaf, 
                                        ladel_int *ancestor, ladel_int *node_type_of_leaf)
{
    if (subtree_root <= node || first_descendant[node] <= max_first_descendant[subtree_root])
    {
        *node_type_of_leaf = NOT_A_LEAF;
        return NONE;
    } else
    {
        max_first_descendant[subtree_root] = first_descendant[node];
        
        ladel_int last_leaf = prev_leaf[subtree_root];
        prev_leaf[subtree_root] = node;
        if (last_leaf == NONE)
        {
            *node_type_of_leaf = FIRST_LEAF;
            return node;
        } else
        {
            *node_type_of_leaf = SUBSEQUENT_LEAF;
            
            ladel_int lca; /*least common ancestor*/
            for (lca = last_leaf; lca != ancestor[lca]; lca = ancestor[lca]);
            
            ladel_int last_path_node, ancestor_of_last_path_node;
            for (last_path_node = last_leaf; last_path_node != lca; last_path_node = ancestor_of_last_path_node)
            {
                ancestor_of_last_path_node = ancestor[last_path_node];
                ancestor[last_path_node] = lca;
            }
            return lca;
        }  
    }
    
}

ladel_int ladel_col_counts(ladel_sparse_matrix *M, ladel_symbolics *sym)
{
    ladel_int *etree = sym->etree, *postorder = sym->postorder, *col_counts = sym->col_counts;
    if (!M || !etree || !postorder || !col_counts) return FAIL;
    
    ladel_int ncol = M->ncol, index, node, post_node, subtree_root, lca, type_of_leaf = NOT_A_LEAF;
    ladel_sparse_matrix *M_lower;
    if (M->symmetry == UPPER) M_lower = ladel_transpose(M, FALSE);
    else if (M->symmetry == LOWER) M_lower = M;
    else return FAIL;

    ladel_int *work = ladel_malloc(4*ncol, sizeof(ladel_int));
    if (!M_lower || !work)
    {
        if (M->symmetry == UPPER) ladel_sparse_free(M_lower);
        ladel_free(work);
        return FAIL;
    }

    ladel_int *prev_leaf = work, *first_descendant = work + ncol, *max_first_descendant = work + 2*ncol, *ancestor = work + 3*ncol;
    for (index = 0; index < 3*ncol; index++) work[index] = NONE;
    for (index = 0; index < ncol; index++) ancestor[index] = index;

    for (node = 0; node < ncol; node++)
    {
        post_node = postorder[node];
        if (first_descendant[post_node] == NONE) 
            col_counts[post_node] = 1;
        else 
            col_counts[post_node] = 0;
        for (; post_node != NONE && first_descendant[post_node] == NONE; post_node = etree[post_node]) 
            first_descendant[post_node] = node;
    }
    for (node = 0; node < ncol; node++)
    {
        post_node = postorder[node];
        if (!IS_ROOT(post_node, etree)) col_counts[etree[post_node]]--;
        for (index = M_lower->p[post_node]; index < M_lower->p[post_node+1]; index++)
        {
            subtree_root = M_lower->i[index];
            lca = ladel_least_common_ancestor(subtree_root, post_node, first_descendant, max_first_descendant, prev_leaf, ancestor, &type_of_leaf);
            if (type_of_leaf != NOT_A_LEAF) col_counts[post_node]++;
            if (type_of_leaf == SUBSEQUENT_LEAF) col_counts[lca]--; /*correct for duplicates*/
        }
        if (!IS_ROOT(post_node, etree)) ancestor[post_node] = etree[post_node];
    }
    for (node = 0; node < ncol; node++)
        if (!IS_ROOT(node, etree)) col_counts[etree[node]] += col_counts[node];

    for (node = 1; node < ncol; node++)
    {
        
        col_counts[node] += col_counts[node-1];
    }
        
    if (M->symmetry == UPPER) ladel_sparse_free(M_lower);
    ladel_free(work);
    return col_counts[ncol-1];
}
