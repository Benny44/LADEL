#include "types.h"
#include "constants.h"

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
            for (lca = prev_leaf; lca != ancestor[lca]; lca = ancestor[lca]);
            
            ladel_int last_path_node, ancestor_of_last_path_node;
            for (last_path_node = prev_leaf; last_path_node != lca; last_path_node = ancestor_of_last_path_node)
            {
                ancestor_of_last_path_node = ancestor[last_path_node];
                ancestor[last_path_node] = lca;
            }
            return lca;
        }  
    }
    
}

ladel_int ladel_col_counts(ladel_sparse_matrix *M, ladel_int *etree, ladel_int *postorder, ladel_int *col_counts)
{
    
}
