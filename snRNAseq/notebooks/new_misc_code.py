import scanpy as sc
import numpy as np
import pandas as pd
import pickle
import matplotlib.pyplot as plt

###### Functions used to generate distinct list of hex colors ######

import random
import matplotlib as mpl

def get_random_color( pastel_factor=0.5):
    return [( x+pastel_factor)/( 1.0+pastel_factor) for x in [random.uniform( 0,1.0) for i in [1,2,3]]]

def color_distance( c1,c2):
    return sum( [abs( x[0]-x[1]) for x in zip( c1,c2)])

def generate_new_color( existing_colors, pastel_factor=0.5):
    max_distance = None
    best_color = None
    for i in range( 0,7_000):
        color = get_random_color( pastel_factor=pastel_factor)
        if not existing_colors:
            return color
        best_distance = min( [color_distance( color,c) for c in existing_colors])
        if not max_distance or best_distance > max_distance:
            max_distance = best_distance
            best_color = color
    return best_color

def random_hex_colors( n_colors, random_seed=0):    #Example:
    #To make your color choice reproducible, uncomment the following line:
    random.seed( random_seed) # is what we used for GABA fig 12345678
    colors = []
    for i in range( 0, n_colors):
        colors.append( generate_new_color( colors, pastel_factor=random.uniform( 0.0,0.10)))
    new_colors = np.array( [mpl.colors.to_hex(  ii) for ii in colors ])
    return( new_colors)


###### Misc Functions ######

# compare two sparse csr matrices for equality
def csr_matrix_equal(a1, a2):
    return( np.array_equal(a1.indptr, a2.indptr) and 
            np.array_equal(a1.indices, a2.indices) and 
            np.array_equal(a1.data, a2.data))

# match barcode order of adata1 to adata2 
def match_bc_order( adata1, adata2):
    bc_list1 = adata1.obs_names.values.tolist()
    bc_list2 = adata2.obs_names.values.tolist()
    ordered_args = np.array([bc_list2.index(ii) for ii in bc_list1])
    adata2 = adata2[ordered_args,:]
    return( adata2)

# match gene name order of adata1 to adata2 
def match_bc_order( adata1, adata2):
    bc_list1 = adata1.var_names.values.tolist()
    bc_list2 = adata2.var_names.values.tolist()
    ordered_args = np.array([bc_list2.index(ii) for ii in bc_list1])
    adata2 = adata2[:,ordered_args]
    return( adata2)

# test elements in A for membership in B
def member_test( A, B):
    b = set(B)
    return( [x in b for x in A])

def save_obj( obj, name):
    with open( name, 'wb') as f:
        pickle.dump( obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj( name):
    with open( name, 'rb') as f:
        return pickle.load(f)
    
    
######################################################
######     UMAT     ##################################
######################################################

def umat_neighbors( stage_order, adata, n_neighbors=5, n_pcs=50, rand=123):
    # if uns['neighbors'] doesnt exist run neighbors on whole lot
    if not('neighbors' in adata.uns):
        sc.pp.neighbors( adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=rand)
    # pass copy of adata if do not want neighbors changed
    n_cells = adata.shape[0]
    conn = np.zeros( (n_cells,n_cells), dtype=float)
    dist = np.zeros( (n_cells,n_cells), dtype=float)
    inds = np.arange( len(conn))
    # neighbors iteratively for adjacent stages
    for itr1, itr2 in zip( stage_order[:-1],stage_order[1:]):
        itr_mk = np.in1d( adata.obs['stage_ids'],[itr1,itr2])
        itr_adata = adata[itr_mk]
        sc.pp.neighbors( itr_adata, n_neighbors=n_neighbors, n_pcs=n_pcs, random_state=rand)
        con_itr = itr_adata.obsp['connectivities'].toarray()
        dis_itr = itr_adata.obsp['distances'].toarray()
        for ii, (c_itr, d_itr) in zip( inds[itr_mk], zip( con_itr, dis_itr)):
            conn[ii,itr_mk] = c_itr
            dist[ii,itr_mk] = d_itr
        print( itr1, itr2, itr_mk.sum())
    adata.obsp['connectivities'] = sp.sparse.csr_matrix( conn)
    adata.obsp['distances']      = sp.sparse.csr_matrix( dist)    
    return( adata)

def umat( stage_order, adata, n_neighbors=5, n_pcs=50, min_dist=0.5, n_comps=2, rand=123):
    adata = umat_neighbors( stage_order, adata, n_neighbors=n_neighbors, n_pcs=n_pcs, rand=rand)
    if ('X_umap' in adata.obsm):
        umap = adata.obsm['X_umap']
    # need to update function to include variables for umap
    sc.tl.umap( adata, random_state=rand, min_dist=min_dist, n_components=n_comps)
    adata.obsm['X_umat'] = adata.obsm['X_umap']
    if ('X_umap' in adata.obsm):
        adata.obsm['X_umap'] = umap
    else:
        del adata.obsm['X_umap']
    return( adata)
