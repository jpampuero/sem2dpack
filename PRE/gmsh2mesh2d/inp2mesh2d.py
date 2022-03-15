"""
A python module that convert the Gmsh generated mesh to .mesh2d file 

Since mesh2d needs boundary conditions, physical boundaries must be defined 
in GMSH and the tag must be consistent with definition in the sem2dpack 
input file

mesh file can be generated from geofile by
gmsh geofile.geo -2 -o mesh -format xxformat -save_all

Chao Liang 11/24/2021 @Chengdu, Sichuan

Note:
** different boundaries must be defined in gmsh as physical curves 
   and labeled as bc1, bc2, bc3,...

** different material must be defined in gmsh as physical surface
   and labeled as mat1, mat2, ...

This module uses meshio library.
see https://github.com/nschloe/meshio, instruction for installation
"""

import meshio
import numpy as np
import sys
import os
import argparse 
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import reverse_cuthill_mckee
import argparse

def find_edge(lines, enod, coord):
    """
    Find the element number and edge number that lines belong to
    lines: [2, number of lines]
    enod : [4, number of elements], sorted
    coord: [2, number of nodes]
    
    return:
    elem_idx: [number of lines][2], 
    edge_idx: [number of lines][2]
    
    Note that one edge can be possessed by two elements
    """
    
    # compute the edge centers for each edge in each element
    edge_centers_x = (coord[0, enod] + coord[0, enod[[1,2,3,0],:]])/2.0
    edge_centers_y = (coord[1, enod] + coord[1, enod[[1,2,3,0],:]])/2.0
    
    # loop through each line find all the edge centers that matches the line center
    nl = lines.shape[1]
    elem_idx = np.ones((2, nl), dtype=np.int64)*-1
    edge_idx = np.ones((2, nl), dtype=np.int64)*-1
    
    for i in range(nl):
        line = lines[:, i]
        xc   = np.mean(coord[0, line])
        yc   = np.mean(coord[1, line])
        dist = (edge_centers_x-xc)**2 + (edge_centers_y - yc)**2
        idx  = np.argwhere(dist==0)
        nelem = idx.shape[0]
        # loop through each element and check both end points are in element
        cnt = 0
        for j in range(nelem):
            idxj = idx[j]
            if (line[0] in enod[:,idxj[1]] and line[1] in enod[:, idxj[1]]):
                elem_idx[cnt, i] = idx[j, 1]
                edge_idx[cnt, i] = idx[j, 0]
                cnt += 1

    return elem_idx, edge_idx

def order_enod(enod, coord):
    """
    Re-order the nodes in counter-clockwise direction
    """
    nel = enod.shape[1]
    enod_order = np.zeros(enod.shape, dtype=np.int64)
    for i in range(nel):
        ei = enod[:, i]
        enod_order[:, i] = order_nodes_single_element(ei, coord[0,:], coord[1,:])
    return enod_order

def  order_nodes_single_element (elem_nodes, x_nodes, z_nodes):
    """
    This is a subroutine borrowed from Elif, 
    it seems the code is reordering the nodes within each element in a counter-clockwise direction
    """
    # print '*Given element nodes   ', elem_nodes[0:4]
    # Checking jacobien (pour voir le sens de numerotation)
    eti1 = int(elem_nodes[0])
    eti2 = int(elem_nodes[1])
    eti3 = int(elem_nodes[2])
    eti4 = int(elem_nodes[3])

    xv1  = float(x_nodes[eti2])-float(x_nodes[eti1])
    xv2  = float(x_nodes[eti3])-float(x_nodes[eti2])
    zv1  = float(z_nodes[eti2])-float(z_nodes[eti1])
    zv2  = float(z_nodes[eti3])-float(z_nodes[eti2])
    jac = -xv2*zv1+ xv1*zv2

    new_nodes = np.zeros(4)
    # To modify this part for more combinations !
    if jac<0:
        # print 'JACOBIAN NEGATIVE !';print
        new_nodes[0] = elem_nodes[3]
        new_nodes[1] = elem_nodes[2]
        new_nodes[2] = elem_nodes[1]
        new_nodes[3] = elem_nodes[0]
        # check if stop exists in python
#        sys.exit()
    else:
        new_nodes[0] = elem_nodes[0]
        new_nodes[1] = elem_nodes[1]
        new_nodes[2] = elem_nodes[2]
        new_nodes[3] = elem_nodes[3]
    ###

    # print '*New element nodes     ', new_nodes[0:4]
    # print

    # Checking if the 1st point is the one we're looking for
    z_list = []
    z_list.append(float(z_nodes[int(new_nodes[0])]))
    z_list.append(float(z_nodes[int(new_nodes[1])]))
    z_list.append(float(z_nodes[int(new_nodes[2])]))
    z_list.append(float(z_nodes[int(new_nodes[3])]))
    z_list = sorted (z_list)

    zmin_loc1 = z_list[0]
    zmin_loc2 = z_list[1]

    # print 'Sorted list of z values            ', z_list
    # print '2 Minimum values on the element:   ', zmin_loc1, zmin_loc2

    count1 = 0; count2=0
    local1 = []; local2 =[]
    for i in range(4):
        zloc = float(z_nodes[  int(new_nodes[i])])
        if ( zloc ==  float(zmin_loc1)  ):
            count1 = count1+1
            local1.append( int(new_nodes[i])  )
            # print 'count1 zloc zmin1'
            # print count1, zloc, local1

        if ( zloc ==  float(zmin_loc2)      ):
            count2 = count2+1
            local2.append( int(new_nodes[i])  )
            # print 'count2 zloc zmin1'
            # print count2, zloc, local2            

        ##
    ###

    # If only 1 point has z_min
    if      count1 == 1   and count2==1:
        x1   = float(x_nodes[ int(local1[0])]) 
        x2   = float(x_nodes[ int(local2[0])])               
        if (x1 < x2):
                loc1 = local1[0]
        else:
                loc1 = local2[0]

    elif    count1 == 2:
        x1   = float(x_nodes[ int(local1[0])]) 
        x2   = float(x_nodes[ int(local1[1])]) 
        if (x1 < x2):
                loc1 = local1[0]
        else:
                loc1 = local1[1]
#modified by flomin
    elif    count2 == 2:
        x1   = float(x_nodes[ int(local2[0])])
        x2   = float(x_nodes[ int(local2[1])])
        loc1 = local1[0]
    else:
        print ('Too many vertices with same z-coordinate')
    ###

    new_nodes2 = np.zeros(4)

    # 1st case
    if      loc1 == int(new_nodes[1]):
        new_nodes2[0] = new_nodes[1]
        new_nodes2[1] = new_nodes[2]
        new_nodes2[2] = new_nodes[3]
        new_nodes2[3] = new_nodes[0]

    elif    loc1 == int(new_nodes[2]):
        new_nodes2[0] = new_nodes[2]
        new_nodes2[1] = new_nodes[3]
        new_nodes2[2] = new_nodes[0]
        new_nodes2[3] = new_nodes[1]

    elif    loc1 == int(new_nodes[3]):
        new_nodes2[0] = new_nodes[3]
        new_nodes2[1] = new_nodes[0]
        new_nodes2[2] = new_nodes[1]
        new_nodes2[3] = new_nodes[2]
    else:
        new_nodes2[0] = new_nodes[0]
        new_nodes2[1] = new_nodes[1]
        new_nodes2[2] = new_nodes[2]
        new_nodes2[3] = new_nodes[3]            
    ###
    return np.int64(new_nodes2)

def read_inp(filename, pairs):
    """
    script to read .inp file (2D, quad mesh).
    return coord, enode, etag, bnds
    
    coord: X,Y coordinates of nodes, [2, total number of nodes]
    enod:  Element connectivity table, [4, total number of elements]
    etag:  Domain or material tag (an integer) for each element
           Vector of length = total number of elements.
    bnds:  Boundaries database
           A list of size of number of boundaries.
           bnds[k](1, be) is the mesh element that contains 
                          the be-th boundary element of the k-th boundary
           bnds[k](2, be) is the index of the corresponding edge of the element.
           
    The numbering of the nodes and edges for each element is as follow
                
            4----------(3)------------3
            |                         |
            |                         |
           (4)                       (2)
            |                         |
            |                         |
            1----------(1)------------2
            
    """
    mesh  = meshio.read(filename)
    
    # X,Y coordinates
    coord = mesh.points.T[:2,:]
    
    # get all the element
    enod = mesh.cells_dict['quad'].T
    
    # reorder enod, zero based node index
    enod = order_enod(enod, coord)
    
    # get all the lines, zero based node index
    lines= mesh.cells_dict['line'].T
    
    # get the element tag
    nelem = enod.shape[1]
    etag  = np.ones(nelem, dtype=np.int32)
    
    # get all the physical surfaces/domains
    for k in mesh.cell_sets_dict:
        mesh.cell_sets_dict[k.lower]=mesh.cell_sets_dict.pop(k)
    
    quad_key_list = []
    for k in mesh.cell_sets_dict:
        if 'quad' in mesh.cell_sets_dict[k]:
            quad_key_list.append(k)
    quad_key_list.sort()
    
    nmat = len(quad_key_list)
    for i in range(nmat):
        mat_i = 'mat'+str(i+1)
        if quad_key_list[i] != mat_i:
            sys.exit('Wrong material name %s'%(quad_key_list[i]))
    
    # set element tags
    for i in range(nmat):
        mat_i = 'mat'+str(i+1)
        etag[mesh.cell_sets_dict[mat_i]['quad']] = i + 1
    
    # construct bnds
    # get all the physical lines as boundaries.
    bnds       = []
    bnds_lines = []
    sort_keys  = np.sort(list(mesh.cell_sets_dict.keys()))
    cnt        = 0
    for k in sort_keys:
        if 'line' in mesh.cell_sets_dict[k]:
            cnt = cnt + 1
            bc_i  = 'bc%d'%(cnt)
            line_idx   = mesh.cell_sets_dict[k]['line']
            bnds_lines = lines[:, line_idx]
            nline      = len(line_idx)
            bnd_i      = np.zeros((2, nline), dtype=np.int64)

            #find the element numbers and edge numbers for a list of lines
            elem_idx, edge_idx = find_edge(bnds_lines, enod, coord)
            #
            if np.any(elem_idx[1,:] != -1) or np.any(edge_idx[1,:] != -1):
                sys.exit('any least one line in boundary %s has two neighboring elements'%(k))
            bnd_i[0,:] = elem_idx[0,:] # use 0-based number
            bnd_i[1,:] = edge_idx[0,:]
            bnds.append(bnd_i)
    
    # so far enods and bnds are still 0-based numbering
    # rcm reordering to reduce bandwidth, reverse_cuthill_mckee
    r = rcm_reorder(enod, etag, bnds, coord, pairs)
    enod = enod[:, r]
    etag = etag[r]
    
    ################# redo the boundaries  ###################
    #if len(pairs)>0:
    bnds       = []
    bnds_lines = []
    sort_keys  = np.sort(list(mesh.cell_sets_dict.keys()))
    cnt        = 0
    for k in sort_keys:
        if 'line' in mesh.cell_sets_dict[k]:
            cnt = cnt + 1
            bc_i  = 'bc%d'%(cnt)
            line_idx   = mesh.cell_sets_dict[k]['line']
            bnds_lines = lines[:, line_idx]
            nline      = len(line_idx)
            bnd_i      = np.zeros((2, nline), dtype=np.int64)

            #find the element numbers and edge numbers for a list of lines
            elem_idx, edge_idx = find_edge(bnds_lines, enod, coord)
            #
            if np.any(elem_idx[1,:] != -1) or np.any(edge_idx[1,:] != -1):
                sys.exit('any least one line in boundary %s has two neighboring elements'%(k))
            bnd_i[0,:] = elem_idx[0,:] # use 0-based number
            bnd_i[1,:] = edge_idx[0,:]
            bnds.append(bnd_i)
    ################# DONE redo the boundaries  ###################
    
    # switch to 1 based numbering for both enod and bnds
    enod += 1
    for i in range(len(bnds)):
        bnds[i] += 1
    
    return coord, enod, etag, bnds
    
def write_mesh2d(filename, coord, enod, bnds, etag):
    """
    write mesh data into .mesh2d file 
    """
    if not filename.endswith('.mesh2d'):
        filename += '.mesh2d'
    f = open(filename,'w')
    print('Writing mesh2d mesh to %s'%(filename))
    f.write('NEL NPEL NNOD NBC\n');
    nel = enod.shape[1]
    npel= enod.shape[0]
    nnod= coord.shape[1]
    nbc = len(bnds)
    f.write('%u %u %u %u\n'%(nel, npel, nnod, nbc))
    f.write('NID X Y\n')
# write all the nodes
    for i in range(nnod):
        f.write('%u %f %f\n'%(i+1, coord[0,i], coord[1,i]))
    f.write('EID NODES TAG\n')
# write element connectivity 
    for i in range(nel):
        line   = [str(i+1)]
        enodi  = enod[:,i]
        for j in enodi:
            line.append(str(j))
        line.append(str(etag[i]))
        f.write("%s\n"%(" ".join(line)))
# write boundaries 
    for i in range(nbc):
        if bnds[i] is None:
            continue
        f.write('BCTAG NBEL\n')
        nbel = bnds[i].shape[1]
        f.write('%u %u\n'%(i+1, nbel))
        f.write('BID EID EDGE\n')
        for j in range(nbel):
            f.write('%u %u %u\n'%(j+1, bnds[i][0,j], bnds[i][1,j]))
    f.close()
    print('Done writing mesh2d mesh to %s'%(filename))
    return 
    
def inp2mesh2d(inpfile, mesh2dfile, pairs):
    print('Start reading %s'%(inpfile))
    coord, enod, etag, bnds = read_inp(inpfile, pairs)
    
    print('Done reading %s'%(inpfile))
    write_mesh2d(mesh2dfile, coord, enod, bnds, etag)
    return

def geo2inp(geofile):
    """
    generate an inp mesh from a geofile
    """
    fileroot = geofile
    inpfile  = fileroot[:-4] + '.inp'
    cmd   = 'gmsh %s -2 -o %s -format inp  -save_all'%(geofile, inpfile)
    os.system(cmd)
    return

def KVtag(bc, dist):
    """
    tagging elements that are KV material
    +1 to original tag value

    compute the minimum distance to each line
    if smaller than
    """

def FindNeigborConnection(enod):
    nel     = enod.shape[1]
    nnode   = np.max(enod) + 1
    node2el = []
    for i in range(nnode):
        node2el.append([])
    cnt = np.zeros(nnode, dtype=np.int)
    
    for i in range(nel):
        nodes = enod[:,i]
        for j in range(len(nodes)):
            i_node = nodes[j]
            if not (i in node2el[i_node]):
                cnt[i_node] += 1
                node2el[i_node].append(int(i))
    
    ne = []
    for i in range(nel):
        ne.append([])
        
    cnt2  = np.zeros(nel,dtype=np.int)
    elcon = []
    for i in range(nnode):
        n2e = node2el[i]
        cnt_i = cnt[i] # number of elements for each node
        for k in range(cnt_i):
            e_i = n2e[k]
            for j in range(k+1, cnt_i):
                e_j = n2e[j]
                if not (e_j in ne[e_i]):
                    cnt2[e_i] += 1
                    cnt2[e_j] += 1
                    ne[e_i].append(e_j)
                    ne[e_j].append(e_i)
                    
                    # append element connection pairs(e_i, e_j)
                    elcon.append([e_i, e_j])
                    elcon.append([e_j, e_i])
    elcon = np.array(elcon,dtype=np.int)
    max_neb= np.max(cnt2)
    return ne, max_neb, elcon

def rcm_reorder(enod, etag, bnds, coord, pairs):
    """
    reorder the element number using rcm algorithm to reduce the fill
    enod, bnds are 0-based index
    pairs are also 0-based index
    """
    ne, max_neb, elcon = FindNeigborConnection(enod)
    
    Nel  = len(ne)
    indi = list(elcon[:,0])
    indj = list(elcon[:,1])
    
    Npair = len(pairs)
    # loop through all the split surfaces
    # and append element connections into the indi and indj
    if Npair>0:
        for i in range(Npair):
            p1 = pairs[i,0]
            p2 = pairs[i,1]
            bnd1 = bnds[p1]
            bnd2 = bnds[p2]
            
            if bnd1.shape[1] != bnd2.shape[1]:
                print("bnd %d and %d have different number of elements, program stop!"%(p1+1, p2+1))
                exit()
            
            bnd1_node1 = enod[bnd1[1, :], bnd1[0, :]]
            bnd1_node2 = enod[(bnd1[1, :] + 1)%4, bnd1[0, :]]
            bnd2_node1 = enod[bnd2[1, :], bnd2[0, :]]
            bnd2_node2 = enod[(bnd2[1, :] + 1)%4, bnd2[0, :]]
            xy_bnd1_c = (coord[:, bnd1_node1] + coord[:, bnd1_node2])/2.0
            xy_bnd2_c = (coord[:, bnd2_node1] + coord[:, bnd2_node2])/2.0
            
            #sort bnd1 and bnd2 using the maximum range of center xy
            Lx   = np.max(xy_bnd1_c[0, :]) - np.min(xy_bnd1_c[0, :])
            Ly   = np.max(xy_bnd1_c[1, :]) - np.min(xy_bnd1_c[1, :])
            if Lx>Ly:
                ind = xy_bnd1_c[0, :].argsort()
            else:
                ind = xy_bnd1_c[1, :].argsort()
            
            xy_bnd1_c = xy_bnd1_c[:, ind]
            bnd1 = bnd1[:, ind]
            
            Lx   = np.max(xy_bnd2_c[0, :]) - np.min(xy_bnd2_c[0, :])
            Ly   = np.max(xy_bnd2_c[1, :]) - np.min(xy_bnd2_c[1, :])
            if Lx>Ly:
                ind = xy_bnd2_c[0, :].argsort()
            else:
                ind = xy_bnd2_c[1, :].argsort()
            
            xy_bnd2_c = xy_bnd2_c[:, ind]
            bnd2 = bnd2[:, ind]
            
            # compair if bnd1 and bnd2 have one-to-one connection
            tiny_xy = 1e-6
            if np.any(np.abs(xy_bnd1_c-xy_bnd2_c)>tiny_xy):
                print("bnd %d and %d do not match in position, program stop!"%(p1+1, p2+1))
            
            # append element index pairs into indi, indj
            for i in range(bnd1.shape[1]):
                indi.append(bnd1[0,i])
                indj.append(bnd2[0,i])
                indi.append(bnd2[0,i])
                indj.append(bnd1[0,i])
            
    v     =  np.ones(len(indi))
    graph = csr_matrix((v, (indi, indj)),shape=(Nel, Nel))
    r = reverse_cuthill_mckee(graph) # rcm permutation
    
    return r

def main():
    """
    run this function when module is executed as a script
    'python inp2mesh2d.py inpfile/geofile mesh2dfile'
    """
    
    # implement better ways to collect command line arguments
    parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
    description=("""\
     ********************** inp2mesh2d.py *****************************
     A python module that convert the Gmsh generated mesh to .mesh2d file  

     reverse_cuthill_mckee algorithm is used to re-order the element so that
     bandwidth of graph can be reduced.

     Since mesh2d needs boundary conditions, physical boundaries must be defined 
     in GMSH and the tag must be consistent with definition in the sem2dpack
     input file
     
     mesh file can be generated from geofile by
     gmsh geofile.geo -2 -o mesh -format xxformat -save_all
     
     Chao Liang 11/24/2021 @Chengdu, Sichuan
     
     Note:
     ** different boundaries must be defined in gmsh as physical curves 
        and labeled as bc1, bc2, bc3,...

        ** different material must be defined in gmsh as physical surface
           and labeled as mat1, mat2, ...

     This module uses meshio library.
     see https://github.com/nschloe/meshio, instruction for installation
     """))
    parser.add_argument("inp_geo_file",type=str,help="""Path where geo file or inp file is stored, when geofile is used, gmsh must be installed!""")
    parser.add_argument("mesh2dfile",type=str,help="Path where mesh2d file is stored")
    parser.add_argument("pairs", type=int, nargs='*',help='split surface pairs')

    args = parser.parse_args()
    inpfile = args.inp_geo_file
    mesh2dfile = args.mesh2dfile
    pairs = args.pairs
    if len(pairs)%2 != 0:
        print('pairs must be even number of integers!')
        exit()
    
    if len(pairs)>0:
        pairs = np.array(pairs, dtype=np.int)
        # reshape and recast to 0-based numbering 
        pairs = np.reshape(pairs, (2,-1)).transpose() - 1

    if inpfile.endswith('.geo'):
        # generate .inp file from .geo file
        geofile = inpfile
        inpfile = geofile[:-4] + '.inp'
        print(geofile, inpfile)
        geo2inp(geofile)

    print('Converting %s to %s'%(inpfile, mesh2dfile))
    # converting inp mesh file to mesh2d
    inp2mesh2d(inpfile, mesh2dfile, pairs)
    return 

if __name__=="__main__":
    main()

