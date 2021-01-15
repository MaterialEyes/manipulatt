# manipulatt - v2021.01.14
#Manipulate unitcell lattice by rotating, tilting and tiling 
import os
import numpy as np
import itertools as it
try:
    import pandas as pd
except ImportError:
    print('Please install pandas before proceeding')
    
try:
    from tqdm.notebook import tqdm
except ImportError:
    print('Please install tqdm before proceeding')
    
import matplotlib.pyplot as plt
from .util import trim_filename, _in_hull, centroid_coordinates

# pymatgen tools
try:
    from pymatgen import MPRester
except ImportError:
    print('Please install pymatgen before proceeding')
from pymatgen.io.xyz import XYZ
from pymatgen.io.ase import AseAtomsAdaptor
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import RotationTransformation

try:
    from matminer.featurizers.site import CrystalNNFingerprint
except ImportError:
    print('Please install matminer before proceeding')
from matminer.featurizers.structure import SiteStatsFingerprint


class supercell(object):
    def __init__(self, MaPKey = 'RyfhcbEkpUtKPxC6'): 
        self.MaPKey = MaPKey
    
    def find_unitcell_mat_proj(self, formula, space_grp=''):
        
        '''
        Generates pymatgen unit cell for the specified chemical formula and space group
    
        Accepts:
                formula:             string; chemical formula to query in materials project (MP) repository
                space_grp:           int (optional) chemical space group number to query MP repository
            
        Returns:
                pymatgen_structure:  pymatgen structure of conventional standard unit cell of specified formula
                                     and space group number

        '''
    
        mpr = MPRester(self.MaPKey)
        query = mpr.query(criteria={"pretty_formula": formula}, properties=["structure","icsd_ids","spacegroup"])
        if space_grp:
            query = [query[i] for i in range(len(query)) if SpacegroupAnalyzer(query[i]['structure']).get_space_group_number() == space_grp]
        selected = query[np.argmin([query[i]['structure'].lattice.volume for i in range(len(query))])]
        pymatgen_structure = SpacegroupAnalyzer(selected["structure"]).get_conventional_standard_structure()
    
        return pymatgen_structure
    
    def construct_supercell_from_unitcell(self, unit_cell, max_dim):
    
        '''
        Constructs supercell of size max_dim[0] x max_dim[1] x max_dim[2] angstrom^3 for any given unit cell
    
        Accepts:
                unit_cell:         pymatgen unit cell structure, can be any rotated and tilted unit cell
                max_dim:           list or tuple of size 3 with the cell dimension of the simulation box
            
        Returns:
                supercell_cropped: numpy array of shape (number of atoms, 4), first 3 columns are atomic coordinates,
                                   4th column is the atomic number of each element present

        '''
    
        assert isinstance(max_dim, (list, tuple)) , "Error: max_dim must be a list or tuple with 3 elements"
        assert(len(max_dim) == 3) , "Error: Number of elemets must be 3"
    
        lattice_matrix = unit_cell.lattice.matrix
        bounding_box = np.matmul([[0,0,0],[1,0,0],[0,1,0],[0,0,1],[1,1,0],[1,0,1],[0,1,1],[1,1,1]],np.diag(max_dim))
        bound_tile = []

        for i in bounding_box:
            bound_tile.append(np.round(np.linalg.solve(lattice_matrix,i)).astype(np.int))
        
        bound_tile = np.asarray(bound_tile)
        bound_tile_min = np.min(bound_tile,axis=0)
        bound_tile_max = np.max(bound_tile,axis=0)
    
        supercell_tile = []
        unit_shape = unit_cell.frac_coords.shape[0]
        atomic_numbers = np.array(unit_cell.atomic_numbers).reshape(unit_shape,1)

        for i in range(bound_tile_min[2]-1,bound_tile_max[2]+1):
            for j in range(bound_tile_min[1]-1,bound_tile_max[1]+1):
                for k in range(bound_tile_min[0]-1,bound_tile_max[0]+1):
                    supercell_tile.append(np.concatenate((np.array(unit_cell.frac_coords +[k,j,i]).reshape(unit_shape,3),atomic_numbers),axis=1))
                    
        supercell_tile = np.array(supercell_tile)
        supercell_tile = supercell_tile.reshape(supercell_tile.shape[0]*supercell_tile.shape[1],supercell_tile.shape[2])
    
        supercell_tile[:,:3] = np.matmul(supercell_tile[:,:3],lattice_matrix.T)
        supercell_tile[:,:3] = supercell_tile[:,:3] - centroid_coordinates(supercell_tile[:,:3])
    
        dim_low = np.array([-max_dim[0]/2,-max_dim[1]/2,-max_dim[2]/2])
        dim_high = np.array([max_dim[0]/2,max_dim[1]/2,max_dim[2]/2])
        inidx = np.all(np.logical_and(supercell_tile[:,:3]>=dim_low,supercell_tile[:,:3]<=dim_high),axis=1)
    
        supercell_cropped = np.delete(supercell_tile,np.where(np.logical_not(inidx))[0],axis=0)    
        supercell_cropped[:,:3] = supercell_cropped[:,:3] - np.min(supercell_cropped[:,:3],axis=0)
    
        print('Min coordinates of the rotated supercell: \n')
        print(np.min(supercell_cropped,axis=0))
        print('Max coordinates of the rotated supercell: \n')
        print(np.max(supercell_cropped,axis=0))
        print('\n')
    
        return supercell_cropped

    def rotate_tilt_unitcell(self, unit_cell, zxzrotation):
    
        '''
        Rotate/tilt unit cell to a specified zxz rotation angles
    
        Accepts:
                unit_cell:         pymatgen strucure of unrotated conventional unit cell 
                zxzrotation:       list or tuple of size 3 with the cell dimension of the simulation box
            
        Returns:
                zxzrotated_obj:    pymagen structure of rotated and tilted standard unit cell

        '''
    
        #assert isinstance(zxzrotation, (list, tuple)) , "Error: zxzrotation must be a list or tuple with 3 elements"
        assert(len(zxzrotation) == 3) , "Error: Number of elemets must be 3"
    
        unit_cell = RotationTransformation([0,0,1],zxzrotation[0]).apply_transformation(unit_cell)
        unit_cell = RotationTransformation([1,0,0],zxzrotation[1]).apply_transformation(unit_cell)
        unit_cell = RotationTransformation([0,0,1],zxzrotation[2]).apply_transformation(unit_cell)
    
        return unit_cell

    def zoneaxis_rotate_tilt_unitcell(self, unit_cell, uvw_project):
        
        '''
        Construct a slab, oriented for a specific direction of viewing, and slice into a cube
        
        Accepts:
                pymatgen_structure:   The input pymatgen structure
                uvw_project:          The direction of projection (sceen to viewer) is a lattice vector ua + vb + wc.
                uvw_upward:           The upward direction is a lattice vector ua + vb + wc (must be normal to along_uvw).
                tilt_ang:             The CCW rotation around 'uvw_project' (applied after uvw_upward is set)
                max_dimension:        A float representing the max edge length of the supercell
            
        Returns:
                                      A pymatgen supercell (cube) oriented for a specific direction of viewing 
        '''

        atomObj = AseAtomsAdaptor.get_atoms(unit_cell)

        # Make copy of atom object and get cell/projection vector info
        atom  = atomObj.copy()

        # Convert u,v,w vector to cartesian
        along_xyz = np.matmul(np.array(atom.get_cell()).T, uvw_project)

        # Rotate coordinates and cell so that 'along_xyz' is coincident with [0,0,1] 
        atom.rotate(along_xyz,[0,0,1], rotate_cell=True)

        return AseAtomsAdaptor.get_structure(atom)

