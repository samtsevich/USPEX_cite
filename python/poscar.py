## ----------------------------------------------------------------------
## START OF FILE
## ----------------------------------------------------------------------
##
## Filename: poscar.py
## Author: Fred Qi
## Created: 2013-08-27 15:29:42(+0800)
##
## ----------------------------------------------------------------------
### CHANGE LOG
## ----------------------------------------------------------------------
## Last-Updated: 2014-09-22 15:44:08(+0400) [by Fred Qi]
##     Update #: 562
## ----------------------------------------------------------------------


def load_poscar(filename):
    """
    Load atom coordinates from a USPEX generated POSCAR file,
    which contains serveral groups of data.
    Reference URL: http://cms.mpi.univie.ac.at/vasp/guide/node59.html
    """

    poscar_file = None
    if isinstance(filename, str):
        poscar_file = open(filename)
    elif hasattr(filename, 'readlines'):
        poscar_file = filename

    if poscar_file is None:
        return None

    data = []
    cnt, cnt_end = (0, 6)
    scaling = 1.0
    seldynamic = False
    direct = True
    for line in poscar_file.readlines():
        # Omit empty lines
        if 0 == len(line.strip()):
            continue
        # 'EA' == line[0:2]:
        if 0 == cnt:
            # The beginning of a POSCAR section, comment line
            atoms = dict()
            atoms['lattice'] = []
        elif 1 == cnt:
            scaling = float(line.strip())
        elif 1 < cnt and 5 > cnt:
            # Record the lattice from L3 ~ L5
            atoms['lattice'].append([float(x) for x in line.split()])
        elif 5 == cnt:
            items = line.split()
            try:
                atoms['num_ions'] = [int(x) for x in items]
                cnt_end = 6 + sum(atoms['num_ions'])
                atoms['coords'] = []
            except ValueError:
                atoms['atom_types'] = items
                continue
        elif 6 == cnt:
            mode = line.strip()[0].upper()
            if 'S' == mode:
                seldynamic = True
                continue
            elif 'C' == mode or 'K' == mode:
                direct = False
            elif 'D' == mode:
                direct = True
        elif 6 < cnt:
            atoms['coords'].append([float(x) for x in line.split()])
            if cnt_end == cnt:
                assert(len(atoms['coords']) == sum(atoms['num_ions']))
                # n = sum(atoms['num_ions'])
                poscar = {'scaling': scaling,
                          'direct': direct,
                          'seldynamic': seldynamic,
                          'num_ions': atoms['num_ions'],
                          'lattice': atoms['lattice'],
                          'positions': atoms['coords']}
                if 'atom_types' in atoms:
                    poscar['atom_types'] = atoms['atom_types']

                data.append(poscar)
                # End of one POSCAR section, reset state counter
                cnt = -1

        cnt += 1

    if isinstance(filename, str):
        poscar_file.close()
    return data


def save_poscar(filename, poscars):
    """Save POSCARS to file."""

    def write_matrix_by_line(thefile, mat):
        for dt in mat:
            ln = ' '.join(('%.11g' % d for d in dt))
            thefile.write(ln + '\n')

    if isinstance(filename, str):
        output = open(filename, 'w')
    elif hasattr(filename, 'write'):
        output = filename

    idx = 1
    for poscar in poscars:
        output.write('COMPOUND%d\n' % idx)
        output.write('%.1f\n' % poscar['scaling'])
        write_matrix_by_line(output, poscar['lattice'])
        ln = ('%d' % n for n in poscar['num_ions'])
        output.write(' '.join(ln) + '\n')
        output.write('Direct\n')
        write_matrix_by_line(output, poscar['positions'])

        idx += 1

    if isinstance(filename, str):
        output.close()


def _get_comp_id(nions, nion=256):
    return nions[0]*nion + nions[1]


def poscar_cluster(poscars):
    """Cluster structures according to their compounds."""
    comps = dict()
    for poscar, idx in zip(poscars, range(len(poscars))):
        comp_id = _get_comp_id(poscar['num_ions'])
        if not comp_id in comps:
            comps[comp_id] = list()
        comps[comp_id].append(idx)

    clusters = [0]*len(poscars)
    for index, idx in zip(comps.itervalues(), range(len(comps))):
        for ii in index:
            clusters[ii] = idx

    return clusters


def poscar_to_mat(poscar_fn, mat_fn, force=False):
    """[Obsoleted!]
    Convert the given POSCAR data to a matlab mat file with matfn.

    """

    import scipy.io as sio
    import utils as cu

    update = cu.need_update(mat_fn, [poscar_fn])
    if force or update:
        data = load_poscar(poscar_fn)
        sio.savemat(mat_fn, {'poscar': data}, oned_as='column')


def fitness_to_mat(fitness_fn, mat_fn):
    """Convert the fitness data to matlab mat file."""
    import scipy.io as sio
    data, gen_idx = load_fitness(fitness_fn)
    sio.savemat(mat_fn, {'fitness': data})


def display_poscar(data):
    """Display a section of POSCAR data."""

    def display_matrix(matrix):
        print '[', matrix[0]
        for line in matrix[1:-1]:
            print ' ', line
        print ' ', matrix[-1], ']'

    keys = ['scaling', 'atom_types', 'num_ions', 'seldynamic', 'direct',
            'lattice', 'positions']

    for k in keys:
        if not k in data:
            continue
        if 'positions' == k or 'lattice' == k:
            print k
            display_matrix(data[k])
        else:
            print k, (10-len(k))*' ', data[k]


def data_conversion():
    """Example usage"""
    data = load_poscar('../data/HfO2dataset/gatheredPOSCARS')
    poscar_to_mat(data, 'HfO2_poscar.mat')
    fitness_to_mat('../data/HfO2dataset/fitness_nospace.dat',
                   'HfO2_fit.mat')
    display_poscar(data[0])

    data = load_poscar('../data/HfSiOvarComp/gatheredPOSCARS')
    poscar_to_mat(data, 'HfSiO_poscar.mat')
    fitness_to_mat('../data/HfSiOvarComp/fitness_nospace1.dat',
                   'HSiO_fit.mat')


if __name__ == "__main__":

    # data_conversion( )
    import sys
    # print sys.argv

    # data, gen_idx = load_fitness( sys.argv[1] )
    # print data, len(data)

    # data = load_poscar(sys.argv[1])
    # display_poscar(data[0])

    if len(sys.argv) > 2:
        poscar_to_mat(sys.argv[1], sys.argv[2], True)

## ----------------------------------------------------------------------
### END OF FILE
## ----------------------------------------------------------------------
