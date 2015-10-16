### utils.py ---
##
## Filename: utils.py
## Author: Fred Qi
## Maintainer:
## Copyright (C) 2008, all rights reserved.
## Created: 2013-12-29 18:16:33(+0400)
##
## Last-Updated: 2014-04-23 14:11:21(+0400) [by Fred Qi]
##     Update #: 1211
######################################################################
##
### Commentary:
##    Commonly used routines.
##
######################################################################
##
### Change Log:
##
##
######################################################################

import yaml
from StringIO import StringIO
import numpy as np
import os.path as path
import subprocess as sp

import poscar as posrdr


def need_update(fn_dst, fn_srcs):
    """Check whether the destination file requires updating."""
    if not path.exists(fn_dst):
        return True

    tm_dst = path.getmtime(fn_dst)

    tm_src = -1
    for fn in fn_srcs:
        tm = path.getmtime(fn)
        tm_src = tm if tm > tm_src else tm_dst

    # print fn_dst + ': ', tm_src, tm_dst
    return tm_src > tm_dst


## ###########################################################################
## Util functions for data conversion
## ###########################################################################


def poscar_to_mat(task, force=False):
    # print 'Convering POSCAR to matlab matrix', task
    import scipy.io as sio
    update = need_update(task['output'], [task['input']])
    if force or update:
        data = posrdr.load_poscar(task['input'])
        sio.savemat(task['output'], {'poscar': data}, oned_as='column')


def poscar_to_fingerprint(task, force=False):
    """Compute the fingerprints of structures specified by POSCARS."""

    m_scr = """poscars = load( '%s' );
addpath('/home/fred/github/icrystal/scripts');
addpath('/home/fred/github/icrystal/scripts/USPEX');
FP = fingerprint_batch( poscars.poscar, %s );
save( '%s', 'FP' );
exit( );
"""
    # the last file in poscar_fns is a matlab mat file
    m_cmd = m_scr % (task['input'], task['params'], task['output'])
    # print m_cmd

    update = need_update(task['output'], [task['input']])
    if force or update:
        sp.call(['matlab', '-nodisplay', '-r "%s"' % m_cmd])


def individual_to_npz(task, force=False):
    """Convert individuals to npz."""
    update = need_update(task['output'], [task['input']])
    if force or update:
        comp, indiv_data = load_individuals_txt(task['input'])
        data = {'Composition': comp}
        cols = ('Enthalpy', 'Q_entropy', 'A_order', 'S_order')
        for ind in range(len(cols)):
            data[cols[ind]] = indiv_data[:, ind]

        np.savez(task['output'], **data)


def individual_to_mat(task, force=False):
    """Convert individuals to npz."""
    import scipy.io as sio
    update = need_update(task['output'], [task['input']])
    if force or update:
        comp, indiv_data = load_individuals_txt(task['input'])
        data = {'comp': comp,
                'enth': indiv_data[:, 0],
                'a_order_list': indiv_data[:, 2]}

        sio.savemat(task['output'], data)


def individual_filter(task, force=False):
    """Convert selected fields of individuals to npz."""
    update = need_update(task['output'], [task['input']])
    if force or update:
        comp, indiv_data = load_individuals_txt(task['input'])
        data = {'Composition': comp}
        cols = ('Enthalpy', 'Q_entropy', 'A_order', 'S_order')
        for ind in range(len(cols)):
            data[cols[ind]] = indiv_data[:, ind]

        sel = filter_individual_naive(indiv_data[:, 0], comp,
                                      indiv_data[:, -2])
        data['selected'] = sel

        np.savez(task['output'], **data)


def composition_filter(task, force=False):
    """Convert selected fields of individuals to npz."""
    update = need_update(task['output'], [task['input']])
    if force or update:
        data = posrdr.load_poscar(task['input'])
        comps = set([dt['num_ions'][0] for dt in data])
        comps = list(comps)
        sel = np.zeros((len(data),), dtype=np.int)
        for dt, idx in zip(data, range(len(data))):
            sel[idx] = comps.index(dt['num_ions'][0])

        kfolds = np.zeros((len(data), len(comps)), dtype=np.bool)
        for idx in range(len(comps)):
            kfolds[:, idx] = np.where(sel == idx, True, False)
            # print idx, np.sum(kfolds[:, idx])

        np.savetxt(task['output'], kfolds, fmt='%d')


def generation_to_kfold(task, force=False):
    """Extract generation info and create kfold from USPEX results."""
    update = need_update(task['output'], [task['input']])
    if force or update:
        create_generation_kfold(task['input'], task['output'])


def average_enthalpy(task, force=False):
    """Average the enthalpy by the number of atoms."""
    enth_outfn = task['output']
    enth_infn, poscar_infn = task['input'], task['auxinput']
    update = need_update(enth_outfn, (enth_infn, poscar_infn))
    if force or update:
        enth, gen = load_fitness(enth_infn)
        poscars = posrdr.load_poscar(poscar_infn)
        for poscar, idx in zip(poscars, range(len(enth))):
            enth[idx] /= poscar['num_ions'][0]
            save_fitness(enth_outfn, enth, gen)


def disp_performance(perfs, methods=('pearson', 'rmse', 'mae', 'median_ae')):
    """Display performance data."""
    flt_fmt = lambda val: '%.4g' % val
    percent_fmt = lambda val: '%.4g%%' % val

    disp_funcs = {'pearson': percent_fmt,
                  'rmse': flt_fmt,
                  'mae': flt_fmt,
                  'median_ae': flt_fmt}

    disp_str = {'pearson': 'Pearson Coef',
                'rmse': 'Root-Mean SE',
                'mae': 'Mean AE',
                'median_ae': 'Median AE'}

    if isinstance(methods, str):
        methods = (methods)

    msgs = None

    if hasattr(methods, '__iter__'):
        msgs = list()
        for meth in methods:
            fmt = disp_funcs[meth]
            avg, folds = perfs[meth]['overall'], perfs[meth]['folds']
            folds_str = ' '.join((fmt(v) for v in folds))
            msg = '%12s: %6s [%s]' % (disp_str[meth], fmt(avg), folds_str)
            msgs.append(msg)

    return msgs


def eval_performance(yo, yt, folds,
                     methods=('pearson', 'rmse', 'mae', 'median_ae')):
    """Evaluate performance with given methods."""
    perf_func = {'pearson': pearson_coef,
                 'rmse': rmse,
                 'mae': mae,
                 'median_ae': median_ae}

    if isinstance(methods, str):
        methods = (methods)

    perfs = None

    if hasattr(methods, '__iter__'):
        perfs = dict()
        nfolds = int(np.max(folds)+1)
        for meth in methods:
            efunc = perf_func[meth]
            avg = efunc(yo, yt)
            perf = np.zeros((nfolds,))
            for idx in range(nfolds):
                test = np.where(folds == float(idx), True, False)
                perf[idx] = efunc(yo[test], yt[test])
            perfs[meth] = {'overall': avg, 'folds': perf}

    return perfs


def pearson_coef(yo, yt):
    """Compute the Pearson coefficent."""

    myo = np.mean(yo)
    myt = np.mean(yt)

    zyo = yo - myo
    zyt = yt - myt

    pn = np.linalg.norm(zyo)*np.linalg.norm(zyt)

    return np.dot(zyo, zyt)*100.0 / (pn+1e-8)


def rmse(yo, yt):
    """Compute the root mean square error of the computation."""
    dts = np.square(yt - yo)
    return np.sqrt(np.mean(dts))


def mae(yo, yt):
    """Compute the mean absolute error of the computation."""
    return np.mean(np.abs(yt - yo))


def median_ae(yo, yt):
    """Compute the median absolute error of the computation."""
    return np.median(np.abs(yt - yo))


def run_matlab(tmplfn, param, force=False):
    """Run matlab command"""

    if not path.exists(tmplfn):
        return

    src_fns = [param['dist_fn'], tmplfn]
    update = need_update(param['out_fn'], src_fns)
    if not update and not force:
        return

    m_tmpl_file = open(tmplfn)
    tmpl = ''.join(m_tmpl_file.readlines())
    m_tmpl_file.close()

    mcmd = tmpl.format(**param)
    sp.call(['matlab', '-nodisplay', '-r "%s\n exit();"' % mcmd])


def gen_random_index(length, valth):
    """Generate random indexes for extracting training and testing samples.
    :Parameters:
       length : number of indicators to be generated.
       valth : the value used as threshold to determine training and testing
               sets.

    :Returns:

       idx_train : indicator vector of the training set.
       idx_test : indicator vector of the testing set.
       percentage : percentage of samples used for training.

    """

    idx_r = np.random.uniform(-1, 1, size=length)
    idx_train = idx_r >= valth
    idx_test = np.logical_not(idx_train)
    percentage = np.sum(idx_train)*100.0/length

    return idx_train, idx_test, percentage


def create_generation_kfold(prop_fn, kfold_fn):
    """Create a persistent k-fold from USPEX generated fitness data."""
    _, tests = load_fitness(prop_fn)
    n, n_folds = len(tests), int(np.max(tests))
    kfolds = np.zeros((n, n_folds), dtype=np.bool)
    for idx in range(n_folds):
        kfolds[:, idx] = np.where(tests > idx, True, False)

    np.savetxt(kfold_fn, kfolds, fmt='%d')


def convert_enthalpy(infn, outfn):
    """Convert the old version enthalpy data to the latest version."""
    import re

    if isinstance(infn, basestring):
        infile = open(infn)
    elif hasattr(infn, 'readlines'):
        infile = infn

    re_txt = '^\\\\n\\\\n(?P<comment> generation)*(?P<num>\d+)\\\\n.*$'
    re_line = re.compile(re_txt)

    nlines = []
    cmt_fmt = '-------- %s%s --------\n'
    for ln in infile.readlines():
        res = re_line.match(ln)
        if res:
            if res.group('comment'):
                nlines.append(cmt_fmt % res.groups())
        else:
            nlines.append(ln)

    if isinstance(infn, basestring):
        infile.close()

    if isinstance(outfn, basestring):
        outfile = open(outfn, 'w')
    elif hasattr(outfn, 'write'):
        outfile = outfn

    outfile.writelines(nlines)

    if isinstance(outfn, basestring):
        outfile.close()


def load_individuals_txt(fname, cols=('Composition',
                                      'Enthalpy', 'Q_entropy',
                                      'A_order', 'S_order')):
    """Load the individual file using split."""

    if isinstance(fname, basestring):
        infile = open(fname)
    elif hasattr(fname, 'readlines'):
        infile = fname

    lines = infile.readlines()

    ln = lines[1].split()
    comp_start, comp_end = ln.index('['), ln.index(']')
    comp_slice = range(comp_start+1, comp_end)
    idx_enthalpy = comp_end + 1

    nsample = len(lines) - 1
    data = np.zeros((nsample, 4))
    comp = np.zeros((nsample, len(comp_slice)))

    for idx in range(nsample):
        ln = lines[idx+1].split()
        comp[idx, :] = [int(dt) for dt in ln[comp_start+1:comp_end]]
        data[idx, 0] = float(ln[idx_enthalpy])
        data[idx, 1:] = [float(dt) for dt in ln[-3:]]

    if isinstance(fname, basestring):
        infile.close()

    return comp, data


def load_individuals_re(fname, cols=('Composition',
                                     'Enthalpy', 'Q_entropy',
                                     'A_order', 'S_order')):
    """Load the individual file using regular expression match."""

    import re

    if isinstance(fname, basestring):
        infile = open(fname)
    elif hasattr(fname, 'readlines'):
        infile = fname
    re_patts = ['\s*(?P<Gen>\d+)\s*(?P<ID>\d+)\s*(?P<Origin>\w+)\s*',
                '\[\s*(?P<Composition>(\d+\s*)+)\s*\]\s*',
                '(?P<Enthalpy>\d+\.\d+)\s*(?P<Volume>\d+\.\d+)\s*',
                '\[\s*(?P<KPOINTS>(\d+\s*)+)\s*\]\s*(?P<SYMM>\d+)\s*',
                '(?P<Q_entropy>\d+\.\d+)\s*(?P<A_order>\d+\.\d+)\s*',
                '(?P<S_order>\d+\.\d+).*']
    lnre = re.compile(''.join(re_patts), re.I)
    lines = infile.readlines()

    nsample = len(lines) - 1
    data = np.zeros((nsample, 4))
    if 'Composition' in cols:
        mres = lnre.match(lines[1])
        ncomp = len(mres.group('Composition').split())
        comp = np.zeros((nsample, ncomp), int)
        for idx in range(nsample):
            mres = lnre.match(lines[idx+1])
            comp_str = mres.group('Composition')
            comp[idx, :] = [int(dt) for dt in comp_str.split()]

    idx_enthalpy = lines[1].split().index(']')+1
    for idx in range(nsample):
        ln = lines[idx+1].split()
        data[idx, 0] = float(ln[idx_enthalpy])
        data[idx, 1:] = [float(dt) for dt in ln[-3:]]

    if isinstance(fname, basestring):
        infile.close()

    if 'Composition' in cols:
        return comp, data

    return data


def save_individuals(npz_fn, selected, comp, data,
                     cols=('Enthalpy', 'Q_entropy', 'A_order', 'S_order')):
    """Save individuals data into npz format."""
    np.savez(npz_fn, selected=selected, comp=comp,
             enthalpy=data[:, 0], Q_entropy=data[:, 1],
             A_order=data[:, 2], S_order=data[:, 3])


def same_composition(compa, compb, tol=0.0001):
    """To test whether the compositions of two structures are the same?
    % added: Dec 26, 2010; USPEX 8.4.2
    % checks whether compositions c1 and c2 are identical, tolerance = 0.0001
    % identity : c1 = n*c2 or c2 = n*c1 where n is integer
    """

    za = np.where(compa == 0, False, True)
    zb = np.where(compb == 0, False, True)

    if not np.all(za == zb):
        return False

    rt = compa[za]*1.0/compb[za]
    rmax, rmin = np.amax(rt), np.amin(rt)

    return rmax-rmin < tol


def filter_individual_naive(enth, comp, aorder):
    """Filtering same structures using a naive mechanism based on enthalpy,
composition, and a_order.

    """

## Perf
# ncalls  tottime  percall  cumtime  percall filename:lineno(function)
#      1    0.000    0.000    5.301    5.301 pycrystal/utils.py:104(individual_filter)
#      1    0.565    0.565    5.281    5.281 pycrystal/utils.py:319(filter_individual_naive)
#  99702    1.779    0.000    4.684    0.000 pycrystal/utils.py:300(same_composition)
# 299107    1.184    0.000    1.184    0.000 {method 'reduce' of 'numpy.ufunc' objects}
# 199404    1.091    0.000    1.091    0.000 {numpy.core.multiarray.where}
#  99702    0.092    0.000    0.767    0.000 numpy/core/fromnumeric.py:1643(all)
#  99702    0.075    0.000    0.563    0.000 numpy/core/fromnumeric.py:1847(amax)
#  99702    0.052    0.000    0.516    0.000 {method 'all' of 'numpy.ndarray' objects}

    nsample = len(enth)
    assert len(comp) == nsample
    assert len(aorder) == nsample

    accepted = np.ones((len(enth), )) == 1
    scomp = np.sum(comp, 1)
    for idx in range(nsample - 1):
        for jj in range(idx+1, nsample):
            if not accepted[jj]:
                continue
            if same_composition(comp[idx, :], comp[jj, :]):
                order_diff = abs(aorder[idx] - aorder[jj])
                enth_diff = abs(enth[idx]/scomp[idx] - enth[jj]/scomp[jj])
                if enth_diff < 0.004 and order_diff < 0.15:
                    accepted[jj] = False

    return accepted


def load_fitness(filename):
    """Load the fitness data from a given txt file."""
    import numpy as np

    if isinstance(filename, str):
        fit_file = open(filename)
        lines = np.array(fit_file.readlines())
        fit_file.close()
    elif hasattr(filename, 'readlines'):
        lines = np.array(filename.readlines())
    else:
        return None, None

    gen_sep = [True if '--' == ln[:2] else False for ln in lines]
    fit_idx = np.logical_not(gen_sep)

    data = np.array([float(ln.split()[-1]) for ln in lines[fit_idx]])
    gen_idx = np.cumsum(gen_sep) - 1

    return data, gen_idx[fit_idx]


def save_fitness(filename, fitness, gen=None, with_index=True):
    """Save fitness with generation and index information."""
    fit_file = None
    if isinstance(filename, str):
        fit_file = open(filename, 'w')
    elif hasattr(filename, 'write'):
        fit_file = filename
    if fit_file is None:
        return False

    if with_index:
        data = zip(range(1, len(fitness)+1), fitness)
        lines = ['%d %.7g\n' % dt for dt in data]
    else:
        lines = ['%.7g\n' % fit for fit in fitness]

    if gen is None:
        fit_file.writelines(lines)
    else:
        line_comment = '------- generation%d -------\n'
        cg = -1
        for g, ln in zip(gen, lines):
            if g != cg:
                fit_file.write(line_comment % (g+1))
                cg = g
            fit_file.write(ln)

    if isinstance(filename, str):
        fit_file.close()


def load_fingerprint(fing_fn):
    """Load fingerprint from matlab data file."""
    import scipy.io as sio
    data = sio.loadmat(fing_fn, struct_as_record=False, squeeze_me=True)
    
    if 1 == len(data['FP'][0].FingerPrint.shape):
        (ncols,) = data['FP'][0].FingerPrint.shape
        data_len = len(data['FP'])
        fing = np.zeros((data_len, ncols))
        for idx in range(data_len):
            fing[idx, :] = data['FP'][idx].FingerPrint.reshape(1, ncols)
    else:
        (srow, ncol) = data['FP'][0].FingerPrint.shape
        nrow = int(np.sqrt(srow))
        data_len, fing_len = len(data['FP']), (srow+nrow)/2*ncol
        fp_sel = np.zeros((srow, )) == 1.0
        for idx in range(srow):
            row, col = idx/nrow, idx % nrow
            if row <= col:
                fp_sel[idx] = True
        # print fp_sel, len(fp_sel)*ncol, np.sum(fp_sel)*ncol, fing_len
        fing = np.zeros((data_len, fing_len))
        for idx in range(data_len):
            rfp = data['FP'][idx].FingerPrint[fp_sel, :]
            fing[idx, :] = rfp.reshape(1, fing_len)

    return fing


def load_material_data(fing_fn=None, prop_fn=None):
    """Load the data of fingerprint and corresponding fitness from disk.

    :Parameters:
       prop_fn : the name of the file storing properties.
       fing_fn : the name of the file storing fingerprints.

    :Returns:
       prop : property vector.
       fing : fingerprint matrix.
       gen_idx : generation index.
    """
    prop, fing, gen_idx = None, None, None

    if None != prop_fn:
        prop, gen_idx = load_fitness(prop_fn)

    if None != fing_fn:
        fing = load_fingerprint(fing_fn)

    # assume the data is cut at tail
    if (not prop_fn is None) and (not fing_fn is None):
        data_len = min(len(prop), len(fing))
        return prop[:data_len], fing[:data_len, :], gen_idx
    elif fing_fn is None:
        return prop, gen_idx


def load_individuals_npz(indiv_fn):
    """Load filtered individuals data from a npz format file."""

    data = np.load(indiv_fn)

    enthalpy = data['enthalpy']

    obs_t = np.vstack((data['Q_entropy'],
                       data['A_order'],
                       data['S_order']))

    return enthalpy, np.transpose(obs_t)


def savetxt_ext(filename, X, fmt='%.18e', delim=',', header=None):
    """Save data into a file with a yaml header.

    """
    last_line = ''
    if header and 'cols' in header:
        last_line = '\ncols: ' + header['cols']
        del header['cols']
    if len(header):
        hdstr = yaml.safe_dump(header, default_flow_style=False)
    else:
        hdstr = '\n'
    np.savetxt(filename, X, fmt=fmt, delimiter=delim,
               header=hdstr+last_line)


def loadtxt_ext(fname):
    """Load data from a file with a yaml header.

    """
    hdlines = None
    if isinstance(fname, basestring):
        infile = open(fname)
    elif hasattr(fname, 'readlines'):
        infile = fname

    lines = infile.readlines()
    hdlines = [ln[2:] for ln in lines if '#' == ln[0]]
    str = StringIO(''.join(hdlines))
    header = yaml.load(str)

    infile.seek(0)
    data = np.loadtxt(infile, delimiter=',')

    if isinstance(fname, basestring):
        infile.close()

    return header, data


######################################################################
### utils.py ends here
