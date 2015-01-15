#!/usr/bin/python
#-----------------------------------------------------------------------
# File    : chkfim.py
# Contents: examples how to use the functions of the fim module
# Author  : Christian Borgelt
# History : 2012.??.?? file created
#           2013.10.31 added carpenter and closed item set mining
#-----------------------------------------------------------------------
from sys        import argv, stderr
from os         import remove, devnull
from subprocess import call, check_output
from time       import time
from random     import seed, random
from fim        import fim, apriori, eclat, fpgrowth, sam, relim
from fim        import carpenter, ista

#-----------------------------------------------------------------------

def run_pipe (prog, null):
    res = []                    # traverse the program output
    for line in check_output(prog +['-'],stderr=null).split('\n')[:-1]:
        r = [int(x) for x in line.split()]
        res.append((frozenset(r[:-1]), tuple(r[-1:])))
    return res

#-----------------------------------------------------------------------

def run_file (prog, null):
    call(prog +['tmp.txt'], stderr=null)
    res = []
    with open('tmp.txt', 'rt') as inp:
        for line in inp:
            r = [int(x) for x in line.split()]
            res.append((frozenset(r[:-1]), tuple(r[-1:])))
    return res

#-----------------------------------------------------------------------

if __name__ == '__main__':
    runs   = int(argv[1]) if len(argv) > 1 else 1
    tracts = [[i+1 for i in range(100) if random() < 0.1]
                   for k in range(1000)]
    with open('data.txt', 'w') as out:
        for t in tracts:
            for i in t: out.write('%d ' % i)
            out.write('\n')

    stderr.write('frequent item sets:\n')
    stderr.write('fim    ... '); t = time()
    for r in range(runs):
        pypats = fim(tracts, supp=-2, zmin=2, report='a')
    stderr.write('done [%.3fs].\n' % (time()-t))
    ref = set([(tuple(sorted(list(s))), x[0]) for s,x in pypats])
    stderr.write('\n')

    for p,f in [('apriori',  apriori),
                ('eclat',    eclat),
                ('fpgrowth', fpgrowth),
                ('sam',      sam),
                ('relim',    relim)]:
        stderr.write(p +' (all)\n')
        stderr.write('python ... '); t = time()
        for r in range(runs):
            pypats = f(tracts, supp=-2, zmin=2, report='a')
        stderr.write('done [%.3fs].\n' % (time()-t))

        prg = [p, '-s-2', '-m2', '-v %a', 'data.txt']
        stderr.write('pipe   ... '); t = time()
        with open(devnull, 'w') as null:
            for r in range(runs):
                expats = run_pipe(prg, null)
        stderr.write('done [%.3fs].\n' % (time()-t))

        stderr.write('file   ... '); t = time()
        with open(devnull, 'w') as null:
            for r in range(runs):
                expats = run_file(prg, null)
        stderr.write('done [%.3fs].\n' % (time()-t))

        pypats = set([(tuple(sorted(list(s))), x[0]) for s,x in pypats])
        expats = set([(tuple(sorted(list(s))), x[0]) for s,x in expats])
        stderr.write('cmp ok\n' if pypats == expats else 'cmp fail\n')
        stderr.write('ref ok\n' if pypats == ref    else 'ref fail\n')
        stderr.write('\n')

    stderr.write('closed item sets:\n')
    stderr.write('fim    ... '); t = time()
    for r in range(runs):
        pypats = fim(tracts, target='c', supp=-2, zmin=2, report='a')
    stderr.write('done [%.3fs].\n' % (time()-t))
    ref = set([(tuple(sorted(list(s))), x[0]) for s,x in pypats])
    stderr.write('\n')

    for p,f in [('apriori',   apriori),
                ('eclat',     eclat),
                ('fpgrowth',  fpgrowth),
                ('sam',       sam),
                ('relim',     relim),
                ('carpenter', carpenter),
                ('ista',      ista)]:
        stderr.write(p +' (closed)\n')
        stderr.write('python ... '); t = time()

        for r in range(runs):
            pypats = f(tracts, target='c', supp=-2, zmin=2, report='a')
        stderr.write('done [%.3fs].\n' % (time()-t))

        prg = [p, '-tc', '-s-2', '-m2', '-v %a', 'data.txt']
        stderr.write('pipe   ... '); t = time()
        with open(devnull, 'w') as null:
            for r in range(runs):
                expats = run_pipe(prg, null)
        stderr.write('done [%.3fs].\n' % (time()-t))

        stderr.write('file   ... '); t = time()
        with open(devnull, 'w') as null:
            for r in range(runs):
                expats = run_file(prg, null)
        stderr.write('done [%.3fs].\n' % (time()-t))

        pypats = set([(tuple(sorted(list(s))), x[0]) for s,x in pypats])
        expats = set([(tuple(sorted(list(s))), x[0]) for s,x in expats])
        stderr.write('cmp ok\n' if pypats == expats else 'cmp fail\n')
        stderr.write('ref ok\n' if pypats == ref    else 'ref fail\n')
        stderr.write('\n')

    stderr.write('maximal item sets:\n')
    stderr.write('fim    ... '); t = time()
    for r in range(runs):
        pypats = fim(tracts, target='m', supp=-2, zmin=2, report='a')
    stderr.write('done [%.3fs].\n' % (time()-t))
    ref = set([(tuple(sorted(list(s))), x[0])
               for s,x in pypats])
    stderr.write('\n')

    for p,f in [('apriori',   apriori),
                ('eclat',     eclat),
                ('fpgrowth',  fpgrowth),
                ('sam',       sam),
                ('relim',     relim),
                ('carpenter', carpenter),
                ('ista',      ista)]:
        stderr.write(p +' (maximal)\n')
        stderr.write('python ... '); t = time()
        for r in range(runs):
            pypats = f(tracts, target='m', supp=-2, zmin=2, report='a')
        stderr.write('done [%.3fs].\n' % (time()-t))

        prg = [p, '-tm', '-s-2', '-m2', '-v %a', 'data.txt']
        stderr.write('pipe   ... '); t = time()
        with open(devnull, 'w') as null:
            for r in range(runs):
                expats = run_pipe(prg, null)
        stderr.write('done [%.3fs].\n' % (time()-t))

        stderr.write('file   ... '); t = time()
        with open(devnull, 'w') as null:
            for r in range(runs):
                expats = run_file(prg, null)
        stderr.write('done [%.3fs].\n' % (time()-t))

        pypats = set([(tuple(sorted(list(s))), x[0])
                      for s,x in pypats])
        expats = set([(tuple(sorted(list(s))), x[0])
                      for s,x in expats])
        stderr.write('cmp ok\n' if pypats == expats else 'cmp fail\n')
        stderr.write('ref ok\n' if pypats == ref    else 'ref fail\n')
        stderr.write('\n')

    remove('data.txt')
    remove('tmp.txt')
