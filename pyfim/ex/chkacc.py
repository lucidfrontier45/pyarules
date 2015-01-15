#!/usr/bin/python

from sys        import argv, stderr
from os         import remove, devnull
from subprocess import call, check_output
from time       import time
from random     import seed, random
from fim        import apriori, apriacc, accretion

#-----------------------------------------------------------------------

def run_pipe (prog, null):
    res = []                    # traverse the program output
    for line in check_output(prog+['-'], stderr=null).split('\n')[:-1]:
        r = [int(x) for x in line.split()]
        res.append((frozenset(r[:-1]), tuple(r[-1:])))
    return res

#-----------------------------------------------------------------------

def run_file (prog, null):
    call(prog+['tmp.txt'], stderr=null)
    res = []
    with open('tmp.txt', 'rt') as inp:
        for line in inp:
            r = [int(x) for x in line.split()]
            res.append((frozenset(r[:-1]), tuple(r[-1:])))
    return res

#-----------------------------------------------------------------------

if __name__ == '__main__':
    runs   = 1
    stat   = 'p'
    prune  = -2
    opts   = {'accretion': ['-e'+stat, '-d1'],
              'apriacc':   ['-e'+stat, '-d1', '-p'+str(prune)],
              'apriori':   ['-e'+stat, '-d1', '-tm', '-an', '-z',
                            '-p'+str(prune)] }
    args   = {'accretion': {'stat': stat, 'siglvl': 1},
              'apriacc':   {'stat': stat, 'siglvl': 1, 'prune': prune},
              'apriori':   {'eval': stat, 'thresh': 1, 'target': 'm',
                            'agg': 'n', 'mode': 'z', 'prune': prune}}

    #gen = ['genpst', '-n100', '-t1000', '-p0.1', 'data.txt']
    #with open(devnull, 'w') as null:
    #    call(gen, stderr=null)
    #with open('data.txt', 'r') as inp:
    #    tracts = [[int(x) for x in line.split()] for line in inp]
    
    if len(argv) > 1: seed(int(argv[1]))
    tracts = [[i+1 for i in range(100) if random() < 0.1]
                   for k in range(1000)]
    with open('data.txt', 'w') as out:
        for t in tracts:
            for i in t: out.write('%d ' % i)
            out.write('\n')

    for p,f in [('accretion', accretion),
                ('apriacc',   apriacc),
                ('apriori',   apriori)]:
        print p
        stderr.write('python ... '); t = time()
        for r in range(runs):
            pypats = f(tracts, supp=-2, zmin=2, **args[p])
        stderr.write('done [%.3fs].\n' % (time()-t))

        prg = [p, '-s-2', '-m2', '-v %a', 'data.txt'] +opts[p]
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

        pypats = [(tuple(sorted(list(s))), x[0]) for s,x in pypats]
        expats = [(tuple(sorted(list(s))), x[0]) for s,x in expats]
        print('ok' if set(pypats) == set(expats) else 'fail')

    remove('data.txt')
    remove('tmp.txt')
