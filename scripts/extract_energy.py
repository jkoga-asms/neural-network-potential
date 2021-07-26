#!/usr/bin/env python

import optparse
import pathlib
import math

def process_structure(strstruct, atom_energies):
    in_structure = False
    refe = 0.0
    natm = 0
    for line in strstruct:
        if line.startswith('-------------'):
            in_structure = True
            continue
        if line.startswith('Cohesive energy'):
            in_structure = False
            continue
        if in_structure:
            natm += 1
            if atom_energies is not None:
                words = line.split()
                refe += atom_energies[words[0]]
        if line.startswith('Total energy'):
            words = line.split()
            energy = (float(words[3]) - refe)
            return natm, energy

def parse_predict_out(pout, atom_energies):
    f = open(pout)
    energy_eval = False
    strstruct = None
    next_line_is_filename = False
    dirs = []
    ces = []
    for line in f:
        line = line.strip()
        if len(line)==0:
            continue
        if line.startswith('Energy evaluation'):
            energy_eval = True
            continue
        if not energy_eval:
            continue
        if line.startswith('Structure number'):
            strstruct = []
        if next_line_is_filename:
            filename = line
            next_line_is_filename = False 
            continue
        if line.startswith('File name'):
            words = line.split()
            if len(words)>=4:
                filename = words[3]
            else:
                next_line_is_filename = True
                continue
        if line.startswith('All forces are'):
            natm, energy = process_structure(strstruct, atom_energies)
            #print('file '+filename+' natm '+str(natm)+' cohesive energy per atm '+str(energy/float(natm)))
            dirs.append(filename)
            ces.append(energy/float(natm)) 
        if strstruct is not None:
            strstruct.append(line)
    f.close()
    return dirs, ces

def get_energy_from_xsf(xsf, atom_energies):
    f = open(xsf)
    tote = 0.0
    refe = 0.0
    natm = 1
    primcoord_read = False
    natm_read = False
    for line in f:
        line = line.strip()
        if line.startswith('# total energy'):
            words = line.split()
            tote = float(words[4])
            continue
        if line.startswith('PRIMCOORD'):
            primcoord_read = True
            continue
        if not primcoord_read:
            continue
        if not natm_read:
            words = line.split() 
            natm = int(words[0])
            natm_read = True
            continue
        if natm_read and atom_energies is not None:
            words = line.split()
            refe += atom_energies[words[0]]
    return (tote-refe)/float(natm)

def run():
    parser = optparse.OptionParser(usage='%prog [options]',description='predict.xのログファイルからエネルギーを抽出し，対応するXSFのエネルギーとならべて出力する')
    parser.add_option('-n', '--nnp', dest='nnp', type=str, default='predict.out', help='predict.xのログを記録したファイル。デフォルト値はpredict.out')
    parser.add_option('-a', '--atomenergy', dest='atomenergy', type=str, default=None, help=\
    '参照エネルギーを計算するための原子のエネルギーをあたえる。element1:energy1,element2:energy2,...という形式で与える。デフォルト値はNoneで，この場合生のエネルギーが出力される')
    parser.add_option('-r', '--result', dest='result', type=str, default='energy.dat', help=\
   '結果を出力するファイル。1行がある原子配置で，各行の1カラム目がNNPのエネルギー，2カラム目が第一原理計算のエネルギーに対応する。デフォルト値はenergy.dat')
    (options,args) = parser.parse_args()
    atom_energies = None
    if options.atomenergy is not None:
        atom_energies = {}
        aes = options.atomenergy.split(',')
        for ae in aes:
            aa = ae.split(':')
            atom_energies[aa[0]] = float(aa[1])

    dirs, nnpes = parse_predict_out(options.nnp, atom_energies)
    f = open(options.result,'w') 
    maxe = -1.e30
    maxediff = -1.e30
    avediff = 0.0
    rmse = 0.0
    maxedir = None
    idir = 0
    for i in range(len(dirs)):
        ennp = nnpes[i]
        eqe = get_energy_from_xsf(dirs[i], atom_energies)
        f.write(str(ennp)+' '+str(eqe)+'\n')
        if math.fabs(ennp-eqe)>maxediff:
            maxediff = math.fabs(ennp-eqe)
            maxedir = dirs[i]
        avediff += math.fabs(ennp-eqe)
        rmse += math.fabs(ennp-eqe) * math.fabs(ennp-eqe)
        idir += 1
        #if maxe < ennp:
        #    maxe = ennp
        #    maxedir = dirs[i]
    f.close()
    if maxedir is not None:
        print('max diff : '+str(maxediff)+' eV/atom, and the corresponding directory : '+maxedir)
        print('average diff : '+str(avediff/float(idir))+' eV/atom')
        print('rmse :  '+str(math.sqrt(rmse/float(idir)))+' eV/atom')

if __name__ == '__main__':
    run()

