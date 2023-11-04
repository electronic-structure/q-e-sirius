import sys
import re

def get_num_atoms(raw_str):
    result = re.search(r'\s+number\sof\satoms\/cell\s+\=\s+(?P<natoms>\S+)', raw_str)
    return int(result.group('natoms'))

def get_stress(raw_str):
    # find last occurance of the pattern
    for result in re.finditer(r'total\s+stress.+\(kbar\)\s+P\=.+\s+.+\s+.+\s+.+', raw_str):
        pass

    lines = result.group(0).splitlines()

    stress = []
    for i in [0, 1, 2]:
        stress.append([])
        vals = lines[i + 1].split()
        for j in [0, 1, 2]:
            stress[i].append(float(vals[j]))

    return stress

def get_pressure(raw_str):
    # find the last occurance of the pattern
    for result in re.finditer(r'\s+total\s+stress.+\(kbar\)\s+P\=\s+(?P<pressure>\S+)', raw_str):
        pass
    return(float(result['pressure']))

def get_forces(raw_str, natoms):
    fmt = r'     Forces acting on atoms \(cartesian axes, Ry\/au\):\s+'
    for i in range(natoms):
        fmt += r'.+\s+'

    for result in re.finditer(fmt, raw_str):
        pass

    lines = result.group(0).splitlines()[2:]
    forces = []
    ia = 1
    for line in lines:
        result = re.search(r'\s+atom\s+(?P<atom>\S+)\s+type\s+\S+\s+force\s+\=\s+(?P<force>.+)', line)
        if result:
            if ia != int(result.group('atom')):
                print("Wrong index of atom")
                sys.exit(-1)
            vals = result.group('force').split()
            f = []
            for i in [0, 1, 2]:
                f.append(float(vals[i]))
            forces.append(f)
            ia += 1

    return forces

def get_energy(raw_str):
    for result in re.finditer(r'!\s+total energy\s+=\s+(?P<energy>\S+) Ry', raw_str):
        pass

    return float(result['energy'])

def main():
    if len(sys.argv) != 3:
        print(f'Usage: python {sys.argv[0]} [FILE1] [FILE2]')
        sys.exit(0)

    with open(sys.argv[1], "r") as fin:
        raw_str1 = fin.read()

    with open(sys.argv[2], "r") as fin:
        raw_str2 = fin.read()

    #nat = get_num_atoms(raw_str1)

    #stress = get_stress(raw_str)
    #forces = get_forces(raw_str, nat)
    #pressure = get_pressure(raw_str)
    energy1 = get_energy(raw_str1)
    energy2 = get_energy(raw_str2)
    if abs(energy1 - energy2) > 1e-6:
        print(f'Total energy is different: ref: {energy1}, actual: {energy2}')
        sys.exit(1)

    sys.exit(0)

if __name__ == "__main__":
    main()
