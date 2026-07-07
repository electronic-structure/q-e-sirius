import yaml
import sys

def main():
    if len(sys.argv) != 3:
        print(f'Usage: python {sys.argv[0]} [FILE1] [FILE2]')
        sys.exit(0)

    with open(sys.argv[1], "r") as fin:
        d1 = yaml.safe_load(fin)

    with open(sys.argv[2], "r") as fin:
        d2 = yaml.safe_load(fin)

    #save alpha parameters in a list and sort them
    v_1 = []
    v_2 = []
    for a in d1:
      for k in d1[a]:
        # Handle if value is a str (with multiple values per line)
        if isinstance(d1[a][k], str):
          # Split string by whitespace if it contains multiple values
          vals_1 = d1[a][k].split()
          vals_2 = d2[a][k].split()
          for v1, v2 in zip(vals_1, vals_2):
            v_1.append(float(v1))
            v_2.append(float(v2))
        else:
          v_1.append(float(d1[a][k]))
          v_2.append(float(d2[a][k]))
    #sort lists (could be removed)    
    #v_1.sort()
    #v_2.sort()

    #compare lists
    for a1, a2 in zip(v_1, v_2):
      if abs(a1 - a2) > 1e-2:
        print(f'LR parameters are different: ref: {a1}, actual: {a2}')
        sys.exit(1)     
      print(f'LR parameters are ref: {a1}, actual: {a2} diff: {round(abs(a1 - a2), 5)} -> ok!')

    sys.exit(0)

if __name__ == "__main__":
    main()
