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
    alpha_1 = []
    alpha_2 = []
    for a in d1:
      for k in d1[a]:
        alpha_1.append(d1[a][k])
        alpha_2.append(d2[a][k])
    alpha_1.sort()
    alpha_2.sort()

    #compare lists
    for a1, a2 in zip(alpha_1, alpha_2):
      if abs(a1 - a2) > 1e-2:
        print(f'LR parameters are different: ref: {a1}, actual: {a2}')
        sys.exit(1)     
      print(f'LR parameters are ref: {a1}, actual: {a2} diff: {abs(a1 - a2)} -> ok!')

    sys.exit(0)

if __name__ == "__main__":
    main()
