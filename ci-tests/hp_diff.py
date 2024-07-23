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

    for a in d1:
      for k in d1[a]:
        if abs(d1[a][k] - d2[a][k]) > 1e-2:
          print(f'LR parameters are different: ref: {d1[a][k]}, actual: {d2[a][k]}')
          sys.exit(1)     
        print(f'LR parameters are ref: {d1[a][k]}, actual: {d2[a][k]} diff: {abs(d1[a][k] - d2[a][k])} -> ok!')

    sys.exit(0)

if __name__ == "__main__":
    main()
