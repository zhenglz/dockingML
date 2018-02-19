from dockml import rewritePDB
import sys

if __name__ == "__main__" :

    inp = sys.argv[1]
    out = sys.argv[2]

    chain = "A"

    rew = rewritePDB(inp)

    rew.pdbRewrite(out,chain, 1, 1)

