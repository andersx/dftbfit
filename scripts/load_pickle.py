import cPickle
import numpy as np

def load_pickle(filename):
    f = open(filename,"rb")
    p = cPickle.load(f)
    f.close()
    return(p)

co = load_pickle("charges_orca.pickle")
cg = load_pickle("charges_gaussian.pickle")



np.set_printoptions(formatter={'float': '{: 0.3f}'.format}, linewidth=1000000)

for i in sorted(co):

    try:
        qdiff = np.array(co[i]) - np.array(cg[i])
        max_qdiff = max(qdiff.min(), qdiff.max(), key=abs)
        print "%-30s  %7.4f" % (i, max_qdiff), qdiff
    except:
        print "%-30s  %7s" % (i, "FAIL!")

