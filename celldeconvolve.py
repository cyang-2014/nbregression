#!/usr/bin/python
import sys, os.path, pickle
import numpy
import matplotlib
matplotlib.use('Agg') # This is to avoid the need for a $DISPLAY port during plotting
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as bp
argidx = 1
njobs = -1
nsamples = 2000
obssel = None
sampledir = modelfile = datafile = samplefile = ""
while argidx < len(sys.argv) and sys.argv[argidx].startswith("-"):
    if sys.argv[argidx] in ("--help", "-h"):
        break
    elif sys.argv[argidx] == "-s":
        nsamples = int(sys.argv[argidx + 1])
        argidx += 1
    elif sys.argv[argidx] == "-m":
        modelfile = sys.argv[argidx + 1]
        argidx += 1
    elif sys.argv[argidx] == "-i":
        datafile = sys.argv[argidx + 1]
        argidx += 1
    elif sys.argv[argidx] == "-d":
        samplefile = sys.argv[argidx + 1]
        argidx += 1
    elif sys.argv[argidx] == "-o":
        sampledir = sys.argv[argidx + 1]
        argidx += 1
    elif sys.argv[argidx] == "-j":
        njobs = int(sys.argv[argidx + 1])
        argidx += 1
    elif sys.argv[argidx] == "-x":
        xs = sys.argv[argidx + 1].split(",")
        obssel= []
        for x in xs:
            x = x.split("-")
            if len(x) == 1:
                obssel.append(int(x))
            else:
                obssel.extend( range(int(x[0]), int(x[1])+1) )
        argidx += 1
    argidx += 1

if modelfile == "" or datafile == "":
    print "Usage:"
    print "stan_celldeconv.py [OPTIONS] -m MODELFILE -i INDATAFILE"
    print "-m MODELFILE      stan model definition"
    print "-i INDATAFILE     contains gene header line, predictor matrix, followed by an empty line and observed data lines"
    print "Options:"
    print "-d OBSFILE        lines of observations are in a separate file (no header line). CEF file expected if extension is '.cef'."
    print "-x L1,L2-L3,L4    process only specified (ranges of) lines of observations"
    print "-s N              number of iterations including 50% warmup (default=" + str(nsamples) + ")"
    print "-o DIR            write sampling data in this directory"
    print "-j N              use N processors for calculations (default=all)"
    sys.exit(0)

betapdffile = os.path.splitext(datafile)[0] + '_beta_plots.pdf'
betavaluefile = os.path.splitext(datafile)[0] + '_beta_medians.tab'

if sampledir != "" and not os.path.exists(sampledir):
    os.makedirs(sampledir)


class StanData:
    def __init__(self, datafile, samplefile = None, obssel = None):
        self.cellids = []
        self.celldata = []
        self.obsids = []
        self.obsdata = []
        fd = open(datafile)
        line = fd.readline().strip()
        while line.startswith("#"):
            line = fd.readline().strip()
        self.genenames = line.split("\t")
        cline = fd.readline().strip()
        while cline and cline != "" and not cline.startswith("#"):
            cdata = cline.split("\t")
            cellid = cdata.pop(0)
            cdata = [float(v) for v in cdata]
            self.cellids.append(cellid)
            self.celldata.append(cdata)
            cline = fd.readline().strip()
        if samplefile:
            fd.close()
            fd = open(samplefile)
        oline = fd.readline()
        if samplefile.lower().endswith(".cef"):
            fields = oline.split('\t')
            for n in range(int(fields[1]) + int(fields[3]) + 2):
                oline = fd.readline()
        while oline.strip().startswith("#"):
            oline = fd.readline()
        if oline.startswith("\t"): # If there is a header in observation section/file, remove it
            oline = fd.readline()
        olineidx = 1
        oline = oline.strip()
        while oline:
            if obssel == None or olineidx in obssel:
                odata = oline.split("\t")
                obsid = odata.pop(0)
                odata = [int(v) for v in odata]
                self.obsids.append(obsid)
                self.obsdata.append(odata)
            olineidx += 1
            oline = fd.readline().strip()
        fd.close()
        print datafile, ":", self.Npredictors(), "predictors X", self.Noutcomes(), "outcomes."
        sf = samplefile if samplefile else datafile
        print sf, ":", self.Nobservations(), "observations to process."

    def Noutcomes(self):
        return len(self.genenames)
    def Npredictors(self):
        return len(self.cellids)
    def predictors(self):
        return numpy.transpose(self.celldata)
    def Nobservations(self):
        return len(self.obsids)
    def iter_observations(self):
        for obsid, odata in zip(self.obsids, self.obsdata):
            yield obsid, odata


standata = StanData(datafile, samplefile, obssel)
N = standata.Noutcomes()
K = standata.Npredictors()
matrix = standata.predictors()
print "Predictor matix (", K, "X", N, "):"
for i in range(min(K, 5)):
    print "\t".join( [str(matrix[j,i]) for j in range(min(N, 10))] )

from pystan import StanModel
smfile = modelfile + ".pkl"
if not os.path.exists(smfile):
    modelname = os.path.splitext(os.path.basename(modelfile))[0] + "_" + os.path.splitext(os.path.basename(datafile))[0]
    stanmodel = StanModel(file=modelfile, model_name=modelname)
    with open(smfile, 'wb') as f:
        pickle.dump(stanmodel, f)
else:
    stanmodel = pickle.load(open(smfile, 'rb'))

def get_median(sample_array):
        sample_array.sort()
        n = len(sample_array)
        medianvalue = sample_array[n/2] if (n % 2) == 0 else (sample_array[n/2] + sample_array[n/2 + 1]) / 2.0
        return medianvalue

pdf = bp.PdfPages(betapdffile)
results = {}
for obsid, odata in standata.iter_observations():
    sample_outfile = os.path.join(sampledir, obsid + "_samples.txt") if sampledir != "" else None
    sdata = { "N":N, "K":K, "x":matrix, "y":odata }
    fit = stanmodel.sampling(data=sdata, iter=nsamples, n_jobs=njobs, sample_file=sample_outfile)
    pars = fit.extract(["beta"])
    betasamples = pars["beta"]
    results[obsid] = []
    samples = []
    for i in range(K):
        ith_samples = [ bs[i] for bs in betasamples ]
        medianvalue = get_median(ith_samples)
        results[obsid].append(medianvalue)
        samples.append(ith_samples)
    f = plt.figure()
    plt.boxplot(samples) #, labels=standata.cellids)
    plt.xticks(numpy.arange(1 + len(standata.cellids)), ["0"] + standata.cellids, rotation=90)
    f.suptitle(obsid)
    pdf.savefig(f)
    #fit.plot()
pdf.close()
f.clf()

fdb = open(betavaluefile, 'w')
fdb.write("Median of beta values:\n")
fdb.write("Sample\t" + "\t".join(standata.cellids) + "\n")
for obsid in standata.obsids:
    result = results[obsid]
    fdb.write(obsid + "\t" + "\t".join( [str(v) for v in result] ) + "\n")
fdb.close()

print "Wrote beta values to", betavaluefile
print "Wrote beta value boxplots to", betapdffile
