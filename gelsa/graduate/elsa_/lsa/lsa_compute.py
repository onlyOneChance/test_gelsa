try:
    import compcore
except ImportError:
    import lsa.compcore as compcore

try:
    import lsalib
except ImportError:
    import lsa.lsalib as lsalib

import pandas as pd
import numpy as np
import scipy as sp
import sys, csv, re, os, time, argparse, string, tempfile

rpy_import = False

def main():
    # __script__ = "lsa_compute"
    # version_desc = os.popen("lsa_version").read().rstrip()
    # version_print = "%s (rev: %s) - copyright Li Charlie Xia, lixia@stanford.edu" % (__script__, version_desc) 
    # print(version_print, file=sys.stderr)

    parser = argparse.ArgumentParser()
    arg_precision_default=1000
    arg_delayLimit_default=0
    
    parser.add_argument("dataFile", metavar="dataFile", type=argparse.FileType('r'), \
        help="the input data file,\n \
        m by (r * s)tab delimited text; top left cell start with \
        '#' to mark this is the header line; \n \
        m is number of variables, r is number of replicates, \
        s it number of time spots; \n \
        first row: #header  s1r1 s1r2 s2r1 s2r2; \
        second row: x  ?.?? ?.?? ?.?? ?.??; for a 1 by (2*2) data")
    parser.add_argument("resultFile", metavar="resultFile", type=argparse.FileType('w'), \
        help="the output result file")
    parser.add_argument("-e", "--extraFile", dest="extraFile", default=None, \
        type=argparse.FileType('r'),
        help="specify an extra datafile, otherwise the first datafile will be used \n \
                and only lower triangle entries of pairwise matrix will be computed")
    parser.add_argument("-d", "--delayLimit", dest="delayLimit", default=arg_delayLimit_default, type=int,\
        help="specify the maximum delay possible, default: {},\n \
                must be an integer >=0 and <spotNum".format(arg_delayLimit_default))
    parser.add_argument("-r", "--repNum", dest="repNum", default=1, type=int,
        help="specify the number of replicates each time spot, default: 1,\n \
                must be provided and valid. ")
    parser.add_argument("-s", "--spotNum", dest="spotNum", default=4, type=int, 
        help="specify the number of time spots, default: 4,\n \
                must be provided and valid. ")
   
    parser.add_argument("-f", "--fillMethod", dest="fillMethod", default='none', \
        choices=['none', 'zero', 'linear', 'quadratic', 'cubic', 'slinear', 'nearest'], \
        help="specify the method to fill missing, default: none,               \n \
                choices: none, zero, linear, quadratic, cubic, slinear, nearest  \n \
                operation AFTER normalization:  \n \
                none: fill up with zeros ;   \n \
                operation BEFORE normalization:  \n \
                zero: fill up with zero order splines;           \n \
                linear: fill up with linear splines;             \n \
                slinear: fill up with slinear;                   \n \
                quadratic: fill up with quadratic spline;             \n \
                cubic: fill up with cubic spline;                \n \
                nearest: fill up with nearest neighbor") 

    parser.add_argument("-t", "--transFunc", dest="transFunc", default='simple', \
        choices=['simple', 'SD', 'Med', 'MAD'],\
        help="specify the method to summarize replicates data, default: simple, \n \
                choices: simple, SD, Med, MAD                                     \n \
                NOTE:                                                             \n \
                simple: simple averaging                                          \n \
                SD: standard deviation weighted averaging                         \n \
                Med: simple Median                                                \n \
                MAD: median absolute deviation weighted median;" )

    parser.add_argument("-n", "--normMethod", dest="normMethod", default='robustZ', \
        choices=['percentile', 'percentileZ', 'pnz', 'robustZ', 'rnz', 'none'], \
        help="must specify the method to normalize data, default: robustZ, \n \
                choices: percentile, none, pnz, percentileZ, robustZ or a float  \n \
                NOTE:                                                   \n \
                percentile: percentile normalization, including zeros (only with perm)\n \
                pnz: percentile normalization, excluding zeros (only with perm) \n  \
                percentileZ: percentile normalization + Z-normalization \n \
                rnz: percentileZ normalization + excluding zeros + robust estimates (theo, mix, perm OK) \n \
                robustZ: percentileZ normalization + robust estimates \n \
                (with perm, mix and theo, and must use this for theo and mix, default) \n")

    parser.add_argument("-p", "--pvalueMethod", dest="pvalueMethod", default="perm", \
        choices=["perm", "theo", "ddlsa", "mix", "bblsa", "stlta"],
        help="specify the method for p-value estimation, \n \
                default: pvalueMethod=perm, i.e. use  permutation \n \
                theo: theoretical approximaton; if used also set -a value. \n \
                mix: use theoretical approximation for pre-screening \
                if promising (<0.05) then use permutation. ")

    parser.add_argument("-x", "--precision", dest="precision", default=arg_precision_default, type=int,\
        help="permutation/precision, specify the permutation \n \
                number or precision=1/permutation for p-value estimation. \n \
                default is {}, must be an integer >0 ".format(arg_precision_default) )
                
    parser.add_argument("-m", "--minOccur", dest="minOccur", default=50, type=int, 
        help="specify the minimum occurence percentile of all times, default: 50,\n")

    parser.add_argument("-b", "--bootNum", dest="bootNum", default=0, type=int, \
        choices=[0, 100, 200, 500, 1000, 2000],
        help="specify the number of bootstraps for 95% confidence \
                interval estimation, default: 100,\n \
                choices: 0, 100, 200, 500, 1000, 2000. \n \
                Setting bootNum=0 avoids bootstrap. \n \
                Bootstrap is not suitable for non-replicated data.")
                
    parser.add_argument("-q", "--qvalu=[/eMethod", dest="qvalueMethod", default='scipy', choices=['scipy'],\
        help="specify the qvalue calculation method, \n \
                scipy: use scipy and storeyQvalue function, default \n \
                ")

    parser.add_argument("-T", "--trendThresh", dest="trendThresh", default=None, \
        type=float, \
        help="if trend series based analysis is desired, use this option \n \
                NOTE: when this is used, must also supply reasonble \n \
                values for -p, -a, -n options")

    parser.add_argument("-a", "--approxVar", dest="approxVar", default=1, type=float,\
        help="if use -p theo and -T, must set this value appropriately, \n \
                precalculated -a {1.25, 0.93, 0.56,0.13 } for i.i.d. standard normal null \n \
                and -T {0, 0.5, 1, 2} respectively. For other distribution \n \
                and -T values, see FAQ and Xia et al. 2013 in reference")

    parser.add_argument("-v", "--progressive", dest="progressive", default=0, type=int, 
        help="specify the number of progressive output to save memory, default: 0,\n \
                2G memory is required for 1M pairwise comparison. ")

    arg_namespace = parser.parse_args()
    fillMethod = vars(arg_namespace)['fillMethod']
    normMethod = vars(arg_namespace)['normMethod']
    qvalueMethod = vars(arg_namespace)['qvalueMethod']
    pvalueMethod = vars(arg_namespace)['pvalueMethod']
    precision = vars(arg_namespace)['precision']
    transFunc = vars(arg_namespace)['transFunc']
    bootNum = vars(arg_namespace)['bootNum']
    approxVar = vars(arg_namespace)['approxVar']
    trendThresh = vars(arg_namespace)['trendThresh']
    progressive = vars(arg_namespace)['progressive']
    delayLimit = vars(arg_namespace)['delayLimit']
    minOccur = vars(arg_namespace)['minOccur']
    dataFile = vars(arg_namespace)['dataFile']
    extraFile = vars(arg_namespace)['extraFile']
    resultFile = vars(arg_namespace)['resultFile']
    repNum = vars(arg_namespace)['repNum']
    spotNum = vars(arg_namespace)['spotNum']

    if trendThresh is None:
        if pvalueMethod in ["stlta"]:
            print("RunIput:pvalueMethod '{}' cannot apply analysis with trendThresh set to None(i.e. trendThresh=None).".format(pvalueMethod))
            return

    elif trendThresh is not None:
        if pvalueMethod in ["ddlsa", "bblsa"]:
            print("RunIput:pvalueMethod '{}' cannot apply analysis with trendThresh set to {}.".format(pvalueMethod, trendThresh))
            return

    try:
        extraFile_name = extraFile.name 
    except AttributeError:
        extraFile_name = ''
    assert trendThresh==None or trendThresh>=0
    if transFunc == 'SD':
        fTransform = lsalib.sdAverage
    elif transFunc == 'Med':
        fTransform = lsalib.simpleMedian
    elif transFunc == 'MAD':
        fTransform = lsalib.madMedian   
    else:
        fTransform = lsalib.simpleAverage 

    if repNum < 5 and transFunc == 'SD':
        print("Not enough replicates for SD-weighted averaging, fall back to simpleAverage", file=sys.stderr)
        transFunc = 'simple'
    if repNum < 5 and transFunc == 'MAD':
        print("Not enough replicates for Median Absolute Deviation, fall back to simpleMedian", file=sys.stderr)
        transFunc = 'Med'

    if normMethod == 'none':
        zNormalize = lsalib.noneNormalize
    elif normMethod == 'percentile':
        zNormalize = lsalib.percentileNormalize
    elif normMethod == 'percentileZ':
        zNormalize = lsalib.percentileZNormalize
    elif normMethod == 'robustZ':
        zNormalize = lsalib.robustZNormalize
    elif normMethod == 'pnz':
        zNormalize = lsalib.noZeroNormalize
    elif normMethod == 'rnz':
        zNormalize = lsalib.robustNoZeroNormalize
    else:
        zNormalize = lsalib.percentileZNormalize

    start_time = time.time()

    col = spotNum
    total_row_0 = 0
    total_row_1 = 0
    block = 5000

    next(dataFile)
    for line in dataFile:
        total_row_0 += 1

    next(extraFile)
    for line in extraFile:
        total_row_1 += 1

    if qvalueMethod in ['R'] and rpy_import:
        qvalue_func = lsalib.R_Qvalue
    else:
        qvalue_func = lsalib.storeyQvalue 

    i_m = 0
    j_m = 0
    start_0 = 1
    end_0 = block
    start_1 = 1
    end_1 = block

    if end_0 >= total_row_0:
        end_0 = total_row_0
    if end_1 >= total_row_1:
        end_1 = total_row_1

    data = compcore.LSA(total_row_0, total_row_1)
    outer_total = total_row_0 // block
    outer_desc = "total_task"
    inner_total = total_row_1 // block
    inner_desc = "inner_task"

    pall_array = []
    while i_m * block < total_row_0:
            skip_header = start_0
            skip_footer = total_row_0 - end_0
            firstData = np.genfromtxt(dataFile.name, comments='#', delimiter='\t',missing_values=['na', '', 'NA'], filling_values=np.nan,usecols=range(1,spotNum*repNum+1), skip_header=skip_header, skip_footer=skip_footer)

            if len(firstData.shape) == 1:
                firstData = np.array([firstData])

            firstFactorLabels = np.genfromtxt(dataFile.name, comments='#', delimiter='\t', usecols=range(0,1), dtype='str', skip_header=skip_header, skip_footer=skip_footer).tolist()
            if type(firstFactorLabels)==str:
                firstFactorLabels=[firstFactorLabels]

            cleanData = []
            factorNum = firstData.shape[0]
            tempData=np.zeros( ( factorNum, repNum, spotNum), dtype='float' ) 
            for i in range(0, factorNum):
                for j in range(0, repNum):
                    try:
                        tempData[i,j] = firstData[i][np.arange(j,spotNum*repNum,repNum)]
                    except IndexError:
                        print("Error: one input file need more than two data row or use -e to specify another input file", file=sys.stderr)
                        quit()
            for i in range(0, factorNum):
                for j in range(0, repNum):
                    tempData[i,j] = lsalib.fillMissing( tempData[i,j], fillMethod )
            cleanData.append(tempData)

            while j_m * block < total_row_1:
                    
                    skip_header = start_1
                    skip_footer = total_row_1 - end_1
                    secondData = np.genfromtxt(extraFile.name, comments='#', delimiter='\t',missing_values=['na', '', 'NA'], filling_values=np.nan,usecols=range(1,spotNum*repNum+1), skip_header=skip_header, skip_footer=skip_footer)

                    if len(secondData.shape) == 1:
                        secondData = np.array([secondData])

                    secondFactorLabels=np.genfromtxt(extraFile.name, comments='#', delimiter='\t', usecols=range(0,1), dtype='str', skip_header=skip_header, skip_footer=skip_footer).tolist()
                    if type(secondFactorLabels)==str:
                        secondFactorLabels=[secondFactorLabels]

                    factorNum = secondData.shape[0]
                    tempData=np.zeros((factorNum,repNum,spotNum),dtype='float')
                    for i in range(0, factorNum):
                        for j in range(0, repNum):
                            try:
                                tempData[i,j] = secondData[i][np.arange(j,spotNum*repNum,repNum)]
                            except IndexError:
                                print("Error: one input file need more than two data row or use -e to specify another input file", file=sys.stderr)
                                quit()
                    for i in range(0, factorNum):
                        for j in range(0, repNum):
                            tempData[i,j] = lsalib.fillMissing( tempData[i,j], fillMethod )
                    cleanData.append(tempData)

                    array = lsalib.palla_applyAnalysis( cleanData[0], cleanData[1], data, col, onDiag=True, delayLimit=delayLimit,bootNum=bootNum, pvalueMethod=pvalueMethod, 
                                                        precisionP=precision, fTransform=fTransform, zNormalize=zNormalize, approxVar=approxVar, resultFile=resultFile, trendThresh=trendThresh, 
                                                        firstFactorLabels=firstFactorLabels, secondFactorLabels=secondFactorLabels, qvalueMethod=qvalueMethod)
                    
                    pall_array.append(array)
                    cleanData.pop()

                    j_m += 1
                    start_1 = start_1 + block
                    end_1 = end_1 + block
                    if end_1 >= total_row_1:
                        end_1 = total_row_1

            i_m += 1
            j_m = 0
            start_1 = 1
            end_1 = block
            if end_1 >= total_row_1:
                end_1 = total_row_1

            start_0 = start_0 + block
            end_0 = end_0 + block
            if end_0 >= total_row_0:
                end_0 = total_row_0

    data_set = np.vstack(pall_array)

    lsaP = data_set[:,7]
    PCC = data_set[:,9]
    SCC = data_set[:,14]
    SPCC = data_set[:,11]
    SSCC = data_set[:,16]

    qvalues = qvalue_func(lsaP).tolist()
    pccqvalues = qvalue_func(PCC).tolist()
    sccqvalues = qvalue_func(SCC).tolist()
    spccqvalues = qvalue_func(SPCC).tolist()
    ssccqvalues = qvalue_func(SSCC).tolist()

    data_0 = np.column_stack((qvalues, pccqvalues, spccqvalues, sccqvalues, ssccqvalues))
    data_set = np.hstack((data_set, data_0))

    df1 = pd.read_csv(dataFile.name, sep='\t', index_col=0)
    rows1 = list(df1.index)
    df2 = pd.read_csv(extraFile.name, sep='\t', index_col=0)
    rows2 = list(df2.index)
    combination_list = [f"{row1},{row2}" for row1 in rows1 for row2 in rows2]
    csv_filename = arg_namespace.resultFile.name
    column_names = ['X', 'Y', 'LS', 'lowCI', 'upCI', 'Xs', 'Ys', 'Len', 'Delay', 'P', 'PCC', 'Ppcc', 'SPCC', 
    'Pspcc', 'Dspcc', 'SCC', 'Pscc', 'SSCC', 'Psscc', 'Dsscc',  'Xi', 'Yi', 'Q', 'Qpcc', 'Qspcc', 'Qscc', 'Qsscc']
    data_with_row_names = np.column_stack((combination_list, data_set))
    np.savetxt(csv_filename, data_with_row_names, delimiter=",", fmt='%s', header=",".join(column_names), comments='')    

    dataFile.close()
    extraFile.close()
    
    for i in range(0,1):
        print(lsaP[i])
    
    # print("finishing up...", file=sys.stderr)
    end_time=time.time()
    print("time elapsed %f seconds" % (end_time - start_time), file=sys.stderr)

if __name__=="__main__":
    main()
