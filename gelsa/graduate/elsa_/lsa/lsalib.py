try:
    import compcore
except ImportError:
    import lsa.compcore as compcore

import warnings
from itertools import chain
import csv, sys, os, random, tempfile, time
import multiprocessing
import numpy as np
import numpy.testing
import scipy as sp
import scipy.interpolate
import scipy.stats
import scipy.linalg
import statsmodels.api as sm

rpy_import = False

disp_decimal=8
kcut_min=100
Rmax_min=10
Rmax_max=50
my_decimal = 2
pipi = np.pi**2
pipi_inv = 1/pipi
Q_lam_step = 0.05
Q_lam_max = 0.95
                   

def fillMissing(tseries, method):
  if method == 'none':
    return tseries
  else:
    y = tseries[np.logical_not(np.isnan(tseries))]
    x = np.array(range(0, len(tseries)), dtype='float')[np.logical_not(np.isnan(tseries))]
    try:
      spline_fit = sp.interpolate.interp1d( x, y, method )
    except:
      return tseries
    yy = np.zeros( len(tseries), dtype='float' )
    for i in range(0, len(tseries)):
      if not np.isnan(tseries[i]):
        yy[i] = tseries[i]
      else:
        try:
          yy[i] = spline_fit(i)
        except ValueError:
          yy[i] = tseries[i]
    return yy

def simpleAverage(tseries):
  Xf = ma_average(tseries, axis=0)
  return Xf

def sdAverage(tseries):
  try:
    sd = np.ma.std(tseries,axis=0,ddof=1)
  except FloatingPointError:
    return simpleAverage(tseries)
  if np.any(sd.mask) or (np.ma.sum(sd==0))>0:
    return simpleAverage(tseries)
  Xf = ma_average(tseries, axis=0)*(1/sd)*(1/np.ma.sum(1/sd))*(1/sd)
  return Xf

def simpleMedian(tseries):
  Xf = ma_median(tseries, axis=0)
  return Xf

def madMedian(tseries):
  Xf = tseries
  mad = ma_median( np.ma.abs(Xf - ma_median(Xf, axis=0)), axis=0 )
  if np.any(mad.mask) or (np.ma.sum(mad==0))>0:
    return simpleMedian(tseries)
  Xf = ma_median(Xf, axis=0)*(1/mad)*(1/np.ma.sum(1/mad))*(1/mad)
  return Xf

def percentileNormalize(tseries):
  ranks = tied_rank(tseries)
  nt = np.ma.masked_invalid(sp.stats.distributions.norm.ppf( ranks/(len(ranks)-np.sum(ranks.mask)+1) ),copy=True)
  nt = nt.filled(fill_value=0)
  return nt

def percentileZNormalize(tseries):
  ranks = tied_rank(tseries)
  nt = np.ma.masked_invalid(sp.stats.distributions.norm.ppf( ranks/(len(ranks)-np.sum(ranks.mask)+1) ),copy=True)
  try:
    zt = (nt - np.ma.mean(nt, axis=0))/np.ma.std(nt)
  except FloatingPointError:
    zt = nt - np.ma.mean(nt, axis=0)
  zt = zt.filled(fill_value=0)
  return zt

def robustZNormalize(tseries):
  ranks = tied_rank(tseries)
  nt = np.ma.masked_invalid(sp.stats.distributions.norm.ppf( ranks/(len(ranks)-np.sum(ranks.mask)+1) ),copy=True)
  mad_sd = 1.4826 * np.ma.median( np.ma.abs(nt - np.ma.median(nt)))
  range_sd = (np.ma.max(nt) - np.ma.min(nt))/4
  sd_est = range_sd if mad_sd == 0 else mad_sd 
  try:
    zt = (nt - np.ma.median(nt))/sd_est
  except FloatingPointError:
    zt = nt - np.ma.median(nt)
  zt = zt.filled(fill_value=0)
  return zt

def noZeroNormalize(tseries):
  nt = np.ma.masked_equal(tseries, 0)
  if type(nt.mask) == np.bool_:
    nt.mask = [nt.mask] * nt.shape[0]
  ranks = tied_rank(nt)
  nt = np.ma.masked_invalid(sp.stats.distributions.norm.ppf( ranks/(len(ranks)-np.sum(ranks.mask)+1) ),copy=True)
  try:
    zt = (nt - np.ma.mean(nt, axis=0))*(1/np.ma.std(nt, axis=0))
  except FloatingPointError:
    zt = nt - np.ma.mean(nt, axis=0)
  zt = np.ma.masked_invalid(zt)
  zt = zt.filled(fill_value=0)
  return zt

def noneNormalize(tseries):
  nt = tseries.filled(fill_value=0)
  return nt

def ji_calc_trend(oSeries, lengthSeries, thresh):
  tSeries = np.zeros(lengthSeries-1, dtype='float')
  for i in range(0, lengthSeries-1):
    if oSeries[i] == 0 and oSeries[i+1] > 0:
      trend = 1
    elif oSeries[i] == 0 and oSeries[i+1] < 0:
      trend = -1
    elif oSeries[i] == 0 and oSeries[i+1] == 0:
      trend = 0
    else:
      trend = (oSeries[i+1]-oSeries[i])/np.abs(oSeries[i])

    if np.isnan(trend):
      tSeries[i] = np.nan
    elif trend >= thresh:
      tSeries[i] = 1
    elif trend <= -thresh:
      tSeries[i] = -1
    else:
      tSeries[i] = 0
  return tSeries

def rpy_spearmanr(Xz, Yz):
  try:
    sr=r('''cor.test''')(Xz,Yz,method='spearman')
    return (sr[3][0],sr[2][0])
  except rpy2.rinterface.RRuntimeError:
    return (np.nan,np.nan)

def rpy_pearsonr(Xz, Yz):
  try:
    sr=r('''cor.test''')(Xz,Yz,method='pearson')
    return (sr[3][0],sr[2][0])
  except rpy2.rinterface.RRuntimeError:
    return (np.nan,np.nan)

def scipy_spearmanr(Xz, Yz):
  try:
    return scipy.stats.spearmanr(Xz, Yz)
  except:
    return (np.nan,np.nan)

def scipy_pearsonr(Xz, Yz):
  try:
    return scipy.stats.pearsonr(Xz, Yz)
  except:
    return (np.nan,np.nan)

def calc_spearmanr(Xz, Yz, sfunc=rpy_spearmanr):
  if not rpy_import:
    sfunc = scipy_spearmanr
  mask = np.logical_or(Xz.mask, Yz.mask)
  Xz.mask = mask
  Yz.mask = mask
  (SCC, P_SCC) = sfunc(Xz.compressed(), Yz.compressed())
  return (SCC, P_SCC)

def calc_pearsonr(Xz, Yz, pfunc=rpy_pearsonr):
  if not rpy_import:
    pfunc = scipy_pearsonr
  mask = np.logical_or(Xz.mask, Yz.mask)
  Xz.mask = mask
  Yz.mask = mask
  (PCC, P_PCC) = pfunc(Xz.compressed(), Yz.compressed())
  return (PCC, P_PCC)

def ma_average(ts, axis=0):
  ts.mask
  ns = np.ma.mean(ts, axis=0)
  if type(ns.mask) == np.bool_: 
    ns.mask = [ns.mask] * ns.shape[axis]
  return ns

def calc_shift_corr(Xz, Yz, D, corfunc=calc_pearsonr): 
  d_max=0
  r_max=0
  p_max=1
  for d in range(-D, D+1):
    if d < 0:
      X_seg = Xz[:(len(Xz)+d)]
      Y_seg = Yz[-d:len(Yz)]
    elif d == 0:
      X_seg = Xz
      Y_seg = Yz
    else:
      X_seg = Xz[d:len(Xz)]
      Y_seg = Yz[:len(Yz)-d]
    assert len(X_seg) == len(Y_seg)
    mask = np.logical_or(X_seg.mask, Y_seg.mask)
    X_seg.mask = mask
    Y_seg.mask = mask
    cor = corfunc(X_seg, Y_seg)
    if np.abs(cor[0]) >= np.abs(r_max):
      r_max = cor[0] 
      d_max = d
      p_max = cor[1]
  return (r_max, p_max, d_max)

def acf(x):
    n = len(x)
    acf = np.zeros(n)
    for lag in range(n):
        data1 = x[0:n-lag]
        data2 = x[lag:n]
        acf[lag] = np.dot(data1 - np.mean(x), data2 - np.mean(x))/ (n)
    return acf

def scale(data):
    data_array = np.array(data)
    mean_val = np.mean(data_array)
    std_dev = np.std(data_array)
    scaled_data = []

    for item in data_array:
        scaled_item = (item - mean_val) / std_dev
        scaled_data.append(scaled_item)

    return np.array(scaled_data)


#===========================================================================
#1.在这里更改stlsa，将相应的方差进行计算
num_states = 3

def discretization(x, t):
    n = len(x)
    b = np.diff(x) / np.abs(x[:-1])
    b[b >= t] = 1
    b[b <= (-1) * t] = -1
    b[(b > -t) & (b < t)] = 0
    for i in range(len(b)):
        if np.isnan(b[i]):
            b[i] = np.sign(x[i + 1])
    return b

def markov_matrix(num_states,chain_x):
    transition_counts = np.zeros((num_states, num_states))
    for i in range(1, len(chain_x)):
        current = chain_x[i - 1] + 1
        back = chain_x[i] + 1
        transition_counts[current, back] += 1
    return transition_counts / transition_counts.sum(axis=1, keepdims=True)

def markovchain_sigma_2d(ax, ay):
    return 1 + 2 * (2 * ax - 1) * (2 * ay - 1) / (1 - (2 * ax - 1) * (2 * ay - 1))

def markovchain_sigma_3d(bx, dx, cx, by, dy, cy):
    return 4 * (dx / (1 - bx - cx + 2 * dx)) * (dy / (1 - by - cy + 2 * dy)) * (1 + 2 * (bx - cx) * (by - cy) / (1 - (bx - cx) * (by - cy)))

def markovchain_sigma_2d_3d(ax, by, dy, cy):
    return 2 * (dy / (1 - by - cy + 2 * dy)) * (1 + 2 * (2 * ax - 1) * (by - cy) / (1 - (2 * ax - 1) * (by - cy)))

def Sigma(series1, series2, fTransform, zNormalize, trendThresh):
    timespots = series1.shape[1]
    x = []
    y = []
    x = ji_calc_trend(zNormalize(fTransform(series1)), timespots, trendThresh)
    y = ji_calc_trend(zNormalize(fTransform(series2)), timespots, trendThresh)
    if trendThresh > 0:
      num_states = 3
    else:
      num_states = 2
    x = x.astype(int)
    y = y.astype(int)
    T_X = markov_matrix(num_states,x)
    T_Y = markov_matrix(num_states,y)
    
    if T_X.shape[0] == 2 and T_Y.shape[0] == 2:
        ax = (T_X[0, 0] + T_X[1, 1]) / 2
        ay = (T_Y[0, 0] + T_Y[1, 1]) / 2
        return markovchain_sigma_2d(ax, ay)
    elif T_X.shape[0] == 3 and T_Y.shape[0] == 3:
        bx = (T_X[0, 0] + T_X[2, 2]) / 2
        dx = (T_X[1, 0] + T_X[1, 2]) / 2
        cx = (T_X[2, 0] + T_X[0, 2]) / 2
        by = (T_Y[0, 0] + T_Y[2, 2]) / 2
        dy = (T_Y[1, 0] + T_Y[1, 2]) / 2
        cy = (T_Y[2, 0] + T_Y[0, 2]) / 2
        return markovchain_sigma_3d(bx, dx, cx, by, dy, cy)
    elif T_X.shape[0] == 3 and T_Y.shape[0] == 2:
        bx = (T_X[0, 0] + T_X[2, 2]) / 2
        dx = (T_X[1, 0] + T_X[1, 2]) / 2
        cx = (T_X[2, 0] + T_X[0, 2]) / 2
        ay = (T_Y[0, 0] + T_Y[1, 1]) / 2
        return markovchain_sigma_2d_3d(ay, bx, dx, cx)
    elif T_X.shape[0] == 2 and T_Y.shape[0] == 3:
        ax = (T_X[0, 0] + T_X[1, 1]) / 2
        by = (T_Y[0, 0] + T_Y[2, 2]) / 2
        dy = (T_Y[1, 0] + T_Y[1, 2]) / 2
        cy = (T_Y[2, 0] + T_Y[0, 2]) / 2
        return markovchain_sigma_2d_3d(ax, by, dy, cy)

#===========================================================================
def Omega(series1, series2, fTransform, zNormalize):
    Xz = zNormalize(fTransform(series1))
    Yz = zNormalize(fTransform(series2))
    y = Yz
    x = Xz
    timepoints = len(x)
    z = x * y
    residual = z - np.mean(z)
    rhohat = np.sum(residual[1:] * residual[:-1]) / np.sum(residual[:-1] ** 2)
    alphahat = 4 * rhohat ** 2 / (1 - rhohat ** 2) ** 2
    bandwidth = int(np.ceil(1.1447 * (alphahat * timepoints) ** (1 / 3)))
    covx = acf(x)
    covy = acf(y)
    covz = covx * covy
    
    omega = 0
    if not np.isnan(bandwidth):
      if bandwidth > 1:
        for x in range(-bandwidth+1, bandwidth):
          abs_x = abs(x)
          omega += covz[abs_x] * (1 - abs(x / bandwidth))
      else:
        omega = covz[0]
      if omega < 0:
        omega = variance = np.var(x) * np.var(y)
    return omega

def readPvalue(x, y, P_table, R, N, x_sd=1., M=1., alpha=1., beta=1., x_decimal=my_decimal, trendThresh=None):
  try:
    if manager.pvalueMethod in ['ddlsa']:
      xOmegay = Omega(x, y, manager.fTransform, manager.zNormalize)
      xi = int(np.around(R*M/(xOmegay*np.sqrt(alpha*beta*N))*(10**x_decimal)))
    elif manager.pvalueMethod in ['stlta']:
      sigma = Sigma(x, y, manager.fTransform, manager.zNormalize, trendThresh)
      xi = int(np.around(R*M/(sigma*np.sqrt(alpha*beta*N))*(10**x_decimal)))
    else:
      xi = int(np.around(R*M/(x_sd*np.sqrt(alpha*beta*N))*(10**x_decimal)))
  except (OverflowError, ValueError) as e:
    return np.nan
    
  if xi in P_table:
    return P_table[xi]
  elif xi>max(P_table.keys()):
    return 0.
  else:
    return np.nan

def theoPvalue(Rmax, Dmax=0, precision=.001, x_decimal=my_decimal):
  Rmax = np.max((Rmax, Rmax_min))
  Rmax = np.min((Rmax, Rmax_max))
  P_table = dict()
  for xi in range(0,Rmax*10**(x_decimal)+1): 
    if xi == 0:
      P_table[xi] = 1
      continue
    x = xi/float(10**(x_decimal))
    xx = x**2
    pipi_over_xx = pipi/xx
    alpha = precision
    B = 2*Dmax+1

    Kcut = np.max((kcut_min, int(np.ceil(.5 - np.log((alpha/(2**B-1))**(1/B)*xx*(1-np.exp(-pipi_over_xx))/8/2)/pipi_over_xx))))
    A = 1/xx
    Rcdf = 0
    P_two_tail = 1.
    for k in range(1,Kcut+1):
      C = (2*k-1)**2
      Rcdf = Rcdf + (A+pipi_inv/C)*np.exp(-C*pipi_over_xx/2)
      P_current = 1 - (8**B)*(Rcdf**B)
      if P_current<0:
        break
      else:
        P_two_tail = P_current
    P_table[xi] = P_two_tail
  return P_table

def permuPvalue(series1, series2, delayLimit, col, precisionP, \
    Smax, fTransform, zNormalize, trendThresh=None):
  # warnings.filterwarnings("ignore", category=RuntimeWarning)
  lengthSeries = series1.shape[1]
  timespots = series1.shape[1]
  if trendThresh != None:
    lengthSeries = timespots - 1
  PP_set = np.zeros(precisionP, dtype='float')
  x = []
  y = []

  if trendThresh == None:
    Xz = zNormalize(fTransform(series1))
  else:
    Xz = ji_calc_trend(zNormalize(fTransform(series1)), timespots, trendThresh)
  Y = np.ma.array(series2)

  for i in range(0, precisionP):
    np.random.shuffle(Y.T)
    if trendThresh == None:
      Yz = zNormalize(fTransform(Y))
    else:
      Yz = ji_calc_trend(zNormalize(fTransform(Y)), timespots, trendThresh)
    y.extend(Yz)
  x = Xz

  data = compcore.LSA(1, precisionP)
  data.assign(1, precisionP, delayLimit, col, x, y)
  data.dp_lsa()
  PP_set = np.array(data.score)
  if Smax >= 0:
    P_two_tail = np.sum(np.abs(PP_set) >= Smax)/np.float64(precisionP)
  else:
    P_two_tail = np.sum(-np.abs(PP_set) <= Smax)/np.float64(precisionP)
  return P_two_tail

def BlockResample(x):
    model = sm.tsa.ARIMA(x, order=(1, 0, 0))
    results = model.fit()
    EstimatedCoefficient = results.params[0]
    lopt = max(int(np.ceil((6**(1/2) * abs(EstimatedCoefficient) / (1 - EstimatedCoefficient**2))**(2/3) * len(x)**(1/3) - 0.5)), 1)
    block_number = len(x) // lopt + 1
    block_index = np.random.randint(max(len(x) - lopt + 1, 1), size=block_number)
    perm_x = np.zeros(len(x))
    
    if block_number == 1:
        middle = np.random.randint(len(x))
        return np.concatenate((x[middle:], x[:middle]))
    else:
        for i in range(block_number):
            perm_x[lopt*(i):lopt*(i+1)] = x[block_index[i]:block_index[i]+lopt]
        return perm_x[:len(x)]

def MovingBlockBootstrap(series1, series2, delayLimit, col, precisionP, Smax, fTransform, zNormalize):
  x = []
  y = []

  Xz = zNormalize(fTransform(series1))
  Yz = zNormalize(fTransform(series2))
  x = Xz
  for i in range(0, precisionP):
    y.extend(BlockResample(Yz))
  data = compcore.LSA(1, precisionP)
  data.assign(1, precisionP, delayLimit, col, x, y)
  data.dp_lsa()
  PP_set = np.array(data.score)
  if Smax >= 0:
    P_two_tail = np.sum(np.abs(PP_set) >= Smax)/np.float64(precisionP)
  else:
    P_two_tail = np.sum(-np.abs(PP_set) <= Smax)/np.float64(precisionP)
  return P_two_tail

def bootstrapCI(series1, series2, Smax, delayLimit, col, bootCI, bootNum, \
    fTransform, zNormalize, trendThresh=None, debug=0):
  if series1.shape[0] == 1:
    return (Smax, Smax, Smax)
  assert series1.shape[1] == series2.shape[1]
  timespots = series1.shape[1]
  lengthSeries = series1.shape[1]
  if trendThresh != None:
    lengthSeries = timespots - 1
  x = []
  y = []
  BS_set = np.zeros(bootNum, dtype='float')
  for i in range(0, bootNum):
    if trendThresh == None:
      Xb = zNormalize(fTransform(np.ma.array([sample_wr(series1[:,j], series1.shape[0]) for j in range(0,series1.shape[1])]).T))
      Yb = zNormalize(fTransform(np.ma.array([sample_wr(series2[:,j], series2.shape[0]) for j in range(0,series2.shape[1])]).T))
    else:
      Xb = ji_calc_trend( zNormalize(fTransform(np.ma.array([ sample_wr(series1[:,j], series1.shape[0]) for j in range(0,series1.shape[1]) ]).T)), timespots, trendThresh )
      Yb = ji_calc_trend( zNormalize(fTransform(np.ma.array([ sample_wr(series2[:,j], series2.shape[0]) for j in range(0,series2.shape[1]) ]).T)), timespots, trendThresh )
    x.extend(Xb)
    y.extend(Yb)

  data = compcore.LSA(bootNum,bootNum)
  data.assign(bootNum,bootNum, delayLimit, col, x, y)
  data.dp_lsa()
  PP_set = data.score
  BS_set = [PP_set[i ** 2 - 1] for i in range(bootNum)]
  BS_set.sort()
  BS_mean = np.mean(BS_set)
  a1 = (1-bootCI)/2.0
  a2 = bootCI+(1-bootCI)/2.0
  if debug in [1, 3]:
    Smax = 2*Smax - BS_mean
  if debug in [2, 3]:
    z0 = sp.stats.distributions.norm.ppf(np.sum(BS_set <= Smax)/float(bootNum))
    a1 = sp.stats.distributions.norm.cdf(2*z0+sp.stats.distributions.norm.ppf(a1))
    a2 = sp.stats.distributions.norm.cdf(2*z0+sp.stats.distributions.norm.ppf(a2))
  return ( BS_mean, BS_set[int(np.floor(bootNum*a1))-1], BS_set[int(np.ceil(bootNum*a2))-1] )

def sample_wr(population, k):
  n = len(population)
  _random, _int = random.random, int
  result = np.array([np.nan] * k)
  for i in range(k):
    j = _int(_random() * n)
    if type(population) == np.ma.MaskedArray:
      if population.mask[j]:
        result[i] = np.nan
      else:
        result[i] = population[j]
    else:
      result[i] = population[j]
  if type(population) == np.ma.MaskedArray:
    result = np.ma.masked_invalid(result)
  return result

def ma_median(ts, axis=0):
  ns = np.ma.median(ts, axis=axis)
  if type(ns.mask) == np.bool_:       #fix broken ma.median, mask=False instead of [False, ...] for all mask
    ns.mask = [ns.mask] * ns.shape[axis]
  return ns

def R_Qvalue(pvalues, lam=np.arange(0,Q_lam_max,Q_lam_step), method='smoother', robust=False, smooth_df=3):
  try:
    pvalues_not_nan = np.logical_not(np.isnan(pvalues))
    pvalues_input = pvalues[pvalues_not_nan] 
    r('''library(qvalue)''')
    
    qvalues=r['''qvalue'''](p=pvalues_input, **{'lambda':lam, 'pi0.method':method, 'robust':robust, 'smooth.df':3})
    qvalues_return = [np.nan]*len(pvalues)
    
    j = 0
    for i in range(0, len(pvalues)):
      if not np.isnan(pvalues[i]):
        qvalues_return[i]=qvalues[2][j]
        j = j+1
  except:
    print >>sys.stderr, "from R: unusable pvalues -> ", pvalues_input
    qvalues_return=[np.nan]*len(pvalues)

  return qvalues_return

def storeyQvalue(pvalues, lam=np.arange(0,Q_lam_max,Q_lam_step), method='smoother', robust=False, smooth_df=3):
  try:
    mpvalues = np.ma.masked_invalid(pvalues,copy=True)
    rpvalues = mpvalues[~mpvalues.mask]
    p_num = len(pvalues)
    rp_num = len(rpvalues)

    if rp_num <= 1:
      return np.array( [np.nan] * p_num, dtype='float')

    rp_max = np.max(rpvalues)
    rp_lam = lam[lam<rp_max]

    if len(rp_lam) <= 1:
      return np.array( [np.nan if np.isnan(pvalues[i]) else 0 for i in range(0,p_num)],  dtype='float')

    pi_set = np.zeros(len(rp_lam), dtype='float')
    for i in range(0, len(rp_lam)): 
      pi_set[i] = np.mean(rpvalues>=rp_lam[i])/(1-rp_lam[i])
      
    if method=='smoother':
      spline_fit = sp.interpolate.interp1d(rp_lam, pi_set, kind=smooth_df)
      pi_0 = spline_fit(np.max(rp_lam))
      pi_0 = np.max( [np.min( [np.min(pi_0), 1]), 0] ) 
      if pi_0 == 0:
        method='bootstrap'
    if method=='bootstrap':
      pi_min = np.min(pi_set)
      mse = np.zeros((100, len(rp_lam)), dtype='float')
      pi_set_boot = np.zeros((100, len(rp_lam)), dtype='float')
      for j in range(0, 100):
        p_boot = sample_wr(rpvalues, rp_num)
        for i in range(0, len(rp_lam)):
          pi_set_boot[j][i] = np.mean(p_boot>=rp_lam[i])/(1-rp_lam[i]) 
        mse[j] = (pi_set_boot[j]-pi_min)**2
      min_mse_j = np.argmin(mse)
      pi_0 = np.min(pi_set_boot[min_mse_j])
      pi_0 = np.max([np.min( [np.min(pi_0), 1]), ])
      if pi_0 == 0:
        pi_0 = Q_lam_step 

    rp_argsort = np.argsort(rpvalues)
    rp_ranks = tied_rank(rpvalues)
    if robust:
      rqvalues = pi_0*rp_num*rpvalues*(1/(rp_ranks*(1-np.power((1-rpvalues),rp_num))))
    else:
      rqvalues = pi_0*rp_num*rpvalues*(1/rp_ranks) 
    rqvalues[rp_argsort[rp_num-1]] = np.min( [rqvalues[rp_argsort[rp_num-1]], 1] ) 
    for i in reversed(range(0,rp_num-1)): 
      rqvalues[rp_argsort[i]] = np.min( [rqvalues[rp_argsort[i]], rqvalues[rp_argsort[i+1]], 1] )

    qvalues=np.array([np.nan]*p_num)
    j=0
    for i in range(0, p_num):
      if not mpvalues.mask[i]:
        qvalues[i]=rqvalues[j]
        j += 1
  except:
    qvalues=np.array( [np.nan] * p_num, dtype='float')
  return qvalues

def tied_rank(values):
  assert type(values) == np.ma.MaskedArray
  V = np.ma.asarray(values)
  nans = (np.nonzero(V.mask)[0]).tolist()
  V = V[~V.mask]
  v_num = {}
  v_cum = {}
  for v in V:
    v_num[v] = v_num.get(v,0) + 1
  suvs = v_num.keys()
  suvs = list(suvs)
  suvs.sort()

  c = 0
  for v in suvs:
    c += v_num[v]
    v_cum[v] = c
  sV = np.array( [ (2*v_cum[V[i]]-v_num[V[i]]+1)/2 for i in range(0, len(V)) ], dtype='float' )
  for idx in nans:
    sV = np.insert(sV, idx, np.nan)
  sV = np.ma.masked_invalid(sV,copy=True)
  return sV 

def myfunc(i):
  num_1 = manager.num_1
  start = i * num_1
  end = (i+1) * num_1
  X = [0] * num_1
  Y = [0] * num_1
  Smax = manager.Smax[start:end]
  Sl = [0.0] * num_1
  Su = [0.0] * num_1
  Xs = manager.Xs[start:end]
  Ys = manager.Ys[start:end]
  Al = manager.Al[start:end]
  delay = manager.delay[start:end]

  lsaP = [0.0] * num_1
  PCC = [0.0] * num_1
  P_PCC = [0.0] * num_1
  SPCC = [0.0] * num_1
  P_SPCC = [0.0] * num_1
  D_SPCC = [0.0] * num_1
  SCC = [0.0] * num_1
  P_SCC = [0.0] * num_1
  SSCC = [0.0] * num_1
  P_SSCC = [0.0] * num_1
  D_SSCC = [0.0] * num_1

  for j in range(0, num_1):
    (X[j], Y[j]) = (i, j)
    if np.all(manager.y_series[j].mask) or np.all(manager.x_series[i].mask) or manager.x_occur[i] or manager.y_occur[j]:
      (Smax[j], Sl[j], Su[j], Xs[j], Ys[j], Al[j], delay[j], lsaP[j], PCC[j], 
      P_PCC[j], SPCC[j], P_SPCC[j], D_SPCC[j], SCC[j], P_SCC[j], SSCC[j], P_SSCC[j], D_SSCC[j]) = \
      (0, 0, 0, -2, -2, 0, 0, 1, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan, np.nan)
    else:
      (PCC[j], P_PCC[j]) = calc_pearsonr(ma_average(manager.x_series[i], axis=0), ma_average(manager.y_series[j], axis=0))
      (SCC[j], P_SCC[j]) = calc_spearmanr(ma_average(manager.x_series[i], axis=0), ma_average(manager.y_series[j], axis=0))
    try:
      (SPCC[j], P_SPCC[j], D_SPCC[j]) = calc_shift_corr(ma_average(manager.x_series[i], axis=0), ma_average(manager.y_series[j], axis=0), manager.delayLimit, calc_pearsonr) 
      (SSCC[j], P_SSCC[j], D_SSCC[j]) = calc_shift_corr(ma_average(manager.x_series[i], axis=0), ma_average(manager.y_series[j], axis=0), manager.delayLimit, calc_spearmanr) 
    except FloatingPointError:
      (SPCC[j], P_SPCC[j], D_SPCC[j]) = (np.nan, np.nan, np.nan)
      (SSCC[j], P_SSCC[j], D_SSCC[j]) = (np.nan, np.nan, np.nan)

    if Al[j] == 0:
      (Smax[j], Sl[j], Su[j], Xs[j], Ys[j], Al[j], delay[j], lsaP[j]) = (0, 0, 0, -1, -1, 0, 0, 1)
    promisingP = 0.05
    if manager.pvalueMethod in ['theo', 'ddlsa', 'mix', 'stlta']:
      Xp = np.ma.array(manager.x_series[i],copy=True)
      Yp = np.ma.array(manager.y_series[j],copy=True)
      lsaP[j] = readPvalue(Xp, Yp,manager.P_table, R=np.abs(Smax[j])*manager.lengthSeries, N=manager.lengthSeries, 
      x_sd=manager.stdX, M=manager.replicates, alpha=1.,beta=1., x_decimal=my_decimal, trendThresh=manager.trendThresh)

    if (manager.pvalueMethod in ['mix'] and manager.lsaP[j]<=promisingP) or (manager.pvalueMethod in ['perm']):
      Xp = np.ma.array(manager.x_series[i],copy=True)
      Yp = np.ma.array(manager.y_series[j],copy=True)
      lsaP[j] = permuPvalue(Xp, Yp, manager.delayLimit, manager.col, manager.precisionP, 
      np.abs(Smax[j]),manager.fTransform, manager.zNormalize, manager.trendThresh)

    if manager.pvalueMethod in ['bblsa']:
      Xp = np.ma.array(manager.x_series[i],copy=True)
      Yp = np.ma.array(manager.y_series[j],copy=True)
      lsaP[j] = MovingBlockBootstrap(Xp, Yp, manager.delayLimit, manager.col, 
      manager.precisionP, np.abs(Smax[j]),manager.fTransform, manager.zNormalize)

    if manager.bootNum > 0:
      Xb = np.ma.array(manager.x_series[i],copy=True)
      Yb = np.ma.array(manager.y_series[j],copy=True)
      (Smax[j], Sl[j], Su[j]) = bootstrapCI(Xb, Yb, Smax[j], manager.delayLimit, manager.col, manager.bootCI, manager.fTransform, manager.zNormalize, manager.trendThresh)
    else:
      (Smax[j], Sl[j], Su[j]) = (Smax[j], Smax[j], Smax[j])

  data = np.column_stack((Smax, Sl, Su, Xs, Ys, Al, delay, lsaP, PCC, P_PCC, SPCC, P_SPCC, D_SPCC, SCC, P_SCC, SSCC, P_SSCC, D_SSCC, X, Y))
  return data

def palla_applyAnalysis(firstData, secondData, data, col, onDiag=True, delayLimit=3, minOccur=.5, bootCI=.95, bootNum=0, pvalueMethod='perm', 
                        precisionP=1000,fTransform=simpleAverage, zNormalize=noZeroNormalize, approxVar=1,resultFile=tempfile.TemporaryFile('w'), trendThresh=None,
                        firstFactorLabels=None, secondFactorLabels=None, qvalueMethod='R', progressive=0):

  warnings.filterwarnings("ignore", category=RuntimeWarning)
  global manager

  firstFactorNum = firstData.shape[0]
  firstRepNum = firstData.shape[1]
  firstSpotNum = firstData.shape[2]
  secondFactorNum = secondData.shape[0]
  secondRepNum = secondData.shape[1]
  secondSpotNum = secondData.shape[2]
  
  if not firstFactorLabels:
    firstFactorLabels= [str(v) for v in range(1, firstFactorNum+1)]
  if not secondFactorLabels:
    secondFactorLabels= [str(v) for v in range(1, secondFactorNum+1)]
  assert secondSpotNum == firstSpotNum 
  assert secondRepNum == firstRepNum

  pairwiseNum = firstFactorNum*secondFactorNum

  Smax = [0] * pairwiseNum
  Xs = [0] * pairwiseNum
  Ys = [0] * pairwiseNum
  Al = [0] * pairwiseNum
  delay = [0] * pairwiseNum

  x_series = []
  y_series = []
  x_occur = []
  y_occur = []

  timespots = secondSpotNum
  lengthSeries = secondSpotNum

  #这个就是所谓的LTA算法条件
  if trendThresh != None:
    lengthSeries = timespots - 1

  replicates = firstRepNum
  stdX = np.sqrt(approxVar)
  
  if qvalueMethod in ['R'] and rpy_import:
    qvalue_func = R_Qvalue
  else:
    qvalue_func = storeyQvalue
  
  P_table = 1
  if pvalueMethod in ['theo','mix','ddlsa','stlta']:
    P_table = theoPvalue(Rmax=lengthSeries, Dmax=delayLimit,precision=1./precisionP, x_decimal=my_decimal)

  x = []
  y = []
  res = []
  
  for i_0 in range(0,firstFactorNum):
      series0 = np.ma.masked_invalid(firstData[i_0], copy=True)
      series0_badOccur = np.sum(np.logical_not(np.isnan(ma_average(series0)),ma_average(series0)==0))/float(timespots) < minOccur

      x_occur.append(series0_badOccur)
      if series0.shape[1] == None:
        series0.shape = (1, series0.shape[0])
      x_series.append(series0)

      timespots = series0.shape[1]

      #如果是LTA计算，那么就要对数据进行另一种处理方式
      if trendThresh != None:
          xSeries = ji_calc_trend(zNormalize(fTransform(series0)),timespots,trendThresh).tolist()
      else:
          xSeries = zNormalize(fTransform(series0)).tolist()
      x.extend(xSeries)

  for i_1 in range(0,secondFactorNum):
      series1 = np.ma.masked_invalid(secondData[i_1], copy=True)
      series1_badOccur = np.sum(np.logical_not(np.isnan(ma_average(series1)),ma_average(series1)==0))/float(timespots) < minOccur
      y_occur.append(series1_badOccur)
      if series1.shape[1] == None:
        series1.shape = (1, series1.shape[0])
      y_series.append(series1)

      timespots = series1.shape[1]
      if trendThresh != None:
          ySeries = ji_calc_trend(zNormalize(fTransform(series1)),timespots,trendThresh).tolist()
      else:
          ySeries = zNormalize(fTransform(series1)).tolist()
      y.extend(ySeries)

  s = delayLimit
  data.assign(firstFactorNum, secondFactorNum, s, col, x, y)
  data.dp_lsa()

  res.append(list(data.score))
  res.append(list(data.x_0))
  res.append(list(data.x_1))
  res.append(list(data.y_0))
  res.append(list(data.y_1))

  Smax = res[0]
  Xs = res[1]
  Ys = res[3]
  Al = [1 + a - b for a, b in zip(res[2], res[1])]
  delay = [b - a for a, b in zip(res[3], res[1])]

  array = []

  with multiprocessing.Manager() as manager:
    num_processes = manager.num_processes = 2*os.cpu_count()
    manager.num_0 = firstFactorNum
    manager.num_1 = secondFactorNum
    manager.col = col
    manager.delayLimit = delayLimit
    manager.lengthSeries = lengthSeries
    manager.pvalueMethod = pvalueMethod
    manager.replicates = replicates
    manager.stdX = stdX
    manager.P_table = P_table
    manager.bootCI = bootCI
    manager.bootNum = bootNum
    manager.fTransform = fTransform
    manager.zNormalize = zNormalize
    manager.trendThresh = trendThresh
    manager.precisionP = precisionP
    manager.qvalue_func = qvalue_func

    manager.Smax = Smax
    manager.Xs = Xs
    manager.Ys = Ys
    manager.Al = Al
    manager.delay = delay

    manager.x_series = x_series
    manager.y_series = y_series
    manager.x_occur = x_occur
    manager.y_occur = y_occur

    pool = multiprocessing.Pool(processes=num_processes)
    results = [pool.apply_async(myfunc, args=(process_id,)) for process_id in range(firstFactorNum)]

    data_segments = []
    for result in results:
      data_segment = result.get()
      data_segments.append(data_segment)
    array = np.vstack(data_segments)
    pool.close()
    pool.join()

  data.lsa_clean()
  return array


if __name__=="__main__":
  print("hello world!")
