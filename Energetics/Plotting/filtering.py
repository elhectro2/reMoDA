import pandas as pd
import statsmodels.api as sm
from scipy.signal import lfilter, savgol_filter
from tsmoothie.smoother import *


def get_y_lfilter(y, n=15, a=1):
    b = [1.0 / n] * n
    return lfilter(b, a, y)


def get_y_sav_gol(y, win_len=101, polyorder=2, mode="nearest"):
    return savgol_filter(y, win_len, polyorder, mode=mode)


def get_y_convolution(y, win_len=30, win_type="ones"):
    smoother = ConvolutionSmoother(window_len=win_len, window_type=win_type)
    smoother.smooth(y)
    return smoother.smooth_data[0]


def get_y_lowess(y, x, fraction=0.3):
    return sm.nonparametric.lowess(y, x, frac=fraction)[:, 1], sm.nonparametric.lowess(y, x, frac=fraction)[:, 0]


def get_y_rolling_mean(y, x, win_len=30):
    df = pd.DataFrame(y, x)
    return df.rolling(win_len).mean()[0]
