import glob
import numpy as np
import pandas as pd
import statsmodels.api as sm
from tqdm import tqdm
import os
from sklearn.impute import SimpleImputer
import re
from statsmodels.miscmodels.ordinal_model import OrderedModel
from statsmodels.stats.multitest import fdrcorrection
from mne.stats import bonferroni_correction
from joblib import Parallel, delayed
from sklearn.linear_model import Lasso
import warnings
warnings.filterwarnings('error')
from sklearn.linear_model import LinearRegression


def results_summary(tgt_out_df):
    oratio_out_lst, p_out_lst = [], []
    for i in range(len(tgt_out_df)):
        oratio = f'{tgt_out_df.oratio.iloc[i]:.2f}'
        lbd = f'{tgt_out_df.or_lbd.iloc[i]:.2f}'
        ubd = f'{tgt_out_df.or_ubd.iloc[i]:.2f}'
        oratio_out_lst.append(oratio + ' [' + lbd + '-' + ubd + ']')
        if tgt_out_df.pval_bfi.iloc[i] < 0.001:
            p_out_lst.append('***')
        elif tgt_out_df.pval_bfi.iloc[i] < 0.01:
            p_out_lst.append('**')
        elif tgt_out_df.pval_bfi.iloc[i] < 0.05:
            p_out_lst.append('*')
        else:
            p_out_lst.append('')
    return (oratio_out_lst, p_out_lst)


def process(nmr_f, tmp_df, cov_f_lst):
    tmp_df.rename(columns={nmr_f: 'x_nmr'}, inplace=True)
    rm_eid_idx = tmp_df.index[tmp_df.x_nmr.isnull() == True]
    tmp_df.drop(rm_eid_idx, axis=0, inplace=True)
    tmp_df.reset_index(inplace=True, drop=True)
    nb_all, nb_case = len(tmp_df), tmp_df.targe.sum()
    prop_case = np.round(nb_case / nb_all * 100, 3)
    Y = tmp_df.targe
    X = tmp_df[cov_f_lst + ['x_nmr']]
    # 创建填充器
    imputer = SimpleImputer(strategy='median')  # 或 'median'、'constant' 等
    # 填充 X_train 的缺失值
    x = pd.DataFrame(imputer.fit_transform(X), columns=X.columns)

    
    try:
        lasso=Lasso(alpha=0.1)
        log_mod = lasso.fit(sm.add_constant(x),Y)
        beta_val=log_mod.coef_
        nine_beta=beta_val[8]
        model_ols = sm.OLS(Y, sm.add_constant(x)).fit(alpha=lasso.alpha, L1_wt=1)
        model_ols = sm.OLS(Y, sm.add_constant(x)).fit()
        #print(model_ols.summary())
        oratio = np.round(np.exp(model_ols.params).loc['x_nmr'], 5)
        pval = model_ols.pvalues.loc['x_nmr']
        ci_mod = model_ols.conf_int(alpha=0.05)
        lbd, ubd = np.round(np.exp(ci_mod.loc['x_nmr'][0]), 5), np.round(np.exp(ci_mod.loc['x_nmr'][1]), 5)
        tmpout = [nmr_f, nb_all, nb_case, prop_case, oratio, lbd, ubd, pval, nine_beta]
        #tmpout = [nmr_f, nb_all, nb_case, prop_case, oratio, lbd, ubd, pval]
    except KeyError as e:
        print(f"出现键错误：{e}，可能 'x_nmr' 不在模型参数中。")
        tmpout = [nmr_f, nb_all, nb_case, prop_case, np.nan, np.nan, np.nan, np.nan, np.nan]
    except Exception as e:
        print(f"出现未知错误：{e}")
        tmpout = [nmr_f, nb_all, nb_case, prop_case, np.nan, np.nan, np.nan, np.nan, np.nan]
    return tmpout

def sort_nicely(l):
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c.replace("_","")) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )
    return l

nb_cpus =10

dpath = '/Users/'
badoutfile = dpath + 'bad_tgt_all.csv'

tgt_df = pd.read_csv(dpath + 'name.csv')
tgt_lst = tgt_df.NAME.tolist()

nmr_df = pd.read_csv(dpath + 'PD.csv')
nmr_f_lst = nmr_df.columns.tolist()[2:]

cov_df = pd.read_csv(dpath + 'cov_de.csv')
mydf = pd.read_csv(dpath + 'mydf.csv')
#mydf = pd.read_csv(dpath + 'mydf_r.csv')

#cov_f_lst = ["Age", "Sex", "TDI", "BMI", "Smoke", "FastingTime", "SBP", "BL_HYPT"]
cov_f_lst = ["Age", "TDI", "BMI", "Smoke", "FastingTime", "SBP", "BL_HYPT"]

bad_tgt = []

for tgt_name in tqdm(tgt_lst):
    tgt_file = dpath+tgt_name+'.csv'
    tmp_tgt_df = pd.read_csv(tgt_file, usecols=['eid', 'targe', 'yrs'] )
    rm_bl_idx = tmp_tgt_df.index[(tmp_tgt_df.targe == 1) & (tmp_tgt_df.yrs > 0)]
    tmp_tgt_df.drop(rm_bl_idx, axis=0, inplace=True)
    tmp_tgt_df.reset_index(inplace=True, drop=True)
    tmp_df = pd.merge(mydf, tmp_tgt_df, how='inner', on=['eid'])
    #tmp_df = tmp_df[tmp_df['Sex'] == 0]
    #tmp_df = tmp_df.reset_index(drop=True)
    
    try:
        if tmp_df.targe.sum()>=50:
            tgt_out_df = Parallel(n_jobs=nb_cpus)(delayed(process)(nmr_f, tmp_df, cov_f_lst) for nmr_f in nmr_f_lst)
            tgt_out_df = pd.DataFrame(tgt_out_df)
            tgt_out_df.columns = ['protein_code','nb_individuals', 'nb_case', 'prop_case(%)', 'oratio', 'or_lbd', 'or_ubd', 'pval_raw','beta']
            _, p_f_bfi = bonferroni_correction(tgt_out_df.pval_raw.fillna(1), alpha=0.05)
            tgt_out_df['pval_bfi'] = p_f_bfi
            tgt_out_df.loc[tgt_out_df['pval_bfi'] >= 1, 'pval_bfi'] = 1
            tgt_out_df['or_output'], tgt_out_df['pval_significant'] = results_summary(tgt_out_df)
            tgt_out_df = tgt_out_df[[ 'protein_code','nb_individuals', 'nb_case', 'prop_case(%)', 'oratio',
                                     'or_lbd', 'or_ubd', 'pval_raw', 'pval_bfi', 'or_output', 'pval_significant']]
            tgt_out_df.to_csv(dpath+tgt_name+'.csv', index=False)
            print(tgt_out_df)
        else:
            bad_tgt.append(tgt_name)
    except:
        bad_tgt.append(tgt_name)


bad_df = pd.DataFrame(bad_tgt)
bad_df.to_csv(badoutfile, index=False)
