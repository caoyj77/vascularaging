import numpy as np
import pandas as pd
import xgboost as xgb
from sklearn.model_selection import KFold
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score # type: ignore

# 模拟生成示例数据（实际需替换）
dpath = '/Users/'
train_tgt_out_df = pd.read_csv(dpath+'proteins.csv') 
#significant_markers = train_tgt_out_df.loc[train_tgt_out_df['pval_bfi'] < 0.05, 'protein_code'].tolist()
significant_markers = train_tgt_out_df.loc[train_tgt_out_df['adj_p'] < 0.05, 'Protein'].tolist()

#train_df = pd.read_csv(dpath+'mydf29.csv') 
#train_df = train_df[train_df['Sex'] == 1]
#train_df = train_df.reset_index(drop=True)
cov = pd.read_csv(dpath+'discovery.csv') 
pd = pd.read_csv(dpath+'PD.csv') 
train_df = pd.merge(cov, pd, on='eid', how='inner')

X = train_df[significant_markers]# 特征数据
y = train_df['Age']  # 年龄标签

# 创建KFold对象
kf = KFold(n_splits=5, shuffle=True, random_state=42)

# 存储每次折叠的均方误差（MSE）和决定系数（R²）
mae_scores = []
r2_scores = []

# 五折交叉验证循环
for train_index, test_index in kf.split(X):
    X_train, X_test = X.iloc[train_index], X.iloc[test_index]
    y_train, y_test = y[train_index], y[test_index]
    #y_test.to_csv(dpath + 'predictxb_test.csv',mode='a',header=False, index=False)
    # 将数据转换为XGBoost所需的DMatrix格式（针对回归任务）
    dtrain = xgb.DMatrix(X_train, label=y_train)
    dtest = xgb.DMatrix(X_test)
    
    # 设置XGBoost回归器参数（可按需调整参数）
    params = {
        'objective':'reg:squarederror',
        'eval_metric': 'rmse'
    }
    
    # 训练模型
    model = xgb.train(params, dtrain)
    
    # 在测试集上进行预测
    y_pred = model.predict(dtest)
    df = pd.DataFrame(y_pred)
    # 计算均方误差和决定系数
    mae = mean_absolute_error(y_test, y_pred)
    r2 = r2_score(y_test, y_pred)
    mae_scores.append(mae)
    r2_scores.append(r2)
    print(f"Fold {len(mae_scores)} MSE: {mae}, R²: {r2}")
    
#X = xgb.DMatrix(X)
#pre=model.predict(X)
#pre.to_csv(dpath + 'km/prediction_xg_d.csv',mode='a',header=False, index=False)

# 计算平均均方误差和平均决定系数
#average_mse = np.mean(mse_scores)
#print("平均均方误差（MAE）:", average_mse)
#print("平均决定系数（R²）:", r2_scores)

table2 = pd.read_csv(dpath+'mydf_r.csv') 

X_pre=table2[significant_markers]
X_pre = xgb.DMatrix(X_pre)
y_pre=model.predict(X_pre) 
y_pre = pd.DataFrame(y_pre) 
eid_col = table2['eid']  # 标签列
y_true = table2['Age']   # 真实年龄 
y_pre = y_pre.values.ravel() 
age_gap=y_pre-y_true
output_df = pd.DataFrame({
    'eid': eid_col,           # 假设table2中存在'eid'列
    'true_age': y_true,      # 真实年龄
    'predicted_age': y_pre,          # 预测年龄
    'age_gap':age_gap
}) 
print(output_df)
output_df.to_csv(dpath + 'prediction.csv',mode='a',header=True, index=False)
y_true=table2['Age']
#mae = mean_absolute_error(y_true, y_pre)
#r2 = r2_score(y_true, y_pre)
#mae_scores.append(mae)
#r2_scores.append(r2)
#print(f"Fold {len(mae_scores)} MAE:{mae}, R²: {r2}")
